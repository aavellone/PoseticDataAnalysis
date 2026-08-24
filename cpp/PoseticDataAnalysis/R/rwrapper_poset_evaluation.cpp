/**
 * @file rwrapper_poset_evaluation.cpp
 * @brief R/C++ interface — POSet evaluation functions (HPC Optimized).
 *
 * Implements all R-callable (`extern "C"`) entry points for:
 * - Exact MRP (Mutual Ranking Probability) computation via tree-of-ideals
 * enumeration.
 * - Approximate MRP via the Bubley-Dyer Markov-chain Monte Carlo sampler.
 * - Generic exact evaluation of user-supplied functions over the linear
 * extensions.
 * - Brügemann-Lerche-Sørensen dominance matrix.
 *
 * Every public function follows the zero-copy C++20 paradigm:
 * 1. Unwrap R arguments into stack-allocated or unique_ptr structures.
 * 2. Execute the C++ computation leveraging move semantics.
 * 3. Return an R named list (or matrix) via `RCreate` helpers.
 *
 * @author Alessandro Avellone
 * @version 3.0 (HPC C++20)
 * @date 2025
 */

#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>

#include "dimensionality_reduction.h"
#include "display_message.h"
#include "function_linear_extension.h"
#include "function_linear_extension_average_height.h"
#include "function_linear_extension_mutual_ranking_probability.h"
#include "function_linear_extension_separation_asymmetric_lower.h"
#include "function_linear_extension_separation_asymmetric_upper.h"
#include "function_linear_extension_separation_symmetric.h"
#include "generic_functions.h"
#include "linear_extension_generator.h"
#include "linear_extension_generator_tree_of_ideals.h"
#include "linear_extension_generator_bubley_dyer.h"
#include "linear_extension_generator_binary_variable.h"
#include "linear_generator_wrapper.h"
#include "tensor.h"
#include "my_exception.h"
#include "poset.h"
#include "poset_wrapper.h"
#include "r_display.h"
#include "r_function_linear_extension.h"
#include "r_separation.h"
#include "random.h"
#include "rwrapper.h"
#include "rwrapper_conversion.h"
#include "separation.h"

/**
 * @brief Runs exact linear-extension enumeration and accumulates results.
 *
 * Creates one @ref FunctionLinearExtension per entry in @p internal_functions,
 * allocates the corresponding result matrix, then calls
 * `POSet::evaluation()` which iterates over every linear extension exactly
 * once using a Tree-of-Ideals generator.
 *
 * @param[in]  poset               The POSet raw pointer (observer).
 * @param[in]  internal_functions  Internal C++ function identifiers.
 * @param[in]  external_functions  Matching R function SEXPs.
 * @param[in]  output_interval     Progress reporting interval in seconds (max()
 * = silent).
 * @param[out] eval_results        Vector of unique_ptrs to result matrices.
 * @param[out] fles                Vector of unique_ptrs to FLE objects.
 * @param[out] functions_name      Display name for each function.
 * @param[out] le_count            Total number of linear extensions enumerated.
 *
 * @throws MyException if an unknown function name is encountered.
 */
void BuildExactEvaluation(POSet *poset,
    const std::vector<std::string> &internal_functions,
    const std::vector<SEXP> &external_functions,
                          std::optional<std::uint64_t> output_interval,
    std::vector<std::unique_ptr<Tensor<double, 2>>> &eval_results,
    std::vector<std::unique_ptr<FunctionLinearExtension>> &fles,
    std::vector<std::string> &functions_name, std::uint64_t &le_count) {

    const std::size_t num_funcs = internal_functions.size();
    eval_results.resize(num_funcs);
    fles.resize(num_funcs);
    functions_name.resize(num_funcs);

    for (std::uint64_t k = 0; k < num_funcs; ++k) {
        const auto &name = internal_functions[k];
        SEXP r_func = external_functions[k];

        auto it = POSetWrap::kFunctionLinearMapType.find(name);
        if (it == POSetWrap::kFunctionLinearMapType.end()) {
            throw MyException(
                "BuildExactEvaluation: unknown function '" + name + "'");
        }

        std::unique_ptr<FunctionLinearExtension> fle = nullptr;
        std::string display_name = "";

        switch (it->second) {
            case POSetWrap::FunctionLinearType::kMutualRankingProbability:
                fle = std::make_unique<FLEMutualRankingProbability>(poset);
                display_name = "MutualRankingProbability";
                break;
            case POSetWrap::FunctionLinearType::kAverageHeight:
                fle = std::make_unique<FLEAverageHeight>(poset);
                display_name = "AverageHeight";
                break;
            case POSetWrap::FunctionLinearType::kSeparationAsymmetricLower:
                fle = std::make_unique<FLESeparationAsymmetricLower>(poset);
                display_name = "asymmetricLower";
                break;
            case POSetWrap::FunctionLinearType::kSeparationAsymmetricUpper:
                fle = std::make_unique<FLESeparationAsymmetricUpper>(poset);
                display_name = "asymmetricUpper";
                break;
            case POSetWrap::FunctionLinearType::kSeparationSymmetric:
                fle = std::make_unique<FLESeparationSymmetric>(poset);
                display_name = "symmetric";
                break;
            case POSetWrap::FunctionLinearType::kRFunction:
                fle = std::make_unique<FLERInterface>(poset, r_func);
                display_name = "RFunction";
                break;
            case POSetWrap::FunctionLinearType::kDominance:
            case POSetWrap::FunctionLinearType::kMannWhitneyDominance:
            case POSetWrap::FunctionLinearType::kMannWhitneyInferentialDominance:
                throw MyException("kDominance, kMannWhitneyDominance, kMannWhitneyInferentialDominance: not allowed here");
                break;
        }

        const auto shape = fle->Shape();
        eval_results[k] =
        std::make_unique<Tensor<double, 2>>(std::array<std::uint64_t, 2>{shape.at(0), shape.at(1)}, 0.0);
        fles[k] = std::move(fle);
        functions_name[k] = display_name;
    }

    // --- Set up the Tree-of-Ideals generator using stack vector ---
    auto lattice_of_ideals = poset->GetLatticeOfIdeals();
    LEGTreeOfIdeals le_generator(poset->size(), *lattice_of_ideals);
    le_generator.Start(0);

    // --- Run enumeration ---
    le_count = 0;
    bool end_process = false;
    std::unique_ptr<DisplayMessage> display;
    if (!output_interval.has_value()) {
        display = std::make_unique<DisplayMessageNull>();
    } else {
        display = std::make_unique<DisplayMessageEvaluationExactR>(le_count, output_interval.value());
    }
    POSet::evaluation(fles, le_generator, eval_results, le_count, end_process,
        display.get(), POSet::EvaluationUpdateStrategy::Average);
}




// ===========================================================================
// R-callable functions
// ===========================================================================

extern "C" {
    /**
     * @brief Computes the exact Mutual Ranking Probability (MRP) matrix.
     *
     * @param poset_r               ExternalPtr to the wrapped POSet object.
     * @param output_ogni_in_sec_r  Numeric defining the progress print interval in
     * seconds.
     * @return SEXP A named list containing 'mrp' (matrix) and 'n' (integer count).
     */
    SEXP ExactMRP(SEXP poset_r, SEXP output_ogni_in_sec_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;

        try {
            const POSetWrap *poset_wrap = RConvert::ToPOSetWrap(poset_r);
            POSet *poset = poset_wrap->GetPOSet();

            std::optional<std::uint64_t> output_ogni = RConvert::ToOptionalUInt(output_ogni_in_sec_r);

            std::unique_ptr<LinearExtensionGenerator> legenerator = nullptr;
            if (poset_wrap->GetType() == POSetWrap::PosetType::kBinaryVariable) {
                legenerator = std::make_unique<LEGBinaryVariable>(poset->size());
            } else {
                auto lattice_of_ideals = poset->GetLatticeOfIdeals();
                legenerator = std::make_unique<LEGTreeOfIdeals>(
                    poset->size(), *lattice_of_ideals);
            }

            legenerator->Start(0);

            // Allocate MRP matrix and run evaluation
            const std::uint64_t n = poset->size();
            std::vector<std::unique_ptr<Tensor<double, 2>>> eval_results;
            eval_results.push_back(std::make_unique<Tensor<double, 2>>(std::array<std::uint64_t, 2>{n, n}, 0.0));

            std::vector<std::unique_ptr<FunctionLinearExtension>> fles;
            fles.push_back(std::make_unique<FLEMutualRankingProbability>(poset));

            std::uint64_t le_count = 0;
            bool end_process = false;
            std::unique_ptr<DisplayMessage> display;
            if (!output_ogni.has_value()) {
                display = std::make_unique<DisplayMessageNull>();
            } else {
                display = std::make_unique<DisplayMessageEvaluationExactR>(le_count, output_ogni.value());
            }

            POSet::evaluation(fles, *legenerator, eval_results, le_count,
                end_process, display.get(),
                POSet::EvaluationUpdateStrategy::Average);

            // Build R result
            const int ni = static_cast<int>(n);
            SEXP mrp_r = guard.Protect(Rf_allocMatrix(REALSXP, ni, ni));
            auto &mrp = eval_results[0];
            for (int row = 0; row < ni; ++row) {
                for (int col = 0; col < ni; ++col) {
                    REAL(mrp_r)
                    [row + ni * col] = (*mrp)(static_cast<std::uint64_t>(row),
                        static_cast<std::uint64_t>(col));
                }
            }

            SEXP dimnames = guard.Protect(Rf_allocVector(VECSXP, 2));
            SEXP r_names = guard.Protect(Rf_allocVector(STRSXP, ni));
            for (int i = 0; i < ni; ++i) {
                std::string_view name_view = poset->GetElementName(static_cast<std::size_t>(i));
                
                SET_STRING_ELT(
                               r_names,
                               static_cast<R_xlen_t>(i),
                               Rf_mkCharLenCE(name_view.data(), static_cast<int>(name_view.size()), CE_UTF8)
                               );
            }
            
            SET_VECTOR_ELT(dimnames, 0, r_names);
            SET_VECTOR_ELT(dimnames, 1, r_names);
            Rf_setAttrib(mrp_r, R_DimNamesSymbol, dimnames);

            SEXP res = guard.Protect(Rf_allocVector(VECSXP, 2));
            SEXP res_names = guard.Protect(Rf_allocVector(STRSXP, 2));

            SET_VECTOR_ELT(res, 0, mrp_r);
            SET_STRING_ELT(res_names, 0, Rf_mkChar("mrp"));

            SET_VECTOR_ELT(
                res, 1, RCreate::FromInt(guard, static_cast<int>(le_count)));
            SET_STRING_ELT(res_names, 1, Rf_mkChar("n"));

            Rf_setAttrib(res, R_NamesSymbol, res_names);
            result_r = res;
        } catch (std::exception &ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }

        return result_r;
    }

    /**
     * @brief Creates and initializes a Bubley-Dyer MRP MCMC generator object.
     *
     * @param poset_r ExternalPtr to the wrapped POSet object.
     * @param seed_r  Numeric random seed for MCMC reproducibility (optional).
     * @return SEXP ExternalPtr mapping the `BubleyDyerMRPGenerator` instance.
     */
    SEXP BuildBubleyDyerMRPGenerator(SEXP poset_r, SEXP seed_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;

        try {
            // 1. Estrazione del raw pointer
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);
            
            // 2. Estrazione del seed
            std::optional<std::uint64_t> seed = RConvert::ToOptionalSeed(seed_r);
            if (!seed.has_value()) {
                seed = Random::GLOBAL.RndNextInt(0, std::numeric_limits<std::uint64_t>::max());
            }
            
            // 3. Costruzione dell'oggetto C++
            auto generator = BubleyDyerMRPGenerator::BuildBubleyDyerMRPGenerator(poset_wrap, seed.value());

            // 4. Creazione del puntatore esterno in R
            SEXP ptr_r = RCreate::WrapExternalPtr(guard, std::move(generator));

            // 5. Creazione dello "scudo" per il Garbage Collector
            R_SetExternalPtrProtected(ptr_r, poset_r);

            // 6. Il seme usato viene restituito a R come stringa: reimmetterlo
            // riproduce esattamente la stessa sessione.
            result_r = RCreate::GeneratorWithSeed(guard, ptr_r, seed.value());
        } catch (std::exception &ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }

        return result_r;
    }

    /**
     * @brief Runs (or continues) the Bubley-Dyer MCMC sampler and returns the
     * approximated MRP matrix.
     *
     * @param bubley_dyer_r        ExternalPtr to the active
     * `BubleyDyerMRPGenerator`.
     * @param quante_r             Number of MCMC steps to execute (numeric).
     * @param errore_r             Target error bound threshold (numeric).
     * @param output_ogni_in_sec_r Print interval for R console feedback.
     * @return SEXP A named list containing the approximate 'mrp' matrix and 'n'
     * steps performed.
     */
    SEXP BubleyDyerMRP(SEXP bubley_dyer_r, SEXP quante_r, SEXP errore_r,
        SEXP output_ogni_in_sec_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;

        try {

            auto bd = RConvert::ToWrapPtr<BubleyDyerMRPGenerator>(
                bubley_dyer_r, "Il generatore MRP Bubley-Dyer");

            std::optional<std::uint64_t> quante = RConvert::ToOptionalUInt(quante_r);
            std::optional<double> errore = RConvert::ToOptionalDouble(errore_r);
            std::optional<std::uint64_t> output_ogni = RConvert::ToOptionalUInt(output_ogni_in_sec_r);

            if (!errore.has_value() && !quante.has_value()) {
                throw MyException("BubleyDyerMRP: provide either le extensions "
                                  "count or a target error.");
            }
            if (errore.has_value() && quante.has_value()) {
                forward_warning_to_r(
                    "BubleyDyerMRP: 'error' is ignored when 'n' is provided.");
            }

            POSet *poset = bd->poset_;
            const int n = static_cast<int>(poset->size());

            // Start or extend the Markov chain
            bool updated = true;
            if (!bd->used_) {
                if (!quante.has_value()) {
                    quante = bd->le_generator_->EvaluateNumberOfIteration(errore.value());
                }
                bd->le_generator_->Start(quante.value());
            } else {
                updated = bd->le_generator_->UpdateCounters(quante, errore);
                if (!updated)
                    Rprintf("BubleyDyerMRP: desired error already reached — no new "
                            "extensions.\n");
                else
                    bd->le_generator_->Next();
            }

            std::uint64_t total_le = bd->le_generator_->CurrentNumberOfLe() - 1;

            if (updated) {
                bd->used_ = true;

                std::vector<std::unique_ptr<Tensor<double, 2>>> eval_results;
                eval_results.push_back(std::move(bd->mrp_));

                std::vector<std::unique_ptr<FunctionLinearExtension>> fles;
                fles.push_back(std::make_unique<FLEMutualRankingProbability>(poset));

                std::uint64_t le_count = 0;
                bool end_process = false;

                const std::uint64_t total_ext =
                    (bd->le_generator_->NumberOfLe() + 1) -
                    bd->le_generator_->CurrentNumberOfLe() + 1;

                std::unique_ptr<DisplayMessage> display;
                if (!output_ogni.has_value()) {
                    display = std::make_unique<DisplayMessageNull>();
                } else {
                    display = std::make_unique<DisplayMessageEvaluationBDR>(
                        le_count, total_ext, output_ogni.value());
                }

                POSet::evaluation(fles, *(bd->le_generator_), eval_results,
                    le_count, end_process, display.get(),
                    POSet::EvaluationUpdateStrategy::Average);

                bd->mrp_ = std::move(eval_results[0]);
                total_le += le_count;
            }

            // Build R result matrix

            SEXP matrix_r = guard.Protect(Rf_allocMatrix(REALSXP, n, n));
            for (int row = 0; row < n; ++row) {
                for (int col = 0; col < n; ++col) {
                    REAL(matrix_r)
                    [row + n * col] = (*(bd->mrp_))(static_cast<std::uint64_t>(row),
                        static_cast<std::uint64_t>(col));
                }
            }

            SEXP dimnames = guard.Protect(Rf_allocVector(VECSXP, 2));
            SEXP r_names = guard.Protect(Rf_allocVector(STRSXP, n));
            for (int i = 0; i < n; ++i) {
                std::string_view name_view = poset->GetElementName(static_cast<std::size_t>(i));
                SET_STRING_ELT(
                               r_names,
                               static_cast<R_xlen_t>(i),
                               Rf_mkCharLenCE(name_view.data(), static_cast<int>(name_view.size()), CE_UTF8)
                               );
            }
            
            SET_VECTOR_ELT(dimnames, 0, r_names);
            SET_VECTOR_ELT(dimnames, 1, r_names);
            Rf_setAttrib(matrix_r, R_DimNamesSymbol, dimnames);

            SEXP res = guard.Protect(Rf_allocVector(VECSXP, 2));
            SEXP res_names = guard.Protect(Rf_allocVector(STRSXP, 2));

            SET_VECTOR_ELT(res, 0, matrix_r);
            SET_STRING_ELT(res_names, 0, Rf_mkChar("mrp"));

            SET_VECTOR_ELT(
                res, 1, RCreate::FromInt(guard, static_cast<int>(total_le)));
            SET_STRING_ELT(res_names, 1, Rf_mkChar("n"));

            Rf_setAttrib(res, R_NamesSymbol, res_names);
            result_r = res;
        } catch (std::exception &ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }

        return result_r;
    }

    /**
     * @brief Creates a Bubley-Dyer evaluation generator for multiple functions.
     *
     * @param poset_r     ExternalPtr to the wrapped POSet object.
     * @param seed_r      Numeric random seed.
     * @param functions_r R List containing the mapped functions to execute.
     * @return SEXP ExternalPtr mapping the `BubleyDyerEvaluationGenerator`
     * instance.
     */
    SEXP BuildBubleyDyerEvaluationGenerator(
        SEXP poset_r, SEXP seed_r, SEXP functions_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;

        try {
            // 1. Estrazione del raw pointer
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);

            // 2. Estrazione del seed
            std::optional<std::uint64_t> seed = RConvert::ToOptionalSeed(seed_r);
            if (!seed.has_value()) {
                seed = Random::GLOBAL.RndNextInt(0, std::numeric_limits<std::uint64_t>::max());
            }
            
            // 3. Estrazione delle funzioni seed
            std::vector<std::string> internal_funcs;
            std::vector<SEXP> external_funcs;
            RConvert::ToFunctions(functions_r, internal_funcs, external_funcs);
            
            // 4. Costruzione dell'oggetto C++
            auto generator =
            BubleyDyerEvaluationGenerator::BuildBubleyDyerEvaluationGenerator(
                                                                              poset_wrap, seed.value(), internal_funcs, external_funcs);
            // 5. Creazione del puntatore esterno in R
            SEXP ptr_r = RCreate::WrapExternalPtr(guard, std::move(generator));

            // 6. Creazione dello "scudo" per il Garbage Collector
            R_SetExternalPtrProtected(ptr_r, poset_r);

            // 7. Il seme usato viene restituito a R come stringa: reimmetterlo
            // riproduce esattamente la stessa sessione.
            result_r = RCreate::GeneratorWithSeed(guard, ptr_r, seed.value());

        } catch (std::exception &ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }

        return result_r;
    }

    /**
     * @brief Runs (or continues) a multi-function Bubley-Dyer evaluation loop.
     *
     * @param generator_r          ExternalPtr to the active
     * `BubleyDyerEvaluationGenerator`.
     * @param quante_r             Number of MCMC iterations to run.
     * @param errore_r             Target error bound threshold.
     * @param output_ogni_in_sec_r Print interval for R console feedback.
     * @return SEXP A named list where each element is the resulting matrix for a
     * requested function.
     */
    SEXP BubleyDyerEvaluation(
        SEXP generator_r, SEXP quante_r, SEXP errore_r, SEXP output_ogni_in_sec_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;

        try {
            // Punto d'ingresso condiviso con BubleyDyerSeparation (che vi delega):
            // il nome nel messaggio d'errore resta percio' generico.
            auto gen = RConvert::ToWrapPtr<BubleyDyerEvaluationGenerator>(
                generator_r, "Il generatore");

            std::optional<std::uint64_t> quante = RConvert::ToOptionalUInt(quante_r);
            std::optional<double> errore = RConvert::ToOptionalDouble(errore_r);
            std::optional<std::uint64_t> output_ogni = RConvert::ToOptionalUInt(output_ogni_in_sec_r);

            if (!errore.has_value() && !quante.has_value()) {
                throw MyException("BubleyDyerEvaluation: provide either le "
                                  "extensions count or a target error.");
            }
            if (errore.has_value() && quante.has_value()) {
                forward_warning_to_r("BubleyDyerEvaluation: 'error' is ignored "
                                     "when 'n' is provided.");
            }

            bool updated = true;
            if (!gen->used_) {
                if (!quante.has_value()) {
                    quante = gen->le_generator_->EvaluateNumberOfIteration(errore.value());
                }
                gen->le_generator_->Start(quante.value());
            } else {
                updated = gen->le_generator_->UpdateCounters(quante, errore);
                if (!updated) {
                    Rprintf("BubleyDyerEvaluation: desired error already reached — "
                            "no new extensions.\n");
                } else {
                    gen->le_generator_->Next();
                }
            }

            std::uint64_t total_le = gen->le_generator_->CurrentNumberOfLe() - 1;
            if (updated) {
                gen->used_ = true;
                std::uint64_t le_count = 0;
                bool end_process = false;
                const std::uint64_t total_ext =
                    (gen->le_generator_->NumberOfLe() + 1) -
                    gen->le_generator_->CurrentNumberOfLe();

                std::unique_ptr<DisplayMessage> display;
                if (!output_ogni.has_value()) {
                    display = std::make_unique<DisplayMessageNull>();
                } else {
                    display = std::make_unique<DisplayMessageEvaluationBDR>(
                        le_count, total_ext, output_ogni.value());
                }

                POSet::evaluation(gen->fles_, *(gen->le_generator_),
                    gen->eval_results_, le_count, end_process, display.get(),
                    POSet::EvaluationUpdateStrategy::Average);
                total_le += le_count;
            }

            const int n_matrix = static_cast<int>(gen->eval_results_.size());
            SEXP res = guard.Protect(Rf_allocVector(VECSXP, n_matrix + 1));
            SEXP res_names = guard.Protect(Rf_allocVector(STRSXP, n_matrix + 1));

            for (int k = 0; k < n_matrix; ++k) {
                const auto &mat = gen->eval_results_[static_cast<std::size_t>(k)];
                const int nrow = static_cast<int>(mat->Extent(0));
                const int ncol = static_cast<int>(mat->Extent(1));

                SEXP r_mat = guard.Protect(Rf_allocMatrix(REALSXP, nrow, ncol));
                for (int r = 0; r < nrow; ++r) {
                    for (int c = 0; c < ncol; ++c) {
                        REAL(r_mat)
                        [r + nrow * c] = (*mat)(static_cast<std::uint64_t>(r),
                            static_cast<std::uint64_t>(c));
                    }
                }
                // Add dimnames
                SEXP dimnames = guard.Protect(Rf_allocVector(VECSXP, 2));

                // Gestione nomi delle righe
                SEXP r_rn = guard.Protect(Rf_allocVector(STRSXP, nrow));
                for (int i = 0; i < nrow; ++i) {
                    std::string_view row_view = gen->fles_[static_cast<std::size_t>(k)]->GetRowNameAt(i);
                    SET_STRING_ELT(r_rn, i,
                                   Rf_mkCharLenCE(row_view.data(), static_cast<int>(row_view.size()), CE_UTF8));
                }
                
                // Gestione nomi delle colonne
                SEXP r_cn = guard.Protect(Rf_allocVector(STRSXP, ncol));
                for (int i = 0; i < ncol; ++i) {
                    std::string_view col_view = gen->fles_[static_cast<std::size_t>(k)]->GetColNameAt(i);
                    
                    SET_STRING_ELT(r_cn, i,
                                   Rf_mkCharLenCE(col_view.data(), static_cast<int>(col_view.size()), CE_UTF8));
                }

                SET_VECTOR_ELT(dimnames, 0, r_rn);
                SET_VECTOR_ELT(dimnames, 1, r_cn);
                Rf_setAttrib(r_mat, R_DimNamesSymbol, dimnames);

                SET_VECTOR_ELT(res, k, r_mat);
                SET_STRING_ELT(res_names, k,
                    Rf_mkChar(
                        gen->functions_name_[static_cast<std::size_t>(k)].c_str()));
            }

            SET_VECTOR_ELT(
                res, n_matrix, RCreate::FromInt(guard, static_cast<int>(total_le)));
            SET_STRING_ELT(res_names, n_matrix, Rf_mkChar("n"));
            Rf_setAttrib(res, R_NamesSymbol, res_names);

            result_r = res;
        } catch (std::exception &ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        return result_r;
    }

    /**
     * @brief Computes the exact evaluation of multiple functions over all linear
     * extensions.
     *
     * @param poset_r              ExternalPtr to the wrapped POSet object.
     * @param output_ogni_in_sec_r Print interval for R console feedback.
     * @param functions_r          R List containing the functions to be mapped and
     * evaluated.
     * @return SEXP A named list of resulting exact matrices.
     */
    SEXP ExactEvaluation(
        SEXP poset_r, SEXP output_ogni_in_sec_r, SEXP functions_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;

        try {
            const POSetWrap *poset_wrap = RConvert::ToPOSetWrap(poset_r);
            std::optional<std::uint64_t> output_ogni = RConvert::ToOptionalUInt(output_ogni_in_sec_r);

            std::vector<std::string> internal_funcs;
            std::vector<SEXP> external_funcs;
            RConvert::ToFunctions(functions_r, internal_funcs, external_funcs);

            const std::size_t nf = internal_funcs.size();
            std::vector<std::unique_ptr<Tensor<double, 2>>> eval_results;
            std::vector<std::unique_ptr<FunctionLinearExtension>> fles;
            std::vector<std::string> func_names;
            std::uint64_t le_count = 0;

            BuildExactEvaluation(poset_wrap->GetPOSet(), internal_funcs,
                external_funcs, output_ogni, eval_results, fles, func_names,
                le_count);

            const int n_matrix = static_cast<int>(nf);
            SEXP res = guard.Protect(Rf_allocVector(VECSXP, n_matrix + 1));
            SEXP res_names = guard.Protect(Rf_allocVector(STRSXP, n_matrix + 1));

            for (int k = 0; k < n_matrix; ++k) {
                const auto &mat = eval_results[static_cast<std::size_t>(k)];
                const int nrow = static_cast<int>(mat->Extent(0));
                const int ncol = static_cast<int>(mat->Extent(1));

                SEXP r_mat = guard.Protect(Rf_allocMatrix(REALSXP, nrow, ncol));
                for (int r = 0; r < nrow; ++r) {
                    for (int c = 0; c < ncol; ++c) {
                        REAL(r_mat)
                        [r + nrow * c] = (*mat)(static_cast<std::uint64_t>(r),
                            static_cast<std::uint64_t>(c));
                    }
                }
                // Add dimnames
                SEXP dimnames = guard.Protect(Rf_allocVector(VECSXP, 2));

                SEXP r_rn = guard.Protect(Rf_allocVector(STRSXP, nrow));
                auto* current_fle = fles[static_cast<std::size_t>(k)].get();
                // Nomi delle righe
                for (int i = 0; i < nrow; ++i) {
                    std::string_view row_view = current_fle->GetRowNameAt(i);
                    SET_STRING_ELT(r_rn, i,
                                   Rf_mkCharLenCE(row_view.data(), static_cast<int>(row_view.size()), CE_UTF8));
                }
                
                // Nomi delle colonne
                SEXP r_cn = guard.Protect(Rf_allocVector(STRSXP, ncol));
                for (int i = 0; i < ncol; ++i) {
                    std::string_view col_view = current_fle->GetColNameAt(i);
                    SET_STRING_ELT(r_cn, i,
                                   Rf_mkCharLenCE(col_view.data(), static_cast<int>(col_view.size()), CE_UTF8));
                }

                SET_VECTOR_ELT(dimnames, 0, r_rn);
                SET_VECTOR_ELT(dimnames, 1, r_cn);
                Rf_setAttrib(r_mat, R_DimNamesSymbol, dimnames);

                SET_VECTOR_ELT(res, k, r_mat);
                SET_STRING_ELT(res_names, k,
                    Rf_mkChar(func_names[static_cast<std::size_t>(k)].c_str()));
            }

            SET_VECTOR_ELT(
                res, n_matrix, RCreate::FromInt(guard, static_cast<int>(le_count)));
            SET_STRING_ELT(res_names, n_matrix, Rf_mkChar("n"));
            Rf_setAttrib(res, R_NamesSymbol, res_names);

            result_r = res;
        } catch (std::exception &ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }

        return result_r;
    }

    /**
     * @brief Calcola la matrice di dominanza di Brüggemann-Lerche-Sørensen (BLS) per un POSet.
     *
     * Questa funzione fa da wrapper (R/C++) per calcolare la dominanza BLS. A seconda
     * del parametro specificato, calcola la versione "Assoluta" o "Relativa" della dominanza,
     * restituendo i risultati sotto forma di matrice quadrata in R con nomi di riga e colonna
     * corrispondenti agli elementi del POSet.
     *
     * @param poset_r SEXP. Puntatore esterno (ExternalPtr) o oggetto R che wrappa l'istanza
     * C++ del `POSet` tramite la classe `POSetWrap`.
     * @param type_r  SEXP. Valore intero che determina il tipo di calcolo:
     * - 0: Calcola la dominanza BLS Assoluta (`BLSDominanceAbsolute()`).
     * - Diverso da 0: Calcola la dominanza BLS Relativa (`BLSDominanceRelative()`).
     *
     * @return SEXP   Una matrice reale (REALSXP) quadrata in R di dimensioni NxN (dove N è il
     * numero di elementi nel POSet). La matrice include i `dimnames` (nomi di riga
     * e colonna) prelevati direttamente dagli elementi del POSet.
     *
     * @exception Lancia un'eccezione a R (tramite `forward_exception_to_r`) se l'estrazione
     * del POSet fallisce, se i tipi R non sono validi o se si verificano errori
     * interni in C++.
     */
    SEXP BruggemannLercheSorensenDominance(SEXP poset_r, SEXP type_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;
        
        try {
            // 1. Estrazione del POSet dall'oggetto R
            const POSetWrap *poset_wrap = RConvert::ToPOSetWrap(poset_r);
            POSet *poset = poset_wrap->GetPOSet();
            
            // 2. Estrazione del tipo (0 = Assoluto, Altro = Relativo)
            int type = RConvert::ToInt(type_r);
            
            // 3. Calcolo della matrice di dominanza BLS
            const auto r = (type == 0) ? poset->BLSDominanceAbsolute() : poset->BLSDominanceRelative();
            
            const int nrow = static_cast<int>(r.Extent(0));
            const int ncol = static_cast<int>(r.Extent(1));
            
            // 4. Allocazione della matrice in R
            SEXP matrix_r = guard.Protect(Rf_allocMatrix(REALSXP, nrow, ncol));
            double* matrix_ptr = REAL(matrix_r); // Ottimizzazione: estrarre il puntatore fuori dal ciclo
            
            // 5. Popolamento della matrice R in formato Column-Major
            for (int col = 0; col < ncol; ++col) {
                const int offset = nrow * col;
                for (int row = 0; row < nrow; ++row) {
                    // FIX CRITICO: Usa matrix_ptr (o REAL(matrix_r)) al posto di REAL(result_r)!
                    matrix_ptr[row + offset] = r(static_cast<std::uint64_t>(row),
                                                 static_cast<std::uint64_t>(col));
                }
            }
            
            // 6. Configurazione dei Nomi di Riga e Colonna (Dimnames)
            SEXP dimnames = guard.Protect(Rf_allocVector(VECSXP, 2));
            SEXP r_names = guard.Protect(Rf_allocVector(STRSXP, nrow));
            
            for (int i = 0; i < nrow; ++i) {
                std::string_view name_view = poset->GetElementName(static_cast<std::size_t>(i));
                
                SET_STRING_ELT(
                               r_names,
                               static_cast<R_xlen_t>(i),
                               Rf_mkCharLenCE(name_view.data(), static_cast<int>(name_view.size()), CE_UTF8)
                               );
            }
            
            // Assegniamo lo stesso vettore sia alle righe che alle colonne
            SET_VECTOR_ELT(dimnames, 0, r_names);
            SET_VECTOR_ELT(dimnames, 1, r_names);
            Rf_setAttrib(matrix_r, R_DimNamesSymbol, dimnames);
            
            // 7. Ritorno del risultato
            result_r = matrix_r;
            
        } catch (std::exception &ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
}
