/**
 * @file rwrapper_separation.cpp
 * @brief R/C++ interface — separation and dimensionality reduction functions.
 *
 * Implements all R-callable (`extern "C"`) entry points for:
 * - **Separation measures**: quantify how well elements can be distinguished
 *   in a POSet (exact, fuzzy, Bubley-Dyer approximations).
 * - **Lexicographic separation**: special case for lexicographic product POSets.
 * - **Fuzzy in-betweenness**: quantifies how often element r lies "between"
 *   p and q across linear extensions.
 * - **Dimensionality reduction**: finds optimal 2D representations of profiles
 *   that preserve separation structure.
 *
 * These functions support advanced analysis of POSet structures including
 * visualization, profile optimization, and fuzzy-logic based reasoning.
 *
 * @author Alessandro Avellone
 * @version 2.1
 * @date 2025
 */

#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>

#include "rwrapper.h"
#include "my_exception.h"
#include "random.h"
#include "display_message.h"
#include "linear_generator_wrapper.h"
#include "linear_extension_generator.h"
#include "function_linear_extension.h"
#include "r_function_linear_extension.h"
#include "generic_functions.h"
#include "separation.h"
#include "dimensionality_reduction.h"
#include "r_separation.h"
#include "r_display.h"
#include "rwrapper_conversion.h"
#include "poset_wrapper.h"
#include "poset.h"


// Forward declaration of BuildExactEvaluation from rwrapper_poset_evaluation.cpp

void BuildExactEvaluation(POSet *poset,
                          const std::vector<std::string> &internal_functions,
                          const std::vector<SEXP> &external_functions,
                          std::optional<std::uint64_t> output_interval,
                          std::vector<std::unique_ptr<Tensor<double, 2>>> &eval_results,
                          std::vector<std::unique_ptr<FunctionLinearExtension>> &fles,
                          std::vector<std::string> &functions_name,
                          std::uint64_t &le_count);



// =========================================================================
// INTERNAL DISPATCH LOGIC (Hidden from other Translation Units)
// =========================================================================

namespace {
    
    /**
     * @brief Defines the natively supported t-norms and t-conorms.
     * @details Used internally to dispatch the correct C++20 Concept template.
     */
    enum class NormConormType {
        kMinimum,          ///< Gödel t-norm and t-conorm (minimum / maximum).
        kProduct,          ///< Product t-norm and Probabilistic sum t-conorm.
        kCustomRFunction   ///< User-defined custom R closures.
    };
    
    /**
     * @brief A lightweight structure to transport the user's norm/conorm selection.
     */
    struct NormConormSelection {
        NormConormType type;            ///< The mapped standard type, or kCustomRFunction.
        SEXP tnorm_raw = R_NilValue;    ///< The raw R object for the t-norm (if custom).
        SEXP tconorm_raw = R_NilValue;  ///< The raw R object for the t-conorm (if custom).
    };
    
    /**
     * @brief Parses R inputs to determine the appropriate t-norm/t-conorm configuration.
     *
     * @param tnorm_r An R object representing the t-norm (string identifier or closure).
     * @param tconorm_r An R object representing the t-conorm (used for closures).
     * @return NormConormSelection A lightweight struct containing the decision.
     * @throw MyException If the string is NA_character_ or an unknown norm name.
     */
    NormConormSelection ParseNormConorm(SEXP tnorm_r, SEXP tconorm_r) {
        NormConormSelection selection;
        
        // 1. Safely check if the input is a valid, non-empty R string
        if (Rf_isString(tnorm_r) && Rf_length(tnorm_r) > 0) {
            SEXP string_sexp = STRING_ELT(tnorm_r, 0);
            
            // 2. Protect against NA_character_
            if (string_sexp == NA_STRING) {
                throw MyException("ParseNormConorm: t-norm cannot be NA_character_");
            }
            
            const std::string name = CHAR(string_sexp);
            
            // 3. Map to native C++ implementations
            if (name == "minimum") {
                selection.type = NormConormType::kMinimum;
                return selection;
            }
            if (name == "product") {
                selection.type = NormConormType::kProduct;
                return selection;
            }
            
            throw MyException("ParseNormConorm: unknown t-norm/t-conorm '" + name + "'");
        }
        
        // 4. Fallback: treat inputs as custom R functions
        selection.type = NormConormType::kCustomRFunction;
        selection.tnorm_raw = tnorm_r;
        selection.tconorm_raw = tconorm_r;
        
        return selection;
    }
    
    /**
     * @brief Parses and validates a character vector of separation functions from R.
     *
     * @details Extracts string elements from an R character vector (STRSXP), ensuring
     * they match the supported internal separation metrics. Safely ignores NA_character_
     * and handles NULL inputs.
     *
     * @param quali_r An R character vector (e.g., c("symmetric", "asymmetricLower")).
     * @return std::vector<std::string> The validated list of function names.
     */
    inline std::vector<std::string> ParseSeparationFunctions(SEXP quali_r) {
        if (Rf_isNull(quali_r)) {
            return {};
        }
        
        const int num_elements = Rf_length(quali_r);
        
        std::vector<std::string> internal_funcs;
        internal_funcs.reserve(static_cast<std::size_t>(num_elements));
        
        for (int i = 0; i < num_elements; ++i) {
            // 3. Estrazione diretta dal vettore di caratteri (STRSXP)
            SEXP string_sexp = STRING_ELT(quali_r, i);
            
            if (string_sexp == NA_STRING) {
                continue;
            }
            
            const std::string name = CHAR(string_sexp);
            
            // 4. Validazione e inserimento
            if (name == "asymmetricLower" || name == "asymmetricUpper" || name == "symmetric") {
                internal_funcs.push_back(std::move(name));
            } else {
                throw MyException("ParseSeparationFunctions: unknown function '" + name + "'");
            }
        }
        
        return internal_funcs;
    }
    
    
} // namespace

// ***********************************************
// Separation Functions
// ***********************************************

extern "C" {
    
    // ---------------------------------------------------------------------------
    
    /**
     * @brief Computes exact separation measures via complete linear-extension enumeration.
     *
     * Delegates to `BuildExactEvaluation` with the requested separation functions.
     *
     * @param poset_r              ExternalPtrAddr of POSetWrap.
     * @param output_ogni_in_sec_r Integer scalar (progress interval in seconds), or empty.
     * @param quali_r              R list of separation function names
     *                             ("asymmetricLower", "asymmetricUpper", "symmetric").
     * @return Named R list with:
     *   - One matrix per requested function (keys = function names).
     *   - `$n` — total number of linear extensions enumerated.
     */
    SEXP ExactSeparation(SEXP poset_r, SEXP output_ogni_in_sec_r, SEXP quali_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;
        
        try {
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);
            
            const std::optional<std::uint64_t> output_ogni = RConvert::ToOptionalUInt(output_ogni_in_sec_r);
            
            auto internal_funcs = ParseSeparationFunctions(quali_r);
            const std::size_t num_funcs = internal_funcs.size();
            
            std::vector<SEXP> external_funcs(num_funcs, R_NilValue);
            std::vector<std::unique_ptr<Tensor<double, 2>>> eval_results(num_funcs);
            std::vector<std::unique_ptr<FunctionLinearExtension>> fles(num_funcs);
            std::vector<std::string> func_names(num_funcs);
            
            std::uint64_t le_count = 0;
            
            BuildExactEvaluation(poset_wrap->GetPOSet(),
                                 internal_funcs,
                                 external_funcs,
                                 output_ogni, eval_results, fles, func_names,
                                 le_count);
            
            const int num_matrices = static_cast<int>(eval_results.size());
            
            
            auto [res, res_names] = RCreate::NamedList(guard, num_matrices + 1);
            
            for (std::size_t func_idx = 0; func_idx < num_funcs; ++func_idx) {
                const auto& mat = eval_results[func_idx];
                const int nrow = static_cast<int>(mat->Extent(0));
                const int ncol = static_cast<int>(mat->Extent(1));
                
                SEXP r_mat = guard.Protect(Rf_allocMatrix(REALSXP, nrow, ncol));
                double* r_mat_ptr = REAL(r_mat);
                
                // Column-major caching
                for (int col = 0; col < ncol; ++col) {
                    for (int row = 0; row < nrow; ++row) {
                        r_mat_ptr[row + nrow * col] =
                        (*mat)(static_cast<std::uint64_t>(row), static_cast<std::uint64_t>(col));
                    }
                }
                
                
                
                SEXP row_names = guard.Protect(Rf_allocVector(STRSXP, nrow));
                for (int row = 0; row < nrow; ++row) {
                    std::string_view name = fles[func_idx]->GetRowNameAt(static_cast<std::size_t>(row));
                    SET_STRING_ELT(row_names, row,
                                   Rf_mkCharLenCE(name.data(), static_cast<int>(name.size()), CE_UTF8));
                }
                
                SEXP col_names = guard.Protect(Rf_allocVector(STRSXP, ncol));
                for (int col = 0; col < ncol; ++col) {
                    std::string_view name = fles[func_idx]->GetColNameAt(static_cast<std::size_t>(col));
                    SET_STRING_ELT(col_names, col,
                                   Rf_mkCharLenCE(name.data(), static_cast<int>(name.size()), CE_UTF8));
                }
                
                SEXP dim_names = guard.Protect(Rf_allocVector(VECSXP, 2));
                SET_VECTOR_ELT(dim_names, 0, row_names);
                SET_VECTOR_ELT(dim_names, 1, col_names);
                Rf_setAttrib(r_mat, R_DimNamesSymbol, dim_names);
                
                
                RCreate::SetListElement(res, res_names, static_cast<int>(func_idx), r_mat, func_names[func_idx].c_str());
            }
            RCreate::SetListElement(res, res_names, num_matrices, RCreate::FromInt(guard, static_cast<int>(le_count)), "n");
            Rf_setAttrib(res, R_NamesSymbol, res_names);
            
            result_r = res;
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("C++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    /**
     * @brief Computes lexicographic separation for a product of discrete variables.
     *
     * @param modalita_r  R list of character vectors; each entry is the set of
     *                    modalities (levels) for one variable.
     * @return Named R list with five matrices:
     *   - `$symmetric`, `$asymmetricLower`, `$asymmetricUpper`,
     *     `$vertical`, `$horizontal`.
     *   Each matrix has row/col names = profile labels (e.g. "level1_level2_...").
     */
    SEXP RunLexSeparation(SEXP modalita_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;
        
        try {
            const int num_vars = Rf_length(modalita_r);
            std::vector<std::vector<std::string>> modalita_string(static_cast<std::size_t>(num_vars));
            std::vector<std::uint64_t>            modalita_counts(static_cast<std::size_t>(num_vars));
            std::unordered_set<std::uint64_t>     unique_sizes;
            
            
            for (int i = 0; i < num_vars; ++i) {
                SEXP mod_r = VECTOR_ELT(modalita_r, i);
                const int mod_size = Rf_length(mod_r);
                
                unique_sizes.insert(static_cast<std::uint64_t>(mod_size));
                modalita_string[i].resize(static_cast<std::size_t>(mod_size));
                modalita_counts[i] = static_cast<std::uint64_t>(mod_size);
                
                for (int j = 0; j < mod_size; ++j) {
                    modalita_string[i][j] = CHAR(STRING_ELT(mod_r, j));
                }
            }
            
            LexSeparationResult result = (unique_sizes.size() > 1)
                        ? LexSeparationDeg(modalita_counts)
                        : LexSeparationEqDeg(static_cast<std::uint64_t>(num_vars), *unique_sizes.begin());
            const int nrow = static_cast<int>(result.sep_all.Extent(0));

            SEXP all_r   = guard.Protect(Rf_allocMatrix(REALSXP, nrow, nrow));
            SEXP lower_r = guard.Protect(Rf_allocMatrix(REALSXP, nrow, nrow));
            SEXP upper_r = guard.Protect(Rf_allocMatrix(REALSXP, nrow, nrow));
            SEXP vert_r  = guard.Protect(Rf_allocMatrix(REALSXP, nrow, nrow));
            SEXP horiz_r = guard.Protect(Rf_allocMatrix(REALSXP, nrow, nrow));
            
            // Estrazione raw pointers
            double* all_ptr   = REAL(all_r);
            double* lower_ptr = REAL(lower_r);
            double* upper_ptr = REAL(upper_r);
            double* vert_ptr  = REAL(vert_r);
            double* horiz_ptr = REAL(horiz_r);
            
            // Loop Column-Major per performance massime
            for (int col = 0; col < nrow; ++col) {
                for (int row = 0; row < nrow; ++row) {
                    const std::uint64_t ru = static_cast<std::uint64_t>(row);
                    const std::uint64_t cu = static_cast<std::uint64_t>(col);
                    
                    const int idx = row + nrow * col;
                    
                    // Accesso diretto ai membri della struct (niente std::get)
                    all_ptr[idx]   = result.sep_all(ru, cu);
                    lower_ptr[idx] = result.sep_lower(ru, cu);
                    upper_ptr[idx] = result.sep_upper(ru, cu);
                    vert_ptr[idx]  = result.sep_vertical(ru, cu);
                    horiz_ptr[idx] = result.sep_horizontal(ru, cu);
                }
            }
            
            SEXP row_names = guard.Protect(Rf_allocVector(STRSXP, nrow));
            SEXP col_names = guard.Protect(Rf_allocVector(STRSXP, nrow));
            
            for (int pid = 0; pid < nrow; ++pid) {
                std::string label;
                const auto& profile = result.profili[static_cast<std::size_t>(pid)];
                
                for (std::size_t vid = 0; vid < modalita_string.size(); ++vid) {
                    const auto mod_id = profile[vid];
                    const auto& mod   = modalita_string[vid][static_cast<std::size_t>(mod_id)];
                    label += (vid < modalita_string.size() - 1) ? mod + "_" : mod;
                }
                
                SEXP r_label = Rf_mkChar(label.c_str());
                SET_STRING_ELT(row_names, pid, r_label);
                SET_STRING_ELT(col_names, pid, r_label);
            }
            
            SEXP names_r = guard.Protect(Rf_allocVector(VECSXP, 2));
            SET_VECTOR_ELT(names_r, 0, row_names);
            SET_VECTOR_ELT(names_r, 1, col_names);
            
            Rf_setAttrib(all_r,   R_DimNamesSymbol, names_r);
            Rf_setAttrib(lower_r, R_DimNamesSymbol, names_r);
            Rf_setAttrib(upper_r, R_DimNamesSymbol, names_r);
            Rf_setAttrib(vert_r,  R_DimNamesSymbol, names_r);
            Rf_setAttrib(horiz_r, R_DimNamesSymbol, names_r);
            
            auto [res, res_names] = RCreate::NamedList(guard, 5);
            RCreate::SetListElement(res, res_names, 0, all_r,   "symmetric");
            RCreate::SetListElement(res, res_names, 1, lower_r, "asymmetricLower");
            RCreate::SetListElement(res, res_names, 2, upper_r, "asymmetricUpper");
            RCreate::SetListElement(res, res_names, 3, vert_r,  "vertical");
            RCreate::SetListElement(res, res_names, 4, horiz_r, "horizontal");
            
            Rf_setAttrib(res, R_NamesSymbol, res_names);
            
            result_r = res;
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    /**
     * @brief Computes the Mutual Ranking Probability (MRP) matrix for a lexicographic product.
     *
     * @param modalita_r  R list of character vectors (variable modalities).
     * @return R real matrix with row/col names = profile labels.
     */
    SEXP RunLexMRP(SEXP modalita_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;
        
        try {
            const int num_vars = Rf_length(modalita_r);
            std::vector<std::vector<std::string>> modalita_string(static_cast<std::size_t>(num_vars));
            std::vector<std::uint64_t>            modalita_counts(static_cast<std::size_t>(num_vars));
            
            for (int i = 0; i < num_vars; ++i) {
                SEXP mod_r = VECTOR_ELT(modalita_r, i);
                const int mod_size = Rf_length(mod_r);
                
                modalita_string[i].resize(static_cast<std::size_t>(mod_size));
                modalita_counts[i] = static_cast<std::uint64_t>(mod_size);
                
                for (int j = 0; j < mod_size; ++j) {
                    // C/R API usa CHAR(), non R_CHAR()
                    modalita_string[i][j] = CHAR(STRING_ELT(mod_r, j));
                }
            }
            
            auto result = LexMrp(modalita_counts);
            const int nrow = static_cast<int>(result.mrp.Extent(0));
            
            result_r = RCreate::FromDoubleMatrix(guard, result.mrp);
            
            SEXP row_names = guard.Protect(Rf_allocVector(STRSXP, nrow));
            SEXP col_names = guard.Protect(Rf_allocVector(STRSXP, nrow));
            
            for (int pid = 0; pid < nrow; ++pid) {
                std::string label;
                const auto& profile = result.profili[static_cast<std::size_t>(pid)];
                
                for (std::size_t vid = 0; vid < modalita_string.size(); ++vid) {
                    const auto mod_id = profile[vid];
                    const auto& mod   = modalita_string[vid][static_cast<std::size_t>(mod_id)];
                    label += (vid < modalita_string.size() - 1) ? mod + "_" : mod;
                }
                
                SEXP r_label = Rf_mkChar(label.c_str());
                SET_STRING_ELT(row_names, pid, r_label);
                SET_STRING_ELT(col_names, pid, r_label);
            }

            RCreate::AttachDimNames(guard, result_r, row_names, col_names);
            
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    /**
     * @brief Builds an R External Pointer wrapping a BubleyDyerSeparationGenerator.
     *
     * @details Extracts a POSet object from an R external pointer, initializes a
     * random number generator based on an optional seed, and creates a
     * BubleyDyerEvaluationGenerator via its factory method. The resulting C++
     * object is wrapped in a new R external pointer with a registered C finalizer
     * to ensure safe memory cleanup when the object goes out of scope in R.
     *
     * @param poset_r R external pointer (EXTPTRSXP) to a valid POSetWrap instance.
     * @param seed_r R integer vector (INTSXP) for the RNG seed. If empty, the global RNG is used.
     * @param quali_r R object defining the separation functions to be parsed.
     * @return SEXP An R external pointer wrapping the newly created generator.
     */
    SEXP BuildBubleyDyerSeparationGenerator(SEXP poset_r, SEXP seed_r, SEXP quali_r) {
        SEXP result_r = R_NilValue;
        RProtectGuard guard;
        
        try {
            // 1. Estrazione del raw pointer
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);

            // 2. Estrazione del seed
            std::optional<std::uint64_t> seed = RConvert::ToOptionalSeed(seed_r);
            if (!seed.has_value()) {
                seed = Random::GLOBAL.RndNextInt(0, std::numeric_limits<std::uint64_t>::max());
            }
            
            // 3. Estrazione delle funzioni
            auto internal_funcs = ParseSeparationFunctions(quali_r);
            std::vector<SEXP> external_funcs(internal_funcs.size(), R_NilValue);
            
            // 4. Costruzione dell'oggetto C++
            auto generator = BubleyDyerEvaluationGenerator::BuildBubleyDyerEvaluationGenerator(poset_wrap, seed.value(), internal_funcs, external_funcs);

            // 5. Creazione del puntatore esterno in R
            SEXP ptr_r = RCreate::WrapExternalPtr(guard, std::move(generator));

            // 6. Creazione dello "scudo" per il Garbage Collector
            R_SetExternalPtrProtected(ptr_r, poset_r);

            // 7. Il seme usato viene restituito a R come stringa: reimmetterlo
            // riproduce esattamente la stessa sessione.
            result_r = RCreate::GeneratorWithSeed(guard, ptr_r, seed.value());
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    // ---------------------------------------------------------------------------
    
    /**
     * @brief Runs Bubley-Dyer separation (thin alias to BubleyDyerEvaluation).
     *
     * @param generator_r          ExternalPtrAddr of BubleyDyerEvaluationGenerator.
     * @param quante_r             Integer scalar (sample size), or empty.
     * @param errore_r             Real scalar (target error), or empty.
     * @param output_ogni_in_sec_r Integer scalar (progress interval), or empty.
     * @return Same return value as `BubleyDyerEvaluation`.
     */
    SEXP BubleyDyerSeparation(SEXP generator_r, SEXP quante_r, SEXP errore_r, SEXP output_ogni_in_sec_r) {
        return BubleyDyerEvaluation(generator_r, quante_r, errore_r, output_ogni_in_sec_r);
    }
    
    /**
     * @brief Computes general fuzzy separation matrices.
     *
     * @details Extracts the dominance matrix and requested separation types.
     * It uses a compile-time template dispatch (via switch statement) to strictly
     * enforce the C++20 NormConormFunc concept, ensuring maximum inlining and
     * zero virtual overhead during HPC loops.
     *
     * @param dominance_matrix_r R numeric matrix (REALSXP) representing dominance.
     * @param tnorm_r R identifier for the t-norm.
     * @param tconorm_r R identifier for the t-conorm.
     * @param quali_r R character vector (STRSXP) of requested outputs (e.g., "symmetric").
     * @return SEXP A named R list containing the requested matrices.
     */
    SEXP FuzzySeparation(SEXP dominance_matrix_r, SEXP tnorm_r, SEXP tconorm_r, SEXP quali_r) {
        SEXP result_r = R_NilValue;
        RProtectGuard guard;
        
        try {
            
            SEXP dimnames_r = Rf_getAttrib(dominance_matrix_r, R_DimNamesSymbol);
            if (Rf_isNull(dimnames_r)) {
                throw MyException("The dominance matrix must have dimnames (row and column names assigned).");
            }
            
            SEXP row_names = R_NilValue;
            SEXP col_names = R_NilValue;
            
            if (dimnames_r != R_NilValue && Rf_length(dimnames_r) >= 2) {
                row_names = VECTOR_ELT(dimnames_r, 0);
                col_names = VECTOR_ELT(dimnames_r, 1);
            }
            
            auto dominance = RConvert::ToDoubleMatrix(dominance_matrix_r);

            bool do_all = false, do_lower = false, do_upper = false;
            const int n_quali = Rf_length(quali_r);
            
            for (int k = 0; k < n_quali; ++k) {
                const std::string name = CHAR(STRING_ELT(quali_r, k));
                
                if (name == "asymmetricLower") {
                    do_lower = true;
                } else if (name == "asymmetricUpper") {
                    do_upper = true;
                } else if (name == "symmetric") {
                    do_all = true;
                } else {
                    throw std::invalid_argument("FuzzySeparation: unknown request '" + name + "'");
                }
            }
            
            NormConormSelection selection = ParseNormConorm(tnorm_r, tconorm_r);
            
            auto run_separation = [&]() -> GeneralSeparationResult {
                switch (selection.type) {
                    case NormConormType::kMinimum:
                        return GeneralSeparation(dominance, MinNormConorm{}, MaxNormConorm{},
                                                 do_all, do_lower, do_upper);
                    case NormConormType::kProduct:
                        return GeneralSeparation(dominance, ProdNormConorm{}, ProbNormConorm{},
                                                 do_all, do_lower, do_upper);
                    case NormConormType::kCustomRFunction: {
                        // PRE-ALLOCAZIONE per la Norma
                        SEXP norm_arg1 = guard.Protect(Rf_allocVector(REALSXP, 1));
                        SEXP norm_arg2 = guard.Protect(Rf_allocVector(REALSXP, 1));
                        SEXP norm_call = guard.Protect(Rf_lang3(tnorm_r, norm_arg1, norm_arg2));
                        RCustomNormConorm customNorm(norm_arg1, norm_arg2, norm_call, R_GlobalEnv);
                        
                        // PRE-ALLOCAZIONE per la Conorma
                        SEXP conorm_arg1 = guard.Protect(Rf_allocVector(REALSXP, 1));
                        SEXP conorm_arg2 = guard.Protect(Rf_allocVector(REALSXP, 1));
                        SEXP conorm_call = guard.Protect(Rf_lang3(tconorm_r, conorm_arg1, conorm_arg2));
                        RCustomNormConorm customConorm(conorm_arg1, conorm_arg2, conorm_call, R_GlobalEnv);
                        
                        // Chiamata zero-overhead col nostro functor template
                        return GeneralSeparation(dominance, customNorm, customConorm,
                                                 do_all, do_lower, do_upper);
                    }
                    default:
                        throw std::invalid_argument("Unsupported Norm/Conorm type.");
                }
            };
            
            GeneralSeparationResult risultato = run_separation();
            
            SEXP all_r = R_NilValue, lower_r = R_NilValue, upper_r = R_NilValue;
            
            if (do_all) {
                all_r = RCreate::FromDoubleMatrix(guard, risultato.sep_all);
                if (row_names != R_NilValue) {
                    RCreate::AttachDimNames(guard, all_r, row_names, col_names);
                }
            }
            if (do_lower) {
                lower_r = RCreate::FromDoubleMatrix(guard, risultato.sep_lower);
                if (row_names != R_NilValue) {
                    RCreate::AttachDimNames(guard, lower_r, row_names, col_names);
                }
            }
            if (do_upper) {
                upper_r = RCreate::FromDoubleMatrix(guard, risultato.sep_upper);
                if (row_names != R_NilValue) {
                    RCreate::AttachDimNames(guard, upper_r, row_names, col_names);
                }
            }
            
            auto [res, res_names] = RCreate::NamedList(guard, 3);
            RCreate::SetListElement(res, res_names, 0, all_r,   "symmetric");
            RCreate::SetListElement(res, res_names, 1, lower_r, "asymmetricLower");
            RCreate::SetListElement(res, res_names, 2, upper_r, "asymmetricUpper");
            Rf_setAttrib(res, R_NamesSymbol, res_names);
            
            result_r = res;
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    /**
     * @brief Computes the Fuzzy In-Betweenness for a given dominance matrix.
     * * @details This function calculates various types of fuzzy in-betweenness
     * (symmetric, asymmetric lower, asymmetric upper) based on user requests.
     * It uses a highly optimized $O(N^3)$ nested loop with mirrored assignments
     * to halve the iteration space. Furthermore, it employs Template Dispatching
     * via generic lambdas to resolve t-norm and t-conorm operations at compile time,
     * guaranteeing zero-overhead inside the deepest loops.
     * * @param dominance_matrix_r The dominance matrix from R.
     * @param tnorm_r The chosen t-norm (e.g., Minimum, Product).
     * @param tconorm_r The chosen t-conorm (e.g., Maximum, Probabilistic Sum).
     * @param quali_r An R character vector specifying which matrices to compute
     * ("symmetric", "asymmetricLower", "asymmetricUpper").
     * @return SEXP A named R list containing 3D arrays (tensors) of the requested
     * in-betweenness calculations.
     */
    SEXP FuzzyInBetweenness(SEXP dominance_matrix_r, SEXP tnorm_r, SEXP tconorm_r, SEXP quali_r) {
        SEXP result_r = R_NilValue;
        RProtectGuard guard;
        
        try {
            const int nrow = INTEGER(Rf_getAttrib(dominance_matrix_r, R_DimSymbol))[0];
            auto dimnames_r = Rf_getAttrib(dominance_matrix_r, R_DimNamesSymbol);
            
            if (Rf_isNull(dimnames_r)) {
                throw MyException("The dominance matrix must have dimnames (row and column names assigned).");
            }
            
            SEXP row_names  = VECTOR_ELT(dimnames_r, 0);
            SEXP col_names  = VECTOR_ELT(dimnames_r, 1);
            
            auto dominance   = RConvert::ToDoubleMatrix(dominance_matrix_r);
            
            bool do_all = false, do_lower = false, do_upper = false;
            
            using Tensor3D = std::vector<std::vector<std::vector<double>>>;
            
            std::optional<Tensor3D> inbet_lower;
            std::optional<Tensor3D> inbet_upper;
            std::optional<Tensor3D> inbet_all;
            
            const int n_quali = Rf_length(quali_r);
            const auto sz_nrow = static_cast<std::size_t>(nrow);
            
            for (int k = 0; k < n_quali; ++k) {
                const std::string name = R_CHAR(STRING_ELT(quali_r, k));
                
                // Inizializzazione pulita tramite .emplace() con 0.0 di default
                if (name == "asymmetricLower") {
                    do_lower = true;
                    inbet_lower.emplace(sz_nrow, std::vector<std::vector<double>>(sz_nrow, std::vector<double>(sz_nrow, 0.0)));
                } else if (name == "asymmetricUpper") {
                    do_upper = true;
                    inbet_upper.emplace(sz_nrow, std::vector<std::vector<double>>(sz_nrow, std::vector<double>(sz_nrow, 0.0)));
                } else if (name == "symmetric") {
                    do_all = true;
                    inbet_all.emplace(sz_nrow, std::vector<std::vector<double>>(sz_nrow, std::vector<double>(sz_nrow, 0.0)));
                } else {
                    throw MyException("FuzzyInBetweenness: unknown '" + name + "'");
                }
            }

            // --- TEMPLATE DISPATCHING ---
            // Definiamo una Lambda Generica.
            // In questo modo GeneralFuzzyInBetweenness capirà a compile-time il tipo esatto di norma!
            auto compute_inbetweenness = [&](auto norm, auto conorm) {
                for (std::size_t pi = 0; pi < sz_nrow; ++pi) {
                    for (std::size_t qi = pi + 1; qi < sz_nrow; ++qi) {
                        for (std::size_t ri = 0; ri < sz_nrow; ++ri) {
                            
                            double finb_prq = 0.0, finb_qrp = 0.0, finbqrp = 0.0;
                            // Chiamata zero-overhead, norm e conorm sono risolti a compile-time!
                            GeneralFuzzyInBetweenness(pi, qi, ri, dominance, norm, conorm, finb_prq, finb_qrp, finbqrp);
                            
                            if (do_lower) {
                                (*inbet_lower)[pi][qi][ri] = finb_prq;
                                (*inbet_lower)[qi][pi][ri] = finb_qrp;
                            }
                            if (do_upper) {
                                (*inbet_upper)[pi][qi][ri] = finb_qrp;
                                (*inbet_upper)[qi][pi][ri] = finb_prq;
                            }
                            if (do_all) {
                                (*inbet_all)[pi][qi][ri] = finbqrp;
                                (*inbet_all)[qi][pi][ri] = finbqrp;
                            }
                        }
                    }
                }
            };
            
            // Switch per instanziare la lambda col tipo corretto
            NormConormSelection selection = ParseNormConorm(tnorm_r, tconorm_r);
            switch (selection.type) {
                case NormConormType::kMinimum:
                    compute_inbetweenness(MinNormConorm{}, MaxNormConorm{});
                    break;
                case NormConormType::kProduct:
                    compute_inbetweenness(ProdNormConorm{}, ProbNormConorm{});
                    break;
                case NormConormType::kCustomRFunction: {
                    // PRE-ALLOCAZIONE: Costruiamo lo scheletro della chiamata per la Norma
                    SEXP norm_arg1 = guard.Protect(Rf_allocVector(REALSXP, 1));
                    SEXP norm_arg2 = guard.Protect(Rf_allocVector(REALSXP, 1));
                    // Crea l'espressione: tnorm_r(norm_arg1, norm_arg2)
                    SEXP norm_call = guard.Protect(Rf_lang3(tnorm_r, norm_arg1, norm_arg2));
                    RCustomNormConorm customNorm(norm_arg1, norm_arg2, norm_call, R_GlobalEnv);
                    
                    // PRE-ALLOCAZIONE: Costruiamo lo scheletro della chiamata per la Conorma
                    SEXP conorm_arg1 = guard.Protect(Rf_allocVector(REALSXP, 1));
                    SEXP conorm_arg2 = guard.Protect(Rf_allocVector(REALSXP, 1));
                    // Crea l'espressione: tconorm_r(conorm_arg1, conorm_arg2)
                    SEXP conorm_call = guard.Protect(Rf_lang3(tconorm_r, conorm_arg1, conorm_arg2));
                    RCustomNormConorm customConorm(conorm_arg1, conorm_arg2, conorm_call, R_GlobalEnv);
                    
                    // Esegue la tua lambda (il triplo loop) usando le funzioni R
                    compute_inbetweenness(customNorm, customConorm);
                    
                    // Nessun UNPROTECT manuale!
                    // Quando il blocco (case) termina, 'guard' esce di scope e il suo distruttore
                    // chiamerà automaticamente UNPROTECT(6).
                    break;
                }
                default:
                    throw std::invalid_argument("Unsupported Norm/Conorm type.");
            }
            
            SEXP all_r = R_NilValue, lower_r = R_NilValue, upper_r = R_NilValue;
            
            if (inbet_all) {
                all_r = RCreate::From3DDoubleVector(guard, *inbet_all);
                RCreate::AttachDimNames(guard, all_r, row_names, col_names);
            }
            if (inbet_lower) {
                lower_r = RCreate::From3DDoubleVector(guard, *inbet_lower);
                RCreate::AttachDimNames(guard, lower_r, row_names, col_names);
            }
            if (inbet_upper) {
                upper_r = RCreate::From3DDoubleVector(guard, *inbet_upper);
                RCreate::AttachDimNames(guard, upper_r, row_names, col_names);
            }
            
            auto [res, res_names] = RCreate::NamedList(guard, 3);
            RCreate::SetListElement(res, res_names, 0, all_r,   "symmetric");
            RCreate::SetListElement(res, res_names, 1, lower_r, "asymmetricLower");
            RCreate::SetListElement(res, res_names, 2, upper_r, "asymmetricUpper");
            Rf_setAttrib(res, R_NamesSymbol, res_names);
            
            result_r = res;
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    /**
     * @brief Wrapper ad alte prestazioni per l'algoritmo di Riduzione della Dimensionalità.
     *
     * @details Questa funzione funge da ponte (bridge) tra l'ambiente R e l'engine C++.
     * Applica diverse ottimizzazioni HPC:
     * 1. **Cache Locality**: La matrice dei profili di R viene letta per colonne (layout nativo di R)
     * per minimizzare i cache-miss durante la generazione delle chiavi bitwise.
     * 2. **Zero-Allocation**: Utilizza `std::string_view` per mappare direttamente le stringhe
     * di R senza duplicare memoria sull'heap.
     * 3. **Efficient Writing**: La matrice dei risultati viene popolata seguendo il layout
     * memory-contiguous di R.
     * 4. **Exception Safety**: Traduce le eccezioni C++ in messaggi di errore leggibili in R
     * tramite `forward_exception_to_r`.
     *
     * @param profile_r Matrice R (INTSXP) dei profili (0/1).
     * @param weights_r Vettore R (REALSXP) dei pesi associati ai profili.
     * @param loss_r Scalare R (STRSXP) indicante la funzione di loss da utilizzare.
     * @param lpom_strategy_r Scalare R (INTSXP) per la strategia LPOM (0 = absolute, 1 = relative).
     * @param output_ogni_in_sec_r Frequenza di aggiornamento del display (optional).
     * @param thread_percentage_r Percentuale di core CPU da utilizzare (0.0 - 1.0).
     *
     * @return SEXP Lista nominata contenente: allLoss, variablesPriority, bestLossValue,
     * bestVariablesPriority, bestRepresentation.
     */
    SEXP RunDimensionalityReduction(SEXP profile_r,
                                    SEXP weights_r,
                                    SEXP loss_r,
                                    SEXP lpom_strategy_r,
                                    SEXP output_ogni_in_sec_r,
                                    SEXP thread_percentage_r) {
        SEXP result_r = R_NilValue;
        RProtectGuard guard;
        
        try {
            // Estrazione dimensioni matrice profili
            SEXP profile_r_dim = Rf_getAttrib(profile_r, R_DimSymbol);
            if (Rf_isNull(profile_r_dim)) {
                throw MyException("profile must be a matrix with dimensions.");
            }
            
            const R_len_t nrow = INTEGER(profile_r_dim)[0];
            const R_len_t ncol = INTEGER(profile_r_dim)[1];
            
            // ---------------------------------------------------------------------
            // Generazione chiavi ottimizzata per la Cache
            // ---------------------------------------------------------------------
            std::vector<std::uint64_t> keys(static_cast<std::size_t>(nrow), 0);
            const int* prof_ptr = INTEGER(profile_r);
            
            for (R_len_t col = 0; col < ncol; ++col) {
                const int offset = nrow * col;
                for (R_len_t row = 0; row < nrow; ++row) {
                    const std::uint64_t v = static_cast<std::uint64_t>(prof_ptr[row + offset]);
                    keys[static_cast<std::size_t>(row)] = (keys[static_cast<std::size_t>(row)] << 1) | v;
                }
            }
            
            std::unordered_map<std::uint64_t, double> weights;
            weights.reserve(static_cast<std::size_t>(nrow));
            
            const double* w_ptr = REAL(weights_r);
            for (R_len_t row = 0; row < nrow; ++row) {
                weights[keys[static_cast<std::size_t>(row)]] = w_ptr[row];
            }
            
            const std::string_view loss_str      = CHAR(STRING_ELT(loss_r, 0));
            const int lpom_strategy = Rf_asInteger(lpom_strategy_r);
            
            const std::uint64_t total_pairs = (ncol * (ncol - 1)) / 2;
            const std::uint64_t permutations_x_pair = static_cast<std::uint64_t>(std::tgamma(ncol - 1));
            const std::uint64_t total_le = permutations_x_pair * total_pairs;
            
            
            std::atomic<std::uint64_t> all_le_count{0};
            
            // Gestione messaggi di progresso
            const std::optional<std::uint64_t> output_ogni = RConvert::ToOptionalUInt(output_ogni_in_sec_r);
            std::unique_ptr<DisplayMessage> display_message;
            if (!output_ogni.has_value()) {
                display_message = std::make_unique<DisplayMessageNull>();
            } else {
                display_message = std::make_unique<DisplayMessageEvaluationBDRThreads>(
                                                                                all_le_count,
                                                                                total_le,
                                                                                output_ogni.value()
                                                                                );
            }
            
            const double thread_percentage = RConvert::ToDouble(thread_percentage_r);
            
            // Esecuzione dell'algoritmo Core
            DimensionalityReductionResult result(all_le_count);
            ExactDimensionalityReduction(weights,
                                         static_cast<std::uint64_t>(ncol),
                                         loss_str,
                                         lpom_strategy,
                                         display_message.get(),
                                         thread_percentage,
                                         result);
            
            // ---------------------------------------------------------------------
            // Esportazione Risultati verso R
            // ---------------------------------------------------------------------
            auto le   = result.GetLeElaborated();
            auto loss = result.GetLossValues();
            
            const int n_perm = static_cast<int>(le.size());
            SEXP r_le = guard.Protect(Rf_allocMatrix(INTSXP, n_perm, ncol));
            int* r_le_ptr = INTEGER(r_le);
            
            for (int c = 0; c < ncol; ++c) {
                const int offset = n_perm * c;
                for (int r = 0; r < n_perm; ++r) {
                    // Convertiamo a 1-indexing per R
                    r_le_ptr[r + offset] = static_cast<int>(le[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)]) + 1;
                }
            }
            
            // Vettore delle loss
            SEXP r_loss = guard.Protect(Rf_allocVector(REALSXP, static_cast<int>(loss.size())));
            std::copy(loss.begin(), loss.end(), REAL(r_loss));
            
            // Valori migliori trovati
            SEXP r_best_loss = RCreate::FromDouble(guard, result.GetBestLoss());
            auto best_le_vec = result.GetBestLe();
            
            // CORREZIONE: Uso di INTSXP per il best_le in modo che sia coerente con r_le
            SEXP r_best_le = guard.Protect(Rf_allocVector(INTSXP, static_cast<int>(best_le_vec.size())));
            int* r_best_le_ptr = INTEGER(r_best_le);
            for (std::size_t k = 0; k < best_le_vec.size(); ++k) {
                r_best_le_ptr[k] = static_cast<int>(best_le_vec[k]) + 1; // 1-indexed
            }
            
            // Costruzione Rappresentazione Ottimale (Named List)
            auto best_prof_res = result.GetBestProfileResults();
            const int n_prof   = static_cast<int>(best_prof_res.size());
            auto [repr, repr_names] = RCreate::NamedList(guard, 5);
            
            SEXP r_profiles = guard.Protect(Rf_allocVector(INTSXP, n_prof));
            SEXP r_x_pos    = guard.Protect(Rf_allocVector(INTSXP, n_prof));
            SEXP r_y_pos    = guard.Protect(Rf_allocVector(INTSXP, n_prof));
            SEXP r_weights  = guard.Protect(Rf_allocVector(REALSXP, n_prof));
            SEXP r_error    = guard.Protect(Rf_allocVector(REALSXP, n_prof));
            
            int idx = 0;
            for (const auto& [prof_id, tup] : best_prof_res) {
                INTEGER(r_profiles)[idx] = static_cast<int>(prof_id);
                INTEGER(r_x_pos)[idx]    = static_cast<int>(std::get<0>(tup));
                INTEGER(r_y_pos)[idx]    = static_cast<int>(std::get<1>(tup));
                REAL(r_weights)[idx]     = std::get<2>(tup);
                REAL(r_error)[idx]       = std::get<3>(tup);
                ++idx;
            }
            
            RCreate::SetListElement(repr, repr_names, 0, r_profiles, "profiles");
            RCreate::SetListElement(repr, repr_names, 1, r_x_pos,    "x");
            RCreate::SetListElement(repr, repr_names, 2, r_y_pos,    "y");
            RCreate::SetListElement(repr, repr_names, 3, r_weights,  "weights");
            RCreate::SetListElement(repr, repr_names, 4, r_error,    "error");
            Rf_setAttrib(repr, R_NamesSymbol, repr_names);
            
            // Lista finale dei risultati
            auto [res, res_names] = RCreate::NamedList(guard, 5);
            RCreate::SetListElement(res, res_names, 0, r_loss,      "allLoss");
            RCreate::SetListElement(res, res_names, 1, r_le,        "variablesPriority");
            RCreate::SetListElement(res, res_names, 2, r_best_loss, "bestLossValue");
            RCreate::SetListElement(res, res_names, 3, r_best_le,   "bestVariablesPriority");
            RCreate::SetListElement(res, res_names, 4, repr,        "bestRepresentation");
            Rf_setAttrib(res, R_NamesSymbol, res_names);
            
            result_r = res;
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("C++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    
    /**
     * @brief Optimizes 2D profile representation with a fixed variable priority.
     *
     * Similar to `RunDimensionalityReduction` but tests only the single provided
     * permutation (much faster).
     *
     * @param profile_r           R integer matrix (profiles × variables).
     * @param weights_r           R real vector (profile weights).
     * @param loss_r              R string scalar (loss function name).
     * @param lpom_strategy_r     R integer scalar (0 for absolute, 1 for relative).
     * @param variable_priority_r Integer vector (1-indexed variable permutation).
     * @return Named R list with:
     * - `$lossValue` — real scalar (loss for this permutation).
     * - `$variablesPriority` — integer vector (same as input, 1-indexed).
     * - `$representation` — named list with profile 2D coordinates and errors.
     */
    SEXP RunBidimentionalPosetRepresentation(SEXP profile_r, SEXP weights_r,
                                             SEXP loss_r, SEXP lpom_strategy_r, SEXP variable_priority_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;
        
        try {
            // 1. CORREZIONE: Controllo di Sicurezza sulle Dimensioni
            SEXP profile_r_dim = Rf_getAttrib(profile_r, R_DimSymbol);
            if (Rf_isNull(profile_r_dim)) {
                throw MyException("profile_r must be a matrix with dimensions.");
            }
            
            // Uso di R_len_t (il tipo corretto per le dimensioni nell'API C di R)
            const R_len_t nrow = INTEGER(profile_r_dim)[0];
            const R_len_t ncol = INTEGER(profile_r_dim)[1];
            
            // ---------------------------------------------------------------------
            // Generazione chiavi ottimizzata per la Cache (Column-Major)
            // ---------------------------------------------------------------------
            std::vector<std::uint64_t> keys(static_cast<std::size_t>(nrow), 0);
            const int* prof_ptr = INTEGER(profile_r);
            
            for (R_len_t col = 0; col < ncol; ++col) {
                const int offset = nrow * col;
                for (R_len_t row = 0; row < nrow; ++row) {
                    // Accesso lineare sequenziale alla memoria nativa di R
                    const std::uint64_t v = static_cast<std::uint64_t>(prof_ptr[row + offset]);
                    keys[static_cast<std::size_t>(row)] = (keys[static_cast<std::size_t>(row)] << 1) | v;
                }
            }
            
            // ---------------------------------------------------------------------
            // Costruzione Mappa Pesi
            // ---------------------------------------------------------------------
            std::unordered_map<std::uint64_t, double> weights;
            weights.reserve(static_cast<std::size_t>(nrow)); // Pre-allocazione essenziale
            
            const double* w_ptr = REAL(weights_r);
            for (R_len_t row = 0; row < nrow; ++row) {
                weights[keys[static_cast<std::size_t>(row)]] = w_ptr[row];
            }
            
            // Zero-Allocation Strings
            const std::string_view loss_str      = CHAR(STRING_ELT(loss_r, 0));
            const int lpom_strategy = Rf_asInteger(lpom_strategy_r);
            
            // Estrazione e conversione (1-indexed -> 0-indexed per C++)
            const R_len_t var_prio_len = Rf_length(variable_priority_r);
            std::vector<std::uint64_t> variable_priority(static_cast<std::size_t>(var_prio_len));
            
            int* var_prio_ptr = INTEGER(variable_priority_r);
            for (R_len_t i = 0; i < var_prio_len; ++i) {
                variable_priority[static_cast<std::size_t>(i)] = static_cast<std::uint64_t>(var_prio_ptr[i] - 1);
            }
            
            // Predisposizione dell'Atomic Counter
            std::atomic<std::uint64_t> all_le_count{0};
            DimensionalityReductionResult result(all_le_count);
            
            // Esecuzione C++ Engine
            BidimentionalPosetRepresentation(weights,
                                             static_cast<std::uint64_t>(ncol),
                                             loss_str,
                                             lpom_strategy,
                                             variable_priority,
                                             result);
            
            // Estrazione dei risultati
            const double best_loss = result.GetBestLoss();
            auto best_le           = result.GetBestLe();
            auto best_prof_res     = result.GetBestProfileResults();
            
            SEXP r_best_loss = RCreate::FromDouble(guard, best_loss);
            
            // 2. CORREZIONE: Uso di INTSXP per l'esportazione degli indici in R
            SEXP r_best_le = guard.Protect(Rf_allocVector(INTSXP, static_cast<int>(best_le.size())));
            int* r_best_le_ptr = INTEGER(r_best_le);
            for (std::size_t k = 0; k < best_le.size(); ++k) {
                r_best_le_ptr[k] = static_cast<int>(best_le[k]) + 1; // Ri-conversione a 1-indexed
            }
            
            // ---------------------------------------------------------------------
            // Esportazione Best Representation in R List
            // ---------------------------------------------------------------------
            const int n_prof = static_cast<int>(best_prof_res.size());
            auto [repr, repr_names] = RCreate::NamedList(guard, 5);
            
            SEXP r_profiles = guard.Protect(Rf_allocVector(INTSXP, n_prof));
            SEXP r_x_pos    = guard.Protect(Rf_allocVector(INTSXP, n_prof));
            SEXP r_y_pos    = guard.Protect(Rf_allocVector(INTSXP, n_prof));
            SEXP r_weights  = guard.Protect(Rf_allocVector(REALSXP, n_prof));
            SEXP r_error    = guard.Protect(Rf_allocVector(REALSXP, n_prof));
            
            int idx = 0;
            for (const auto& [prof_id, tup] : best_prof_res) {
                INTEGER(r_profiles)[idx] = static_cast<int>(prof_id);
                INTEGER(r_x_pos)[idx]    = static_cast<int>(std::get<0>(tup));
                INTEGER(r_y_pos)[idx]    = static_cast<int>(std::get<1>(tup));
                REAL(r_weights)[idx]     = std::get<2>(tup);
                REAL(r_error)[idx]       = std::get<3>(tup);
                ++idx;
            }
            
            RCreate::SetListElement(repr, repr_names, 0, r_profiles, "profiles");
            RCreate::SetListElement(repr, repr_names, 1, r_x_pos,    "x");
            RCreate::SetListElement(repr, repr_names, 2, r_y_pos,    "y");
            RCreate::SetListElement(repr, repr_names, 3, r_weights,  "weights");
            RCreate::SetListElement(repr, repr_names, 4, r_error,    "error");
            Rf_setAttrib(repr, R_NamesSymbol, repr_names);
            
            // Lista di Risposta Finale
            auto [res, res_names] = RCreate::NamedList(guard, 3);
            RCreate::SetListElement(res, res_names, 0, r_best_loss, "lossValue");
            RCreate::SetListElement(res, res_names, 1, r_best_le,   "variablesPriority");
            RCreate::SetListElement(res, res_names, 2, repr,        "representation");
            Rf_setAttrib(res, R_NamesSymbol, res_names);
            
            result_r = res;
            
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
}

