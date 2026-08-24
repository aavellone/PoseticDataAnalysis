/**
 * @file linear_generator_wrapper.cpp
 * @brief Implementation of wrappers optimized for C++20 and HPC environment.
 *
 * @details Defines the creation logic (Factories) and extraction cycles.
 * The code leverages Return Value Optimization (RVO) and Move Semantics
 * to eliminate unnecessary copies of large data structures.
 */

#include "linear_generator_wrapper.h"

#include <chrono>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <optional>

#include "function_linear_extension.h"
#include "function_linear_extension_average_height.h"
#include "function_linear_extension_mutual_ranking_probability.h"
#include "function_linear_extension_separation_asymmetric_lower.h"
#include "function_linear_extension_separation_asymmetric_upper.h"
#include "function_linear_extension_separation_symmetric.h"
#include "r_function_linear_extension.h"
#include "linear_extension_generator.h"
#include "my_exception.h"
#include "poset_wrapper.h"
#include "random.h"
#include "poset.h"
#include "rwrapper_conversion.h"

/**
 * @brief Factory for Bubley-Dyer.
 * @details Initializes an MCMC generator by injecting a dedicated instance of Random.
 */
std::unique_ptr<LinearGeneratorWrap> LinearGeneratorWrap::BuildBubleyDyerGenerator(
            const POSetWrap* poset_wrap,
            std::uint64_t seed) {
    
    auto result = std::make_unique<LinearGeneratorWrap>();
    result->type_ = LEGType::kBubleyDyer;
    
    // Costruisce in via esclusiva il generatore RNG per questa istanza
    auto rnd = std::make_unique<Random>(seed);
    
    // Inietta il POSet e cede l'ownership del RNG al generatore Bubley-Dyer
    result->poset_ = poset_wrap->GetPOSet();
    result->finite_ = false;
    result->generator_ = std::make_unique<LEGBubleyDyer>(result->poset_, std::move(rnd));
    result->generator_->Start(0);
    return result;
}

/**
 * @brief Factory for Ideals-based Generators.
 * @details Determines the algorithm type (Ideals or Binary Variables) by analyzing the POSet.
 */
std::unique_ptr<LinearGeneratorWrap> LinearGeneratorWrap::BuildLEGenerator(
                                                                           const POSetWrap* poset_wrap) {
    
    auto result = std::make_unique<LinearGeneratorWrap>();
    
    if (poset_wrap->GetType() == POSetWrap::PosetType::kBinaryVariable) {
        result->type_ = LEGType::kBinaryVariable;
    } else {
        result->type_ = LEGType::kTreeOfIdeals;
    }

    result->poset_ = poset_wrap->GetPOSet();
    result->finite_ = true;

    result->generator_ = poset_wrap->GetPOSet()->CreateLinearExtensionGenerator();
    result->generator_->Start(0);
    return result;
}

/**
 * @brief Public interface for generic extraction.
 */
std::vector<std::vector<std::string>> LinearGeneratorWrap::GetFromLE(
                                                                     bool from_start,
                                                                     std::optional<std::uint64_t> quante,
                                                                     std::optional<std::uint64_t> output_ogni_in_sec) {
    return GetFromLEInternal(from_start, quante, output_ogni_in_sec);
}

/**
 * @brief Returns the number of available extensions.
 */
std::uint64_t LinearGeneratorWrap::LESize() const {
    return generator_->LeSize();
}

/**
 * @brief Internal implementation of extraction with feedback.
 * @details Uses DisplayMessage to communicate progress status to R.
 * Implements memory optimizations via reserve().
 */
std::vector<std::vector<std::string>> LinearGeneratorWrap::GetFromLEInternal(
                                                                             bool from_start,
                                                                             std::optional<std::uint64_t> quante,
                                                                             std::optional<std::uint64_t> output_ogni_in_sec) {
    if (from_start) {
        // Start con il limite massimo
        generator_->Start(std::numeric_limits<std::uint64_t>::max());
        finite_ = false;
    }
    
    std::vector<std::vector<std::string>> result;
    
    if (finite_) return result;
    
    // Preallocazione se sappiamo esattamente quante ne vogliamo estrarre
    if (quante.has_value() && quante.value() > 0) {
        result.reserve(quante.value());
    }
    
    std::uint64_t le_count = 0;
    
    std::unique_ptr<DisplayMessage> display_message;
    if (!output_ogni_in_sec.has_value()) {
        display_message = std::make_unique<DisplayMessageNull>();
    } else {
        // Passa le_count per reference per l'aggiornamento real-time su R
        display_message = std::make_unique<DisplayMessageLEGetR>(
                                                                 le_count, quante, output_ogni_in_sec.value());
    }
    
    display_message->Start();
    
    for (le_count = 0; !quante.has_value() || le_count < quante.value(); ) {
        display_message->Display();
        
        auto le = generator_->Get();
        std::vector<std::string> le_str(le.size());
        
        for (std::uint64_t r = 0; r < le.size(); ++r) {
            le_str[r] = poset_->GetElementName(le.GetVal(r));
        }
        
        // Move semantics: aggiunge il vettore in O(1) senza allocare nuove stringhe
        result.push_back(std::move(le_str));
        ++le_count;
        // Condizione di uscita se lo spazio è esaurito
        if (!generator_->HasNext()) {
            finite_ = true;
            std::string message = "Reached maximum number of linear extentions: " + std::to_string(le_count) + ".";
            Rprintf("%s\n", message.c_str());
            display_message->UpdateTotal(le_count);
            break;
        }
        generator_->Next();
    }
    
    display_message->Stop();
    
    return result; 
}

/**
 * @brief Bubley-Dyer extraction with iterative calculation based on error.
 * @throws MyException If quantity and error parameters are inconsistent.
 */
std::vector<std::vector<std::string>> LinearGeneratorWrap::GetFromBubleyDyer(bool from_start,
                                                                             std::optional<std::uint64_t> quante,
                                                                             std::optional<double> errore,
                                                                             std::optional<std::uint64_t> output_ogni_in_sec) {
    // 1. Validazione input: o specifichi le quantità o l'errore
    if (!errore.has_value() && !quante.has_value()) {
        throw MyException("The number of le extensions to generate or the error must be provided.");
    }
    
    // 2. Warning se vengono specificati entrambi
    if (errore.has_value() && quante.has_value()) {
        std::string err_str = "The number of le generated will be: " + std::to_string(quante.value()) + ". Parameter error is ignored";
        forward_warning_to_r(err_str); 
        
    }
    // 3. Calcolo dinamico in base all'errore
    else if (!quante.has_value()) {
        // Cast sicuro del generatore di base a LEGBubleyDyer
        auto* bd_generator = static_cast<LEGBubleyDyer*>(generator_.get());
        
        // Calcola il numero di iterazioni basato sull'errore puntato
        auto q = bd_generator->EvaluateNumberOfIteration(errore.value());
        
        if (from_start) {
            quante = q;
        } else {
            auto nle = bd_generator->CurrentNumberOfLe();
            if ((q + 1) < nle) {
                throw MyException("GetFromBubleyDyer: Internal error.");
            }
            quante = (q + 1) - nle;
        }
        
        std::string message = std::to_string(quante.value()) + " linear extentions will be generated.";
        Rprintf("%s\n", message.c_str());
    }
    
    auto result = GetFromLEInternal(from_start, quante.value(), output_ogni_in_sec);
    return result;
}

/**
 * @brief Factory for MRP. Allocates the output matrix based on the POSet dimension.
 */
std::unique_ptr<BubleyDyerMRPGenerator> BubleyDyerMRPGenerator::BuildBubleyDyerMRPGenerator(
                                                                                            const POSetWrap* poset_wrap,
                                                                                            std::uint64_t seed) {
    // Allocazione stack-safe tramite unique_ptr (zero overhead)
    auto result = std::make_unique<BubleyDyerMRPGenerator>();
    
    // Estrae il raw pointer dal wrapper per usarlo come osservatore
    POSet* poset = poset_wrap->GetPOSet();
    result->poset_ = poset;
    
    // Inizializza la matrice MRP quadrata basata sulla dimensione del POSet
    result->mrp_ = std::make_unique<Tensor<double, 2>>(std::array<std::uint64_t, 2>{poset->size(), poset->size()}, 0.0);
    
    // Costruisce in via esclusiva il generatore RNG per questa istanza
    auto rnd = std::make_unique<Random>(seed);
    
    // Istanzia il generatore cedendo l'ownership di 'rnd'
    result->le_generator_ = std::make_unique<LEGBubleyDyer>(result->poset_, std::move(rnd));
    
    // Inizializzazione essenziale dello stato del generatore
    result->le_generator_->Start(0);
    
    return result;
}

/**
 * @brief Factory for Multiple Evaluations.
 * @details Maps the function name strings into the respective concrete FLE classes.
 */
std::unique_ptr<BubleyDyerEvaluationGenerator>
BubleyDyerEvaluationGenerator::BuildBubleyDyerEvaluationGenerator(const POSetWrap* poset_wrap,
                                                                  std::uint64_t seed,
                                                                  const std::vector<std::string>& internal_functions,
                                                                  const std::vector<SEXP>& external_functions)
{
    auto result = std::make_unique<BubleyDyerEvaluationGenerator>();
    
    const std::size_t num_funcs = internal_functions.size();
    
    // Pre-allocazione per evitare riallocazioni
    result->eval_results_.resize(num_funcs);
    result->fles_.resize(num_funcs);
    result->functions_name_.resize(num_funcs);

    // Estrazione raw pointer
    result->poset_ = poset_wrap->GetPOSet();
    
    
    for (std::uint64_t k = 0; k < num_funcs; ++k) {
        const auto& nome_funzione = internal_functions[k];
        auto r_funzione = external_functions[k];
        
        auto it = POSetWrap::kFunctionLinearMapType.find(nome_funzione);
        if (it == POSetWrap::kFunctionLinearMapType.end()) {
            throw MyException("Wrong function name: " + nome_funzione);
        }
        
        std::unique_ptr<FunctionLinearExtension> fle = nullptr;
        std::string function_name = "";

        auto function_type = it->second;
        
        switch (function_type) {
            case POSetWrap::FunctionLinearType::kMutualRankingProbability:
                fle = std::make_unique<FLEMutualRankingProbability>(result->poset_);
                function_name = "MutualRankingProbability";
                break;
            case POSetWrap::FunctionLinearType::kAverageHeight:
                fle = std::make_unique<FLEAverageHeight>(result->poset_);
                function_name = "AverageHeight";
                break;
            case POSetWrap::FunctionLinearType::kSeparationAsymmetricLower:
                fle = std::make_unique<FLESeparationAsymmetricLower>(result->poset_);
                function_name = "asymmetricLower";
                break;
            case POSetWrap::FunctionLinearType::kSeparationAsymmetricUpper:
                fle = std::make_unique<FLESeparationAsymmetricUpper>(result->poset_);
                function_name = "asymmetricUpper";
                break;
            case POSetWrap::FunctionLinearType::kSeparationSymmetric:
                fle = std::make_unique<FLESeparationSymmetric>(result->poset_);
               function_name = "symmetric";
                break;
            case POSetWrap::FunctionLinearType::kRFunction:
                fle = std::make_unique<FLERInterface>(result->poset_, r_funzione);
                function_name = "RFunction";
                break;
            case POSetWrap::FunctionLinearType::kDominance:
            case POSetWrap::FunctionLinearType::kMannWhitneyDominance:
            case POSetWrap::FunctionLinearType::kMannWhitneyInferentialDominance:
                throw MyException("kDominance, kMannWhitneyDominance, kMannWhitneyInferentialDominance: not allowed here");
                break;
        }
        auto shape = fle->Shape();
        std::uint_fast64_t nrow = shape.at(0);
        std::uint_fast64_t ncol = shape.at(1);
        result->eval_results_.at(k) = std::make_unique<Tensor<double, 2>>(std::array<std::uint64_t, 2>{nrow, ncol});
        
        result->fles_.at(k) = std::move(fle);
        result->functions_name_.at(k) = function_name;
    }
    
    auto rnd = std::make_unique<Random>(seed);
    result->le_generator_ = std::make_unique<LEGBubleyDyer>(result->poset_, std::move(rnd));
    result->le_generator_->Start(0);
    return result;
}


