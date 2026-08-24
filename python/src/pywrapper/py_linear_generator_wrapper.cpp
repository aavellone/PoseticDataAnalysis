/**
 * @file py_linear_generator_wrapper.cpp
 * @brief Implementazione R-free dei wrapper dei generatori (versione Python).
 *
 * @details Copia adattata di linear_generator_wrapper.cpp:
 *  - nessun header R; il progresso e i messaggi passano da py_display.h;
 *  - la factory di valutazione mappa solo le metriche C++ interne per nome.
 */

#include "py_linear_generator_wrapper.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "py_display.h"

#include "function_linear_extension.h"
#include "function_linear_extension_average_height.h"
#include "function_linear_extension_mutual_ranking_probability.h"
#include "function_linear_extension_separation_asymmetric_lower.h"
#include "function_linear_extension_separation_asymmetric_upper.h"
#include "function_linear_extension_separation_symmetric.h"
#include "linear_extension_generator.h"
#include "my_exception.h"
#include "poset.h"
#include "poset_wrapper.h"
#include "random.h"

std::unique_ptr<LinearGeneratorWrap> LinearGeneratorWrap::BuildBubleyDyerGenerator(
    const POSetWrap* poset_wrap, std::uint64_t seed) {
    auto result = std::make_unique<LinearGeneratorWrap>();
    result->type_ = LEGType::kBubleyDyer;

    auto rnd = std::make_unique<Random>(seed);
    result->poset_ = poset_wrap->GetPOSet();
    result->finite_ = false;
    result->generator_ = std::make_unique<LEGBubleyDyer>(result->poset_, std::move(rnd));
    result->generator_->Start(0);
    return result;
}

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

std::vector<std::vector<std::string>> LinearGeneratorWrap::GetFromLE(
    bool from_start, std::optional<std::uint64_t> quante,
    std::optional<std::uint64_t> output_ogni_in_sec) {
    return GetFromLEInternal(from_start, quante, output_ogni_in_sec);
}

std::uint64_t LinearGeneratorWrap::LESize() const {
    return generator_->LeSize();
}

std::vector<std::vector<std::string>> LinearGeneratorWrap::GetFromLEInternal(
    bool from_start, std::optional<std::uint64_t> quante,
    std::optional<std::uint64_t> output_ogni_in_sec) {
    if (from_start) {
        generator_->Start(std::numeric_limits<std::uint64_t>::max());
        finite_ = false;
    }

    std::vector<std::vector<std::string>> result;
    if (finite_) {
        return result;
    }

    if (quante.has_value() && quante.value() > 0) {
        result.reserve(quante.value());
    }

    std::uint64_t le_count = 0;

    std::unique_ptr<DisplayMessage> display_message;
    if (!output_ogni_in_sec.has_value()) {
        display_message = std::make_unique<DisplayMessageNull>();
    } else {
        display_message = std::make_unique<DisplayMessageStdout>(
            le_count, quante, output_ogni_in_sec.value(), "generated");
    }

    display_message->Start();

    for (le_count = 0; !quante.has_value() || le_count < quante.value();) {
        display_message->Display();

        auto le = generator_->Get();
        std::vector<std::string> le_str(le.size());
        for (std::uint64_t r = 0; r < le.size(); ++r) {
            le_str[r] = poset_->GetElementName(le.GetVal(r));
        }
        result.push_back(std::move(le_str));
        ++le_count;

        if (!generator_->HasNext()) {
            finite_ = true;
            PyDisplay::Line("Reached maximum number of linear extensions: " +
                            std::to_string(le_count) + ".");
            display_message->UpdateTotal(le_count);
            break;
        }
        generator_->Next();
    }

    display_message->Stop();
    return result;
}

std::vector<std::vector<std::string>> LinearGeneratorWrap::GetFromBubleyDyer(
    bool from_start, std::optional<std::uint64_t> quante, std::optional<double> errore,
    std::optional<std::uint64_t> output_ogni_in_sec) {
    if (!errore.has_value() && !quante.has_value()) {
        throw MyException("The number of le extensions to generate or the error must be provided.");
    }

    if (errore.has_value() && quante.has_value()) {
        PyDisplay::Warning("The number of le generated will be: " +
                           std::to_string(quante.value()) + ". Parameter error is ignored");
    } else if (!quante.has_value()) {
        auto* bd_generator = static_cast<LEGBubleyDyer*>(generator_.get());
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
        PyDisplay::Line(std::to_string(quante.value()) + " linear extensions will be generated.");
    }

    return GetFromLEInternal(from_start, quante.value(), output_ogni_in_sec);
}

std::unique_ptr<BubleyDyerMRPGenerator> BubleyDyerMRPGenerator::BuildBubleyDyerMRPGenerator(
    const POSetWrap* poset_wrap, std::uint64_t seed) {
    auto result = std::make_unique<BubleyDyerMRPGenerator>();

    POSet* poset = poset_wrap->GetPOSet();
    result->poset_ = poset;
    result->mrp_ = std::make_unique<Tensor<double, 2>>(
        std::array<std::uint64_t, 2>{poset->size(), poset->size()}, 0.0);

    auto rnd = std::make_unique<Random>(seed);
    result->le_generator_ = std::make_unique<LEGBubleyDyer>(result->poset_, std::move(rnd));
    result->le_generator_->Start(0);
    return result;
}

std::unique_ptr<BubleyDyerEvaluationGenerator>
BubleyDyerEvaluationGenerator::BuildBubleyDyerEvaluationGenerator(
    const POSetWrap* poset_wrap, std::uint64_t seed,
    const std::vector<std::string>& internal_functions) {
    auto result = std::make_unique<BubleyDyerEvaluationGenerator>();

    const std::size_t num_funcs = internal_functions.size();
    result->eval_results_.resize(num_funcs);
    result->fles_.resize(num_funcs);
    result->functions_name_.resize(num_funcs);
    result->poset_ = poset_wrap->GetPOSet();

    for (std::uint64_t k = 0; k < num_funcs; ++k) {
        const auto& nome_funzione = internal_functions[k];

        auto it = POSetWrap::kFunctionLinearMapType.find(nome_funzione);
        if (it == POSetWrap::kFunctionLinearMapType.end()) {
            throw MyException("Wrong function name: " + nome_funzione);
        }

        std::unique_ptr<FunctionLinearExtension> fle = nullptr;
        std::string function_name = "";

        switch (it->second) {
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
                throw MyException(
                    "Custom (RFunction) metrics are not supported in the Python "
                    "binding; use the built-in metric names.");
            case POSetWrap::FunctionLinearType::kDominance:
            case POSetWrap::FunctionLinearType::kMannWhitneyDominance:
            case POSetWrap::FunctionLinearType::kMannWhitneyInferentialDominance:
                throw MyException(
                    "kDominance, kMannWhitneyDominance, kMannWhitneyInferentialDominance: "
                    "not allowed here");
        }

        auto shape = fle->Shape();
        std::uint_fast64_t nrow = shape.at(0);
        std::uint_fast64_t ncol = shape.at(1);
        result->eval_results_.at(k) =
            std::make_unique<Tensor<double, 2>>(std::array<std::uint64_t, 2>{nrow, ncol});
        result->fles_.at(k) = std::move(fle);
        result->functions_name_.at(k) = function_name;
    }

    auto rnd = std::make_unique<Random>(seed);
    result->le_generator_ = std::make_unique<LEGBubleyDyer>(result->poset_, std::move(rnd));
    result->le_generator_->Start(0);
    return result;
}
