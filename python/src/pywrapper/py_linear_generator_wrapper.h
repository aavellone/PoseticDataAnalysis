/**
 * @file py_linear_generator_wrapper.h
 * @brief Wrapper HPC per i generatori di estensioni lineari (versione Python).
 *
 * @details Copia R-free di linear_generator_wrapper.h: gestisce il ciclo di
 * vita dei generatori e la produzione di estensioni lineari. Rispetto alla
 * versione R, la factory di valutazione NON accetta funzioni-FLE scritte in R
 * (solo le metriche C++ interne per nome).
 */
#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "function_linear_extension.h"
#include "linear_extension_generator.h"
#include "linear_extension_generator_bubley_dyer.h"
#include "tensor.h"
#include "poset_wrapper.h"

/**
 * @class LinearGeneratorWrap
 * @brief Manager per l'esecuzione degli algoritmi di generazione di estensioni lineari.
 */
class LinearGeneratorWrap {
  public:
    enum class LEGType {
        kBubleyDyer,
        kTreeOfIdeals,
        kBinaryVariable,
        kFromLinearPosets,
        kCount
    };

    LinearGeneratorWrap() = default;
    ~LinearGeneratorWrap() = default;

    LinearGeneratorWrap(const LinearGeneratorWrap&) = delete;
    LinearGeneratorWrap& operator=(const LinearGeneratorWrap&) = delete;
    LinearGeneratorWrap(LinearGeneratorWrap&&) noexcept = default;
    LinearGeneratorWrap& operator=(LinearGeneratorWrap&&) noexcept = default;

    [[nodiscard]] static std::unique_ptr<LinearGeneratorWrap> BuildBubleyDyerGenerator(
        const POSetWrap* poset, std::uint64_t seed);

    [[nodiscard]] static std::unique_ptr<LinearGeneratorWrap> BuildLEGenerator(
        const POSetWrap* poset);

    [[nodiscard]] std::vector<std::vector<std::string>> GetFromLE(
        bool from_start, std::optional<std::uint64_t> quante,
        std::optional<std::uint64_t> output_ogni_in_sec);

    [[nodiscard]] std::vector<std::vector<std::string>> GetFromBubleyDyer(
        bool from_start, std::optional<std::uint64_t> quante,
        std::optional<double> errore, std::optional<std::uint64_t> output_ogni_in_sec);

    [[nodiscard]] std::uint64_t LESize() const;

    [[nodiscard]] LinearExtensionGenerator* GetGenerator() const noexcept {
        return generator_.get();
    }

  private:
    [[nodiscard]] std::vector<std::vector<std::string>> GetFromLEInternal(
        bool from_start, std::optional<std::uint64_t> quante,
        std::optional<std::uint64_t> output_ogni_in_sec);

    std::unique_ptr<LinearExtensionGenerator> generator_;
    POSet* poset_{nullptr};
    LEGType type_;
    bool finite_{false};
};

/**
 * @class BubleyDyerMRPGenerator
 * @brief Wrapper specializzato per il calcolo delle Mutual Ranking Probabilities.
 */
class BubleyDyerMRPGenerator {
  public:
    std::unique_ptr<LEGBubleyDyer> le_generator_{nullptr};
    std::unique_ptr<Tensor<double, 2>> mrp_{nullptr};
    bool used_{false};
    POSet* poset_{nullptr};

    BubleyDyerMRPGenerator() = default;
    ~BubleyDyerMRPGenerator() = default;

    [[nodiscard]] static std::unique_ptr<BubleyDyerMRPGenerator> BuildBubleyDyerMRPGenerator(
        const POSetWrap* poset, std::uint64_t seed);
};

/**
 * @class BubleyDyerEvaluationGenerator
 * @brief Valutazione multipla di metriche sulle estensioni lineari (MCMC).
 * @details Versione Python: solo metriche C++ interne (nessuna funzione R).
 */
class BubleyDyerEvaluationGenerator {
  public:
    std::unique_ptr<LEGBubleyDyer> le_generator_{nullptr};
    std::vector<std::unique_ptr<Tensor<double, 2>>> eval_results_;
    std::vector<std::unique_ptr<FunctionLinearExtension>> fles_;
    std::vector<std::string> functions_name_;
    bool used_{false};
    POSet* poset_{nullptr};

    BubleyDyerEvaluationGenerator() = default;
    ~BubleyDyerEvaluationGenerator() = default;

    [[nodiscard]] static std::unique_ptr<BubleyDyerEvaluationGenerator>
    BuildBubleyDyerEvaluationGenerator(const POSetWrap* poset, std::uint64_t seed,
                                       const std::vector<std::string>& internal_functions);
};
