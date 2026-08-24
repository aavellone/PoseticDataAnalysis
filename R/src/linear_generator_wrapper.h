/**
 * @file linear_generator_wrapper.h
 * @brief HPC wrapper for Linear Extension generators.
 *
 * @details Manages the lifecycle of generators and the production of linear extensions.
 */
#pragma once

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "function_linear_extension.h"
#include "linear_extension_generator.h"
#include "linear_extension_generator_bubley_dyer.h"
#include "tensor.h"
#include "poset_wrapper.h"

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif
#include <Rinternals.h>
#include <R.h>

/**
 * @class LinearGeneratorWrap
 * @brief Main manager for the execution of linear extension generation algorithms.
 *
 * @details This class acts as a unified interface to the different generation engines
 * (Bubley-Dyer, Ideals, etc.). It implements move semantics to ensure that the
 * generator's ownership is unique and non-duplicable.
 */
class LinearGeneratorWrap {
public:
    /**
     * @enum LEGType
     * @brief Available generation algorithm types.
     */
    enum class LEGType {
        kBubleyDyer,        ///< Algorithm based on the Bubley-Dyer Markov Chain.
        kTreeOfIdeals,      ///< Exact algorithm based on the tree of ideals.
        kBinaryVariable,    ///< Specific optimization for binary variables.
        kFromLinearPosets,  ///< Generation starting from already linear posets.
        kCount              ///< Sentinel for the number of types.
    };
    
    /** @brief Default constructor. */
    LinearGeneratorWrap() = default;
    /** @brief Default destructor. */
    ~LinearGeneratorWrap() = default;
    
    /** @name Rule of Five
     * Copying is disabled to preserve the integrity of the unique pointer.
     */
    ///@{
    LinearGeneratorWrap(const LinearGeneratorWrap&) = delete;
    LinearGeneratorWrap& operator=(const LinearGeneratorWrap&) = delete;
    LinearGeneratorWrap(LinearGeneratorWrap&&) noexcept = default;
    LinearGeneratorWrap& operator=(LinearGeneratorWrap&&) noexcept = default;
    ///@}
    
    /**
     * @brief Factory method to instantiate a Bubley-Dyer generator.
     * @param poset Pointer to the source POSet wrapper.
     * @param seed Seed for random number generator initialization.
     * @return std::unique_ptr<LinearGeneratorWrap> Configured instance of the wrapper.
     */
    [[nodiscard]] static std::unique_ptr<LinearGeneratorWrap> BuildBubleyDyerGenerator(
                                                                         const POSetWrap* poset,
                                                                         std::uint64_t seed);
    
    /**
     * @brief Factory method to instantiate an exact generator based on the POSet structure.
     * @param poset Pointer to the source POSet wrapper.
     * @return std::unique_ptr<LinearGeneratorWrap> Configured instance of the wrapper.
     */
    [[nodiscard]] static std::unique_ptr<LinearGeneratorWrap> BuildLEGenerator(
                                                                 const POSetWrap* poset);
    
    /**
     * @brief Extracts a set of linear extensions.
     * @param from_start If true, resets the iteration from the beginning of the generator's state.
     * @param quante Maximum number of extensions to produce.
     * @param output_ogni_in_sec Interval in seconds for feedback to the R interface.
     * @return std::vector<std::vector<std::string>> Matrix of strings containing the extensions.
     */
    [[nodiscard]] std::vector<std::vector<std::string>> GetFromLE(
                                                                  bool from_start,
                                                                  std::optional<std::uint64_t> quante,
                                                                  std::optional<std::uint64_t> output_ogni_in_sec);
    /**
     * @brief Performs sampling via Bubley-Dyer with statistical criteria.
     * @param from_start If true, performs initial setup procedures (burn-in).
     * @param quante Number of requested samples.
     * @param errore Pointer to the desired error value (optional).
     * @param output_ogni_in_sec Frequency of video feedback.
     * @return std::vector<std::vector<std::string>> Generated samples.
     */
    [[nodiscard]] std::vector<std::vector<std::string>> GetFromBubleyDyer(
                                                                          bool from_start,
                                                                          std::optional<std::uint64_t> quante,
                                                                          std::optional<double> errore,
                                                                          std::optional<std::uint64_t> output_ogni_in_sec);
    
    /**
     * @brief Returns the total number of possible linear extensions.
     * @return std::uint64_t Total count (or estimated for sampling).
     */
    [[nodiscard]] std::uint64_t LESize() const;
    
    /**
     * @brief Accesses the internal generator.
     * @return LinearExtensionGenerator* Raw pointer to the managed object (observer).
     */
    [[nodiscard]] LinearExtensionGenerator* GetGenerator() const noexcept { return generator_.get(); }
private:
    /**
     * @brief Core implementation of the data extraction logic.
     * @details Manages memory preallocation and the element retrieval loop.
     */
    [[nodiscard]] std::vector<std::vector<std::string>> GetFromLEInternal(
                                                                          bool from_start,
                                                                          std::optional<std::uint64_t> quante,
                                                                          std::optional<std::uint64_t> output_ogni_in_sec);
    
    std::unique_ptr<LinearExtensionGenerator> generator_; ///< Underlying calculation engine.
    POSet* poset_{nullptr};                               ///< Reference to the POSet (observer).
    LEGType type_;                                        ///< Algorithm currently in use.
    bool finite_{false};                                  ///< Flag indicating if generation is finished.
};

/**
 * @class BubleyDyerMRPGenerator
 * @brief Specialized wrapper for calculating Mutual Ranking Probabilities (MRP).
 * @details Integrates the generator and the results matrix to optimize statistical calculations.
 */
class BubleyDyerMRPGenerator {
public:
    std::unique_ptr<LEGBubleyDyer> le_generator_{nullptr}; ///< MCMC Generator.
    std::unique_ptr<Tensor<double, 2>> mrp_{nullptr};                  ///< Square matrix of MRP results.
    bool used_{false};                                               ///< Usage state of the generator.
    POSet* poset_{nullptr};                                          ///< Associated POSet.
    
    BubleyDyerMRPGenerator() = default;
    ~BubleyDyerMRPGenerator() = default;
    
    /**
     * @brief Constructs an MRP generator with initialization of the results matrix.
     * @param poset POSet wrapper.
     * @param seed Seed for the random generator.
     * @return std::unique_ptr<BubleyDyerMRPGenerator> Pointer to the created instance.
     */
    [[nodiscard]] static std::unique_ptr<BubleyDyerMRPGenerator> BuildBubleyDyerMRPGenerator(
                                                                                             const POSetWrap* poset,
                                                                                             std::uint64_t seed);
};

/**
 * @class BubleyDyerEvaluationGenerator
 * @brief Wrapper for multiple evaluation of metrics on linear extensions.
 * @details Allows applying multiple evaluation functions (FLE) simultaneously
 * during the same MCMC generation cycle, maximizing efficiency.
 */
class BubleyDyerEvaluationGenerator {
public:
    std::unique_ptr<LEGBubleyDyer> le_generator_{nullptr}; ///< MCMC Generator.
    std::vector<std::unique_ptr<Tensor<double, 2>>> eval_results_;    ///< Vector of matrices for results.
    std::vector<std::unique_ptr<FunctionLinearExtension>> fles_;    ///< Set of functions to apply.
    std::vector<std::string> functions_name_;                       ///< Identifying names of the functions.
    bool used_{false};                                              ///< Usage state.
    POSet* poset_{nullptr};                                         ///< Associated POSet.
    
    BubleyDyerEvaluationGenerator() = default;
    ~BubleyDyerEvaluationGenerator() = default;
    
    /**
     * @brief Constructs a multiple evaluator by mapping internal and external (R) functions.
     * @param poset POSet wrapper.
     * @param seed Random seed.
     * @param internal_functions Names of the C++ metrics.
     * @param external_functions Pointers to R functions (SEXP) for custom evaluations.
     * @return std::unique_ptr<BubleyDyerEvaluationGenerator> Pointer to the instance.
     */
    [[nodiscard]] static std::unique_ptr<BubleyDyerEvaluationGenerator>
    BuildBubleyDyerEvaluationGenerator(const POSetWrap* poset,
                                       std::uint64_t seed,
                                       const std::vector<std::string>& internal_functions,
                                       const std::vector<SEXP>& external_functions);
};
