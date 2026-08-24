/**
 * @file linear_extension_generator_bubley_dyer.h
 * @brief Header file for the Bubley-Dyer linear extension generator.
 */

#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <variant>

#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <string>

#include "linear_extension.h"
#include "linear_extension_generator.h"

// Forward declarations to reduce compilation dependencies
class POSet;
class Random;

/**
 * @class LEGBubleyDyer
 * @brief Linear extension generator based on the Bubley-Dyer algorithm.
 *
 * @details Uses a Markov Chain Monte Carlo (MCMC) approach to quickly and
 * almost uniformly generate linear extensions of a Partially Ordered Set (POSet).
 * * @par References:
 * - Bubley, R., & Dyer, M. (1999). "Faster random generation of linear extensions."
 * Discrete Mathematics, 201(1-3), 81-88.
 */
class LEGBubleyDyer final : public LinearExtensionGenerator {
public:
    /**
     * @brief Constructs the LEGBubleyDyer generator.
     *
     * @param poset Pointer to the Partially Ordered Set (POSet) structure.
     * @param rnd Unique pointer to the Random number generator instance.
     */
    LEGBubleyDyer(POSet* poset, std::unique_ptr<Random> rnd);

    /**
     * @brief Constructs the generator with a caller-supplied starting point.
     *
     * @details Initializes initial_le_ from @p initial_le instead of
     * poset->FirstLE(): consente a catene MCMC diverse di partire da estensioni
     * lineari diverse (utile per la valutazione parallela multi-catena).
     *
     * @param poset Pointer to the Partially Ordered Set (POSet) structure.
     * @param rnd Unique pointer to the Random number generator instance.
     * @param initial_le Starting linear extension (moved in).
     * @throws MyException se @p initial_le non è una permutazione coerente di
     *         dimensione poset->size() o viola l'ordine parziale del poset.
     */
    LEGBubleyDyer(POSet* poset, std::unique_ptr<Random> rnd, LinearExtension initial_le);


    /**
     * @brief Default virtual destructor.
     */
    ~LEGBubleyDyer() override = default;
    
    // ABILITA IL MOVIMENTO (Essenziale per std::variant)
    LEGBubleyDyer(LEGBubleyDyer&&) noexcept = default;
    LEGBubleyDyer& operator=(LEGBubleyDyer&&) noexcept = default;
    
    // DISABILITA LA COPIA (Previene bug di performance e collisioni su unique_ptr)
    LEGBubleyDyer(const LEGBubleyDyer&) = delete;
    LEGBubleyDyer& operator=(const LEGBubleyDyer&) = delete;
    
    /**
     * @brief Initializes the generator and the MCMC state.
     *
     * @param quante The maximum number of linear extensions to generate.
     */
    void Start(std::uint64_t quante = 0) override;
    
    /**
     * @brief Advances to the next linear extension via the Markov Chain transition.
     */
    void Next() override;
    
    /**
     * @brief Checks if there are more linear extensions to generate.
     *
     * @return true if the current count is less than the maximum allowed.
     */
    [[nodiscard]] bool HasNext() noexcept override {
        // Uses current_number_le_ from the base LinearExtensionGenerator class
        return current_number_le_ < max_number_le_;
    }
    
    /**
     * @brief Returns the maximum number of linear extensions to generate.
     *
     * @return The maximum limit of extensions set for this generator.
     */
    [[nodiscard]] std::uint64_t NumberOfLe() const noexcept override {
        return max_number_le_;
    }
    
    /**
     * @brief Evaluates the required number of MCMC iterations for a given error bound.
     *
     * @param error The acceptable error margin for the uniform distribution.
     * @return The number of MCMC steps/iterations needed.
     */
    [[nodiscard]] std::uint64_t EvaluateNumberOfIteration(double error) const;
    
    /**
     * @brief Updates internal generation counters based on given parameters.
     *
     * @param quante Optional target number of linear extensions to generate.
     * @param error Optional acceptable error bound for the MCMC process.
     * @return true if the counters were successfully updated, false otherwise.
     */
    [[nodiscard]] bool UpdateCounters(std::optional<std::uint64_t> quante,
                                      std::optional<double> error);
    
private:
    POSet* poset_;                          ///< Pointer to the associated POSet.
    std::unique_ptr<Random> rnd_;           ///< Random number generator.
    
    bool to_update_ = false;                ///< Flag indicating if an update is pending.
    bool is_switched_ = false;              ///< Flag indicating if a state switch occurred.
    std::uint64_t position_to_update_ = 0;  ///< Index of the position to be updated.
    
    /// Max number of linear extensions (defaults to max uint64_t).
    std::uint64_t max_number_le_ = std::numeric_limits<std::uint64_t>::max();
    
    LinearExtension initial_le_;            ///< Buffer storing the initial linear extension.
};
