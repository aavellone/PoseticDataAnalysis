/**
 * @file leg_lexicographic_from_linear.h
 * @brief Header for the lexicographical linear extension generator based on
 * priority groups.
 */

#pragma once

#include <cstdint>
#include <vector>

#include "linear_extension_generator.h"

/**
 * @class LEGLexicographicFromLinear
 * @brief Generates linear extensions by varying the lexicographical priority of
 * groups.
 * @details Specifically designed for product POSets of linear chains.
 */
class LEGLexicographicFromLinear final : public LinearExtensionGenerator {
public:
    /**
     * @brief Constructor.
     * @param group_sizes Vector containing the size of each individual group.
     */
    explicit LEGLexicographicFromLinear(std::vector<std::uint64_t> group_sizes);
    
    ~LEGLexicographicFromLinear() override = default;
    
    // LinearExtensionGenerator Interface
    void Start(std::uint64_t unused_id = 0) override;
    void Next() override;
    [[nodiscard]] bool HasNext() noexcept override;
    [[nodiscard]] std::uint64_t NumberOfLe() const noexcept override;
    
private:
    /**
     * @brief Calculates the total number of elements (sum of group sizes).
     * @param sizes Vector of group dimensions.
     * @return Total number of elements as uint64_t.
     */
    [[nodiscard]] static std::uint64_t CalculateTotalSize(
                                                          const std::vector<std::uint64_t> &sizes) noexcept;
    
    /**
     * @brief Updates the current_linear_extension_ buffer based on current
     * priorities.
     */
    void RefreshData() noexcept;
    
    // Variabili membro: convenzione Google prevede snake_case con underscore
    // finale
    std::vector<std::uint64_t> group_sizes_; ///< Sizes of each group.
    std::vector<std::uint64_t> priority_; ///< Current permutation of group priorities.
    std::vector<std::uint64_t> strides_;          ///< Multipliers to flatten multi-dimensional coordinates.
    std::uint64_t total_permutations_; ///< Total permutations (n_groups!).
    std::vector<std::uint64_t> coords_;         ///< Buffer persistente per le coordinate.
    std::vector<std::uint64_t> wrap_offsets_;   ///< Pre-calcolo di (size-1) * stride.
    std::uint64_t total_elements_;
};
