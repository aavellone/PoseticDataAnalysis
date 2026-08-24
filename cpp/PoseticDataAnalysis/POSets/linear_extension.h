/**
 * @file linear_extension.h
 * @brief Header file for the LinearExtension core data structure.
 */

#pragma once

#include <cstdint>
#include <vector>
#include <format>
#include <string>
#include <ranges>

/**
 * @class LinearExtension
 * @brief Represents a linear extension (total order) in a lattice/poset.
 *
 * @details Core data structure for HPC. Optimized to guarantee an O(1)
 * read/write access complexity. It maintains a contiguous memory bidirectional
 * mapping (position <-> element), maximizing L1/L2 cache prefetching.
 * Completely decoupled from textual I/O logic for maximum throughput.
 */
class LinearExtension {
private:
    /// @brief Vector mapping the position index to the element ID.
    std::vector<std::uint64_t> by_position_;
    
    /// @brief Vector mapping the element ID to its physical position.
    std::vector<std::uint64_t> by_element_;
    
public:
    /**
     * @brief Constructor with static pre-allocation.
     *
     * @param size The total number of elements in the linear extension.
     * @details The 'explicit' modifier prevents implicit conversions from integer types.
     * Vectors are pre-allocated using the member initializer list for optimal performance.
     */
    explicit LinearExtension(std::uint64_t size)
    : by_position_(size),
    by_element_(size) {}
    
    /**
     * @brief Returns the total number of elements.
     *
     * @return std::uint64_t The size of the internal vectors.
     * @details Kept in snake_case intentionally (as an exception to Google Style)
     * to maintain interface uniformity with standard STL containers.
     */
    [[nodiscard]] std::uint64_t size() const noexcept {
        return by_position_.size();
    }
    
    /**
     * @brief Retrieves the element present at a specific position.
     *
     * @param index The index of the requested position.
     * @return std::uint64_t The numerical identifier of the element.
     * @details Branchless O(1) read access.
     */
    [[nodiscard]] std::uint64_t GetVal(std::uint64_t index) const noexcept {
        return by_position_[index];
    }
    
    /**
     * @brief Retrieves the physical position occupied by a specific element.
     *
     * @param element The numerical identifier of the target element.
     * @return std::uint64_t The index where the element is located.
     * @details Branchless O(1) read access.
     */
    [[nodiscard]] std::uint64_t GetPos(std::uint64_t element) const noexcept {
        return by_element_[element];
    }
    
    /**
     * @brief Updates the linear extension by registering an element at a given index.
     *
     * @param index The target position to overwrite.
     * @param value The element value to store.
     * @details Updates both internal mapping structures simultaneously in O(1).
     */
    void Set(std::uint64_t index, std::uint64_t value) noexcept {
        by_position_[index] = value;
        by_element_[value] = index;
    }
    
    /**
     * @brief Genera una rappresentazione testuale della sequenza.
     * @return Una stringa nel formato "[e1, e2, e3, ...]"
     */
    [[nodiscard]] std::string ToString() const {
        std::string result = "[";
        // Usiamo un range view per iterare sulla posizione
        for (size_t i = 0; i < by_position_.size(); ++i) {
            // std::format gestisce la formattazione in modo efficiente
            result += std::format("{}{}", by_position_[i],
                                  (i == by_position_.size() - 1 ? "" : ", "));
        }
        result += "]";
        return result;
    }
};
