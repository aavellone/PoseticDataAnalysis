/**
 * @file genericFunctions.h
 * @brief High-Performance Computing (HPC) operations for binary relations.
 *
 * @details Provides mathematical property checks (reflexivity, transitivity, etc.)
 * and closure operations on binary relations. Optimized for HPC environments
 * using contiguous byte matrices (Tensor<std::uint8_t, 2>) to maximize CPU
 * cache hits and eliminate pointer chasing overhead.
 *
 * @author Alessandro Avellone
 * @version 3.0
 * @date 2026
 */

#pragma once

#include "tensor.h"

#include <cstdint>


namespace generic {
    
    /**
     * @brief Checks if a binary relation is reflexive.
     *
     * @param num_elements The total number of elements in the relation domain.
     * @param adj The adjacency matrix representing the relation.
     * @return true if every element is related to itself, false otherwise.
     * @note Complexity: O(N) where N is num_elements.
     */
    bool IsReflexive(std::uint64_t num_elements, const Tensor<std::uint8_t, 2>& adj);
    
    /**
     * @brief Checks if a binary relation is symmetric.
     *
     * @param num_elements The total number of elements in the relation domain.
     * @param adj The adjacency matrix representing the relation.
     * @return true if for every pair (a, b), the pair (b, a) exists; false otherwise.
     * @note Complexity: O(N^2). Explores only the upper triangle.
     */
    bool IsSymmetric(std::uint64_t num_elements, const Tensor<std::uint8_t, 2>& adj);
    
    /**
     * @brief Checks if a binary relation is antisymmetric.
     *
     * @param num_elements The total number of elements in the relation domain.
     * @param adj The adjacency matrix representing the relation.
     * @return true if there are no mutual relations between distinct elements.
     * @note Complexity: O(N^2). Explores only the upper triangle.
     */
    bool IsAntisymmetric(std::uint64_t num_elements, const Tensor<std::uint8_t, 2>& adj);
    
    /**
     * @brief Checks if a binary relation is transitive.
     *
     * @param num_elements The total number of elements in the relation domain.
     * @param adj The adjacency matrix representing the relation.
     * @return true if (a, b) and (b, c) always imply (a, c); false otherwise.
     * @note Complexity: O(N^3). Optimized loop order for L1 cache prefetching.
     */
    bool IsTransitive(std::uint64_t num_elements, const Tensor<std::uint8_t, 2>& adj);
    
    /**
     * @brief Checks if a binary relation is a preorder.
     *
     * A preorder is reflexive and transitive.
     *
     * @param num_elements The total number of elements in the relation domain.
     * @param adj The adjacency matrix representing the relation.
     * @return true if the relation is a preorder, false otherwise.
     */
    bool IsPreorder(std::uint64_t num_elements, const Tensor<std::uint8_t, 2>& adj);
    
    /**
     * @brief Checks if a binary relation is a partial order.
     *
     * A partial order is reflexive, antisymmetric, and transitive.
     *
     * @param num_elements The total number of elements in the relation domain.
     * @param adj The adjacency matrix representing the relation.
     * @return true if the relation is a partial order, false otherwise.
     */
    bool IsPartialOrder(std::uint64_t num_elements, const Tensor<std::uint8_t, 2>& adj);
    
    /**
     * @brief Computes the reflexive closure of a relation in-place.
     *
     * Mutates the input adjacency matrix to ensure every element relates to itself.
     *
     * @param num_elements The total number of elements in the relation domain.
     * @param adj The adjacency matrix representing the relation (modified).
     * @note Zero memory allocations. Complexity: O(N).
     */
    void ReflexiveClosureInPlace(std::uint64_t num_elements, Tensor<std::uint8_t, 2>& adj);
    
    /**
     * @brief Computes the transitive closure of a relation in-place.
     *
     * Utilizes a boolean-optimized version of the Floyd-Warshall algorithm.
     * Mutates the input adjacency matrix.
     *
     * @param num_elements The total number of elements in the relation domain.
     * @param adj The adjacency matrix representing the relation (modified).
     * @note Zero memory allocations. Complexity: O(N^3).
     */
    void TransitiveClosureInPlace(std::uint64_t num_elements, Tensor<std::uint8_t, 2>& adj);
    
}  // namespace generic
