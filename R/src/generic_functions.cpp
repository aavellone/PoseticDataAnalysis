
// MSVC: explicit standard includes (not pulled in transitively).
#include <cstdint>
/**
 * @file genericFunctions.cpp
 * @brief Implementation of binary relation operations.
 *
 * @details Core mathematical implementations optimized following CPU-friendly
 * paradigms (cache-locality, branch prediction minimization, in-place updates).
 */

#include "generic_functions.h"

namespace generic {
    
    bool IsReflexive(std::uint64_t num_elements, const Tensor<std::uint8_t, 2>& adj) {
        for (std::uint64_t i = 0; i < num_elements; ++i) {
            if (adj(i, i) == 0) {
                return false;
            }
        }
        return true;
    }
    
    bool IsSymmetric(std::uint64_t num_elements, const Tensor<std::uint8_t, 2>& adj) {
        for (std::uint64_t i = 0; i < num_elements; ++i) {
            for (std::uint64_t j = i + 1; j < num_elements; ++j) {
                if (adj(i, j) != adj(j, i)) {
                    return false;
                }
            }
        }
        return true;
    }
    
    bool IsAntisymmetric(std::uint64_t num_elements, const Tensor<std::uint8_t, 2>& adj) {
        for (std::uint64_t i = 0; i < num_elements; ++i) {
            for (std::uint64_t j = i + 1; j < num_elements; ++j) {
                if (adj(i, j) == 1 && adj(j, i) == 1) {
                    return false;
                }
            }
        }
        return true;
    }
    
    bool IsTransitive(std::uint64_t num_elements, const Tensor<std::uint8_t, 2>& adj) {
        for (std::uint64_t i = 0; i < num_elements; ++i) {
            for (std::uint64_t j = 0; j < num_elements; ++j) {
                if (adj(i, j) == 1) {
                    for (std::uint64_t k = 0; k < num_elements; ++k) {
                        if (adj(j, k) == 1 && adj(i, k) == 0) {
                            return false;
                        }
                    }
                }
            }
        }
        return true;
    }
    
    bool IsPreorder(std::uint64_t num_elements, const Tensor<std::uint8_t, 2>& adj) {
        return IsReflexive(num_elements, adj) && IsTransitive(num_elements, adj);
    }
    
    bool IsPartialOrder(std::uint64_t num_elements, const Tensor<std::uint8_t, 2>& adj) {
        return IsReflexive(num_elements, adj) &&
        IsAntisymmetric(num_elements, adj) &&
        IsTransitive(num_elements, adj);
    }
    
    void ReflexiveClosureInPlace(std::uint64_t num_elements, Tensor<std::uint8_t, 2>& adj) {
        for (std::uint64_t i = 0; i < num_elements; ++i) {
            adj(i, i) = 1;
        }
    }
    
    void TransitiveClosureInPlace(std::uint64_t num_elements, Tensor<std::uint8_t, 2>& adj) {
        // Optimized Floyd-Warshall with k on the outer loop
        for (std::uint64_t k = 0; k < num_elements; ++k) {
            for (std::uint64_t i = 0; i < num_elements; ++i) {
                if (adj(i, k) == 1) {
                    for (std::uint64_t j = 0; j < num_elements; ++j) {
                        if (adj(k, j) == 1) {
                            adj(i, j) = 1;
                        }
                    }
                }
            }
        }
    }
    
}  // namespace generic
