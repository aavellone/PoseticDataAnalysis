#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <cstdint>

#include "tensor.h"
#include "my_exception.h"

#include <vector>
#include <algorithm>

// =========================================================================
// 1. MIN-MAX TRANSITIVE CLOSURE (Template Implementation)
// =========================================================================

/**
 * @brief Computes the min-max transitive closure of a fuzzy relation.
 * @param relation A square matrix representing the fuzzy relation.
 * @return Tensor<double, 2> The transitive closure of the input relation.
 */
inline Tensor<double, 2> MinMaxTransitiveClosure(Tensor<double, 2> relation) {
    const std::uint64_t n = relation.Extent(0);
    
    if (relation.Extent(1) != n) {
        throw MyException("The matrix does not represent a fuzzy relation (must be square).");
    }
    
    double* const data = relation.data();
    
    for (std::uint64_t i = 0; i < n; ++i) {
        const double* const row_i = data + (i * n);
        
        for (std::uint64_t j = 0; j < n; ++j) {
            if (j == i) continue;
            
            double* const row_j = data + (j * n);
            const double r_ji = row_j[i];
            
            // 4. Vettorizzazione SIMD Cross-Platform (Compatibile con GCC, Clang e MSVC)
#if defined(__clang__)
            // Soluzione nativa per Clang (macOS / Xcode)
#pragma clang loop vectorize(enable)
#elif defined(__GNUC__)
            // Soluzione nativa per GCC (Linux / Rtools su Windows)
#pragma GCC ivdep
#elif defined(_MSC_VER)
            // Soluzione nativa per il compilatore Microsoft
#pragma loop(ivdep)
#endif
            for (std::uint64_t k = 0; k < n; ++k) {
                row_j[k] = std::max(row_j[k], std::min(r_ji, row_i[k]));
            }
        }
    }
    return relation;
}

// =========================================================================
// 2. FIRST ORDER DOMINANCE ANALYSIS (Declarations)
// =========================================================================

/**
 * @brief Data structure containing the results of the First Order Dominance (FOD) analysis.
 */
struct FODAnalysis {
    Tensor<double, 2> fod_closed;      ///< Transitively closed matrix (R: FOD.CLOSED)
    Tensor<double, 2> approx_cells;    ///< Cell-by-cell deviation matrix (R: APPROX.CELLS)
    double approx_tot;               ///< Global approximation index (R: APPROX.TOT)
    double approx_tot_corr;          ///< Global index corrected with respect to the diagonal (R: APPROX.TOT.CORR)
    std::vector<double> unique_alphas; ///< List of unique alpha cut-off values (R: alpha_tmp)
};

/**
 * @brief Performs the core mathematical phase of the First Order Dominance (FOD) analysis.
 * * @param fod_matrix The original fuzzy dominance matrix (R: FOD.MATRIX).
 * @param tolerance Tolerance for floating-point comparisons to prevent instability (default: 1e-9).
 * @return FODAnalysis A structure containing the key metrics and matrices of the analysis.
 */
FODAnalysis AnalyzeFOD(const Tensor<double, 2>& fod_matrix, double tolerance = 1e-9);
