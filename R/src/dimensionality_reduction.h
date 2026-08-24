/**
 * @file dimensionality_reduction.h
 * @brief High-Performance Computing (HPC) structures and algorithms for dimensionality reduction.
 *
 * @details This module provides advanced tools for analyzing and reducing data dimensionality
 * using Linear Extensions (LE) and loss functions based on Mutual Ranking Probability (MRP).
 * The architecture is optimized for HPC using C++20 standards, including bitwise operations,
 * template-based zero-overhead abstractions, and thread-safe result collection.
 *
 * @par Bibliographic References
 * - Sørensen, P.B., Lerche, D.B., Gyldenkærne, S., Thomsen, M., Fauser, P., Mogensen, B.B.,
 * Kronvang, B., Brüggemann, R., et al. (2004). "Order Theory in Environmental Sciences:
 * Integrative approaches." NERI Technical Report No. 479. National Environmental Research Institute, Denmark.
 * - De Loof, K., De Baets, B., & De Meyer, H. (2006). "Properties of mutual rank
 * probabilities in partially ordered sets".
 */

#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <algorithm>
#include <cstddef>


#include "tensor.h"
#include "display_message.h"
#include "loss_function_mrp.h"

#include <cstdint>
#include <limits>
#include <map>
#include <unordered_map>
#include <mutex>
#include <tuple>
#include <vector>
#include <span>
#include <string_view>
#include <bit>
#include <atomic>
#include <memory>

/**
 * @enum LPOMStrategy
 * @brief Defines the Local Partial Order Model (LPOM) strategy for ranking probability estimation.
 */
enum class LPOMStrategy {
    /** Uses Absolute Ideals (Equation 13' in NERI 479). Faster, assumes independent degrees of freedom. */
    Absolute,
    /** Uses Relative Ideals (Equation 13 in NERI 479). Rigorous, based on set differences. */
    Relative
};

/**
 * @class DimensionalityReductionResult
 * @brief Thread-safe collector for calculated linear extensions, loss values, and best results.
 */
class DimensionalityReductionResult {
public:
    /**
     * @brief Constructs a result collector linked to a global atomic counter.
     */
    explicit DimensionalityReductionResult(std::atomic<std::uint64_t>& all_le_elaborated);
    
    /**
     * @brief Adds a new linear extension and its associated loss value to the collection.
     */
    void AddLe(std::vector<std::uint64_t> le, double loss);
    
    /**
     * @brief Merges the results from another thread-local collector into this one.
     * @param other The thread-local DimensionalityReductionResult to merge.
     */
    void AddLe(DimensionalityReductionResult& other);
    
    /**
     * @brief Builds the final rank profiles mapping for the best permutation found.
     */
    void BuildBestProfileResults(
                                 const std::unordered_map<std::uint64_t, double>& weights,
                                 LossFunctionMRPV2& loss_func,
                                 Tensor<double, 2>& matrice,
                                 const std::vector<std::vector<std::uint64_t>>& list_a,
                                 const std::vector<std::vector<std::uint64_t>>& list_b);
    
    // Getters (implementazione inline in fondo al file per pulizia)
    [[nodiscard]] std::vector<std::vector<std::uint64_t>>& GetLeElaborated();
    [[nodiscard]] std::vector<double>& GetLossValues();
    [[nodiscard]] double GetBestLoss() const;
    [[nodiscard]] const std::vector<std::uint64_t>& GetBestLe() const;
    
    /**
     * @brief Returns the detailed mapping for the best permutation.
     * @return Map: {ProfileID -> tuple<Rank, RankInv, Weight, Error/Prob>}
     */
    [[nodiscard]] const std::map<std::uint64_t, std::tuple<std::uint64_t, std::uint64_t, double, double>>&
    GetBestProfileResults() const;
    
    /**
     * @brief Factory method to create thread-local workers sharing the same progress counter.
     * @return A unique pointer to a new thread-local DimensionalityReductionResult instance.
     */
    [[nodiscard]] std::unique_ptr<DimensionalityReductionResult> SpawnWorkerResult() const;
    
private:
    std::vector<std::vector<std::uint64_t>> le_elaborated_;
    std::vector<double> loss_values_;
    
    std::vector<std::uint64_t> best_permutation_;
    
    ///< Detailed results for the best profile: {ProfileID -> (Rank, RankInv, Weight, Prob)}
    std::map<std::uint64_t, std::tuple<std::uint64_t, std::uint64_t, double, double>> best_profile_results_;
    
    double best_loss_ = std::numeric_limits<double>::max();
    std::atomic<std::uint64_t>& all_le_elaborated_;
};

// =========================================================================
// HPC Probabilistic Engines (Bitwise O(1) in Boolean Lattice)
// =========================================================================

/**
 * @brief Computes the number of elements strictly above @p n in the Boolean
 *        lattice induced by the permutation bitmask representation.
 *
 * For an element @p n encoded as a bitmask of @p cardinality bits, this
 * counts the "up-set" cardinality using the positions of the zero bits.
 * Let z_0 < z_1 < ... < z_{k-1} be the positions of the zero bits of @p n;
 * the formula is:
 *
 *   npup(n) = k + Σ_{i<j} 2^(z_j - z_i - 1)
 *
 * @note This is NOT simply 2^(cardinality - popcount(n)) - 1, because the
 *       combinatorial term Σ 2^(z_j - z_i - 1) encodes the structure of the
 *       intersection of the two linear extensions (le and le_inv), not the
 *       full Boolean lattice. The two coincide only when all zero bits are
 *       contiguous.
 *
 * @param n           Bitmask representation of the element (after permutation).
 * @param cardinality Number of relevant bits (size of the permutation).
 * @return            Number of elements strictly above @p n in the induced poset.
 */
[[nodiscard]] inline std::uint64_t npup_fast(std::uint64_t n, std::uint64_t cardinality) noexcept {
    const uint64_t mask = (1ULL << cardinality) - 1ULL;
    uint64_t zeros = ~n & mask;          // posizioni degli zeri come bitmask
    uint64_t result = 0;
    uint64_t lower_zeros = 0;            // bitmask degli zeri già processati
    
    while (zeros) {
        const uint64_t lsb = zeros & -zeros;          // zero corrente z_j = 2^pos
        const uint64_t pos_j = std::countr_zero(lsb);
        
        result += 1; // contributo diretto
        
        // Σ_{z_i in lower_zeros} 2^(pos_j - pos_i - 1)
        uint64_t prev = lower_zeros;
        while (prev) {
            const uint64_t lsb_i = prev & -prev;
            const uint64_t pos_i = std::countr_zero(lsb_i);
            result += 1ULL << (pos_j - pos_i - 1);
            prev &= prev - 1;
        }
        
        lower_zeros |= lsb;
        zeros &= zeros - 1;
    }
    return result;
}

/**
 * @brief Estimates P(row > col) using Absolute Ideals — Eq. 13' of NERI 479.
 *
 * Eq. 13' is defined as:
 * @code
 *   probQ+(x>y) = 1 / (1 + Q(x) / Q(y))
 * @endcode
 * where Q(x) = (Nu(x) + 1) / (Nd(x) + 1), with Nu(x) and Nd(x) counting
 * ALL elements above and below x in the poset (including those also above/below y).
 *
 * In the paper, Nu(x) and Nd(x) are counted over the actual elements of the
 * poset. Here they are computed via npup_fast on the permutation bitmask:
 *   - Nu(x) ≈ npup_fast(row, num_vars)
 *   - Nd(x) ≈ npup_fast(~row & mask, num_vars)
 *
 * @note The deviation from the exact Eq. 13' arises for the same reason as
 *       msi_relative: npup_fast operates on the combinatorial structure of
 *       the bitmask (zero-bit positions and their pairwise gaps), which is a
 *       proxy for the true up/down-set sizes of the poset induced by the
 *       intersection of le and le_inv. Specifically, the term
 *       Σ_{i<j} 2^(z_j - z_i - 1) has no direct counterpart in the paper's
 *       definition of Nu — it emerges from the encoding of the two linear
 *       extensions into a single bitmask, and causes npup_fast to diverge
 *       from the simple 2^(num_vars - popcount(n)) - 1 that the full Boolean
 *       lattice would give.
 *
 * @param row      Bitmask of the first element (after permutation).
 * @param col      Bitmask of the second element (after permutation).
 * @param num_vars Number of relevant bits (size of the permutation).
 * @return         Estimated probability that row ranks above col.
 */
[[nodiscard]] inline double msi_absolute(std::uint64_t row, std::uint64_t col, std::uint64_t num_vars) noexcept {
    const double q_row = static_cast<double>(1 + npup_fast(row, num_vars))
    / static_cast<double>(1 + npup_fast(~row, num_vars));
    const double q_col = static_cast<double>(1 + npup_fast(col, num_vars))
    / static_cast<double>(1 + npup_fast(~col, num_vars));
    return q_row / (q_row + q_col);
}

/**
 * @brief Estimates P(row > col) using Relative Ideals — Eq. 13 of NERI 479.
 *
 * Eq. 13 is defined as:
 * @code
 *   probQ(x>y) = 1 / (1 + Q(x o y) / Q(y o x))
 * @endcode
 * where Q(x o y) = (Nu(x o y) + 1) / (Nd(x o y) + 1), and:
 *   - Nu(x o y) = number of elements above x but NOT above y
 *   - Nd(x o y) = number of elements below x but NOT below y
 *
 * In the paper, Nu and Nd are counted over the actual elements of the poset.
 * Here they are approximated on the Boolean lattice induced by the permutation
 * bitmask, using npup_fast:
 *   - Nu(row o col) ≈ npup_fast(row) - npup_fast(row | col)
 *   - Nd(row o col) ≈ npup_fast(~row) - npup_fast(~(row & col))
 *
 * @note The deviation from the exact Eq. 13 arises because npup_fast counts
 *       elements in the combinatorial structure of the bitmask (positions of
 *       zero bits and their pairwise gaps), which approximates — but does not
 *       exactly reproduce — the true Nu/Nd of the poset defined by the
 *       intersection of le and le_inv. In particular, elements that are
 *       incomparable in the real poset but ordered in the bitmask lattice are
 *       incorrectly counted.
 *
 * @param row      Bitmask of the first element (after permutation).
 * @param col      Bitmask of the second element (after permutation).
 * @param num_vars Number of relevant bits (size of the permutation).
 * @return         Estimated probability that row ranks above col.
 */
[[nodiscard]] inline double msi_relative(std::uint64_t row, std::uint64_t col, std::uint64_t num_vars) noexcept {
    const uint64_t mask = (1ULL << num_vars) - 1ULL;
    
    // Nu(xoy) = npup(row) - npup(row|col)   [sopra row ma non sopra col]
    // Nd(xoy) = npdown(row) - npdown(row&col) [sotto row ma non sotto col]
    const double nu_row     = static_cast<double>(npup_fast(row,        num_vars));
    const double nu_col     = static_cast<double>(npup_fast(col,        num_vars));
    const double nu_row_or  = static_cast<double>(npup_fast(row | col,  num_vars));
    const double nd_row     = static_cast<double>(npup_fast(~row & mask, num_vars));
    const double nd_col     = static_cast<double>(npup_fast(~col & mask, num_vars));
    const double nd_row_and = static_cast<double>(npup_fast(~(row & col) & mask, num_vars));
    
    const double nu_rel_row = nu_row - nu_row_or;   // Nu(row o col)
    const double nu_rel_col = nu_col - nu_row_or;   // Nu(col o row)
    const double nd_rel_row = nd_row - nd_row_and;  // Nd(row o col)
    const double nd_rel_col = nd_col - nd_row_and;  // Nd(col o row)
    
    // Q(xoy) = (Nu(xoy)+1) / (Nd(xoy)+1)
    const double q_row = (nu_rel_row + 1.0) / (nd_rel_row + 1.0);
    const double q_col = (nu_rel_col + 1.0) / (nd_rel_col + 1.0);
    
    // probQ(row>col) = 1 / (1 + Q(row o col) / Q(col o row))
    //                = Q(col o row) / (Q(row o col) + Q(col o row))
    return q_col / (q_row + q_col);
}

// =========================================================================
// Core Algorithms
// =========================================================================

/**
 * @brief Constructs a Linear Extension (LE) based on a bit-permutation.
 */
void DimensionalityReductionBuildLE(
                                    std::uint64_t numero_variabili,
                                    std::span<const std::uint64_t> permutazione,
                                    std::span<const std::uint64_t> permutazione_inv,
                                    std::vector<std::uint64_t>& le,
                                    std::vector<std::uint64_t>& le_inv);

/**
 * @brief Builds the MRP matrix using deterministic LE intersection and LPOM fallback.
 */
template <LPOMStrategy Strategy>
void DimensionalityReductionBuildMRPIntersection(
                                                 std::span<const std::uint64_t> permutazione,
                                                 std::span<const std::uint64_t> le,
                                                 std::span<const std::uint64_t> le_inv,
                                                 Tensor<double, 2>& mrp_le_intersection,
                                                 const std::vector<std::uint_fast64_t>& elements_used) {
    
    //const std::uint64_t n_elements = mrp_le_intersection.Rows();
    const std::size_t perm_size = permutazione.size();
    
    std::vector<std::uint64_t> shift_positions(perm_size);
    for (std::size_t i = 0; i < perm_size; ++i) {
        shift_positions[i] = perm_size - permutazione[i] - 1;
    }
    
    auto permuta_bit = [&](std::uint64_t n) {
        std::uint64_t ris = 0;
        for (const std::uint64_t pos : shift_positions) {
            ris = (ris << 1) | ((n >> pos) & 1ULL);
        }
        return ris;
    };
    
    // ESECUZIONE SOLO SUI PROFILI ATTIVI
    for (std::size_t re = 0; re < elements_used.size(); ++re) {
        const std::uint64_t row = elements_used[re];
        
        mrp_le_intersection(row, row) = 1.0;
        
        const std::uint64_t r_le = le[row];
        const std::uint64_t r_le_inv = le_inv[row];
        const std::uint64_t p_row = permuta_bit(row);
        
        for (std::size_t ce = re + 1; ce < elements_used.size(); ++ce) {
            const std::uint64_t col = elements_used[ce];
        
            const std::uint64_t c_le = le[col];
            const std::uint64_t c_le_inv = le_inv[col];
            
            if (r_le < c_le && r_le_inv < c_le_inv) {
                mrp_le_intersection(row, col) = 1.0;
                mrp_le_intersection(col, row) = 0.0;
            } else if (c_le < r_le && c_le_inv < r_le_inv) {
                mrp_le_intersection(row, col) = 0.0;
                mrp_le_intersection(col, row) = 1.0;
            } else {
                const std::uint64_t p_col = permuta_bit(col);
                double v = 0.0;
                if constexpr (Strategy == LPOMStrategy::Absolute) {
                    v = msi_absolute(p_row, p_col, perm_size);
                } else {
                    v = msi_relative(p_row, p_col, perm_size);
                }
                mrp_le_intersection(row, col) = v;
                mrp_le_intersection(col, row) = 1.0 - v;
            }
        }
    }
}

/**
 * @brief Entry point for the Exact Dimensionality Reduction algorithm.
 */
void ExactDimensionalityReduction(
                                  const std::unordered_map<std::uint64_t, double>& weights,
                                  std::uint64_t numero_variabili,
                                  std::string_view loss_str,
                                  int lpom_type,
                                  DisplayMessage* display_message,
                                  double thread_percentage,
                                  DimensionalityReductionResult& result);

void BidimentionalPosetRepresentation(const std::unordered_map<std::uint64_t, double>& weights,
                                      std::uint64_t numero_variabili,
                                      std::string_view loss_str,
                                      int lpom_type,
                                      std::vector<std::uint64_t>& variable_priority,
                                      DimensionalityReductionResult& result);


// =========================================================================
// Inline Definitions
// =========================================================================

inline std::vector<std::vector<std::uint64_t>>& DimensionalityReductionResult::GetLeElaborated() { return le_elaborated_; }
inline std::vector<double>& DimensionalityReductionResult::GetLossValues() { return loss_values_; }
inline double DimensionalityReductionResult::GetBestLoss() const { return best_loss_; }
inline const std::vector<std::uint64_t>& DimensionalityReductionResult::GetBestLe() const { return best_permutation_; }
inline const std::map<std::uint64_t, std::tuple<std::uint64_t, std::uint64_t, double, double>>&
DimensionalityReductionResult::GetBestProfileResults() const { return best_profile_results_; }
inline std::unique_ptr<DimensionalityReductionResult> DimensionalityReductionResult::SpawnWorkerResult() const {
    return std::make_unique<DimensionalityReductionResult>(all_le_elaborated_);
}

