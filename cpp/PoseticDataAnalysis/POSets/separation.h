/**
 * @file separation.h
 * @brief HPC-optimized algorithms for fuzzy in-betweenness and lexical separation.
 *
 * This module provides functions to calculate generalized fuzzy in-betweenness,
 * cumulative separations, and lexical separations using C++20 features such as
 * Concepts and zero-copy move semantics.
 * * @see Fodor, J., & Roubens, M. (1994). "Fuzzy Preference Modelling and Multicriteria Decision Support".
 * @see Klement, E. P., Mesiar, R., & Pap, E. (2000). "Triangular Norms".
 */

#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <array>

#include <algorithm>
#include <concepts>
#include <cstdint>
#include <tuple>
#include <vector>

#include "tensor.h"

// ---------------------------------------------------------------------------
// Concepts and Functors for Norms/Conorms
// ---------------------------------------------------------------------------

/**
 * @concept NormConormFunc
 * @brief Concept defining a valid t-norm or t-conorm callable.
 *
 * Ensures the provided callable takes two doubles and returns a double.
 * This replaces virtual inheritance for maximum inlining.
 */
template <typename F>
concept NormConormFunc = std::regular_invocable<F, double, double> &&
std::convertible_to<std::invoke_result_t<F, double, double>, double>;

/**
 * @brief Functor for the Minimum t-norm / t-conorm (Gödel t-norm).
 */
struct MinNormConorm {
    constexpr double operator()(double a, double b) const noexcept { return std::min(a, b); }
};

/**
 * @brief Functor for the Maximum t-norm / t-conorm.
 */
struct MaxNormConorm {
    constexpr double operator()(double a, double b) const noexcept { return std::max(a, b); }
};

/**
 * @brief Functor for the Product t-norm.
 */
struct ProdNormConorm {
    constexpr double operator()(double a, double b) const noexcept { return a * b; }
};

/**
 * @brief Functor for the Probabilistic Sum t-conorm.
 */
struct ProbNormConorm {
    constexpr double operator()(double a, double b) const noexcept { return a + b - a * b; }
};

// ---------------------------------------------------------------------------
// Structs for Return Types (Zero-copy Move Semantics)
// ---------------------------------------------------------------------------

/**
 * @brief Results of a general separation computation.
 */
struct GeneralSeparationResult {
    Tensor<double, 2> sep_all;   /**< Total separation matrix */
    Tensor<double, 2> sep_lower; /**< Lower separation matrix */
    Tensor<double, 2> sep_upper; /**< Upper separation matrix */
};

/**
 * @brief Results of a lexical separation computation.
 */
struct LexSeparationResult {
    Tensor<double, 2> sep_all;        /**< Total separation matrix */
    Tensor<double, 2> sep_lower;      /**< Lower separation matrix */
    Tensor<double, 2> sep_upper;      /**< Upper separation matrix */
    Tensor<double, 2> sep_vertical;   /**< Vertical separation matrix */
    Tensor<double, 2> sep_horizontal; /**< Horizontal separation matrix */
    std::vector<std::vector<std::uint64_t>> profili; /**< Generated profiles mapping */
};

/**
 * @brief Results of a Lexical Mixed-Radix Profile (MRP) computation.
 */
struct LexMrpResult {
    Tensor<double, 2> mrp;                             /**< MRP matrix */
    std::vector<std::vector<std::uint64_t>> profili; /**< Generated profiles mapping */
};

// ---------------------------------------------------------------------------
// Inline Template Implementations (For Maximum HPC Inlining)
// ---------------------------------------------------------------------------

/**
 * @brief Calculates the general fuzzy in-betweenness for a triplet of nodes.
 *
 * Measures the strength of the indirect path via an intermediate node r in a
 * fuzzy directed graph (dominance matrix).
 * * @see Rosenfeld, A. (1975). "Fuzzy Graphs", for foundational fuzzy path formulations.
 *
 * @tparam TTimes The t-norm functor type.
 * @tparam TPlus The t-conorm functor type.
 * @param pi Index of the first node (p).
 * @param qi Index of the second node (q).
 * @param ri Index of the intermediate node (r).
 * @param dominance The dominance matrix representing fuzzy preference/dominance degrees.
 * @param times The t-norm instance for intersection (e.g., path bottleneck).
 * @param plus The t-conorm instance for union.
 * @param finb_prq [out] Fuzzy in-betweenness p -> r -> q.
 * @param finb_qrp [out] Fuzzy in-betweenness q -> r -> p.
 * @param finbqrp [out] Combined fuzzy in-betweenness.
 */
template <NormConormFunc TTimes, NormConormFunc TPlus>
inline void GeneralFuzzyInBetweenness(std::uint64_t pi, std::uint64_t qi, std::uint64_t ri,
                                      const Tensor<double, 2>& dominance,
                                      TTimes times, TPlus plus,
                                      double& finb_prq, double& finb_qrp, double& finbqrp) noexcept {
    double sdom_pr = (pi != ri) ? dominance(pi, ri) : 0.0;
    double sdom_rq = (qi != ri) ? dominance(ri, qi) : 0.0;
    finb_prq = times(times(dominance(pi, qi), sdom_pr), sdom_rq);
    
    double sdom_qr = (qi != ri) ? dominance(qi, ri) : 0.0;
    double sdom_rp = (pi != ri) ? dominance(ri, pi) : 0.0;
    finb_qrp = times(times(dominance(qi, pi), sdom_qr), sdom_rp);
    
    finbqrp = plus(finb_prq, finb_qrp);
}

/**
 * @brief Calculates the cumulative fuzzy in-betweenness for a pair over all possible intermediates.
 */
template <NormConormFunc TTimes, NormConormFunc TPlus>
inline void CumulativeFuzzyInBetweenness(std::uint64_t pi, std::uint64_t qi,
                                         const Tensor<double, 2>& dominance,
                                         TTimes times, TPlus plus,
                                         double& cfinbpq, double& cfinb_pq, double& cfinb_qp) noexcept {
    cfinbpq = 0.0;
    cfinb_pq = 0.0;
    cfinb_qp = 0.0;
    std::uint64_t n_rows = dominance.Extent(0);
    
    for (std::uint64_t ri = 0; ri < n_rows; ++ri) {
        double finb_prq = 0.0, finb_qrp = 0.0, finbqrp = 0.0;
        GeneralFuzzyInBetweenness(pi, qi, ri, dominance, times, plus, finb_prq, finb_qrp, finbqrp);
        cfinb_pq += finb_prq;
        cfinb_qp += finb_qrp;
        cfinbpq += finbqrp;
    }
}

/**
 * @brief Computes general separations for all node pairs based on dominance and cumulative in-betweenness.
 * * @note Uses `Tensor<double, 2>(std::array<std::uint64_t, 2>{0, 0}, 0.0)` for unrequested matrices to prevent constructor ambiguity.
 */
template <NormConormFunc TTimes, NormConormFunc TPlus>
GeneralSeparationResult GeneralSeparation(const Tensor<double, 2>& dominance,
                                          TTimes times, TPlus plus,
                                          bool do_all, bool do_lower, bool do_upper) {
    std::uint64_t n = dominance.Extent(0);
    GeneralSeparationResult result{
        .sep_all   = do_all   ? Tensor<double, 2>(std::array<std::uint64_t, 2>{n, n}, kUninitialized) : Tensor<double, 2>(std::array<std::uint64_t, 2>{0, 0}, 0.0),
        .sep_lower = do_lower ? Tensor<double, 2>(std::array<std::uint64_t, 2>{n, n}, kUninitialized) : Tensor<double, 2>(std::array<std::uint64_t, 2>{0, 0}, 0.0),
        .sep_upper = do_upper ? Tensor<double, 2>(std::array<std::uint64_t, 2>{n, n}, kUninitialized) : Tensor<double, 2>(std::array<std::uint64_t, 2>{0, 0}, 0.0)
    };
    
    for (std::uint64_t pi = 0; pi < n; ++pi) {
        for (std::uint64_t qi = pi + 1; qi < n; ++qi) {
            double cfinbpq = 0.0, cfinb_pq = 0.0, cfinb_qp = 0.0;
            CumulativeFuzzyInBetweenness(pi, qi, dominance, times, plus, cfinbpq, cfinb_pq, cfinb_qp);
            
            if (do_all) {
                double val = dominance(pi, qi) + dominance(qi, pi) + cfinbpq;
                result.sep_all(pi, qi) = val;
                result.sep_all(qi, pi) = val;
            }
            if (do_lower) {
                result.sep_lower(pi, qi) = dominance(pi, qi) + cfinb_pq;
                result.sep_lower(qi, pi) = dominance(qi, pi) + cfinb_qp;
            }
            if (do_upper) {
                result.sep_upper(pi, qi) = dominance(qi, pi) + cfinb_qp;
                result.sep_upper(qi, pi) = dominance(pi, qi) + cfinb_pq;
            }
        }
    }
    return result;
}

// ---------------------------------------------------------------------------
// Lexical Separations Declarations
// ---------------------------------------------------------------------------

/**
 * @brief Computes lexical separation when all variables share the same number of modalities.
 *
 * @see Ehrgott, M. (2005). "Multicriteria Optimization" for lexicographic orders.
 * @param numero_variabili The total number of variables.
 * @param numero_modalita The constant number of modalities for each variable.
 * @return LexSeparationResult containing matrices and generated profiles.
 */
LexSeparationResult LexSeparationEqDeg(std::uint64_t numero_variabili, std::uint64_t numero_modalita);

/**
 * @brief Computes lexical separation with distinct degrees for each variable.
 *
 * @param numero_modalita A vector defining the number of modalities per variable.
 * @return LexSeparationResult containing matrices and generated profiles.
 */
LexSeparationResult LexSeparationDeg(const std::vector<std::uint64_t>& numero_modalita);

/**
 * @brief Computes Lexical Mixed-Radix Profiles (MRP).
 *
 * @see Knuth, D. E. (1997). "The Art of Computer Programming, Volume 2: Seminumerical Algorithms"
 * (Section 4.1 for Mixed-Radix Numeral Systems).
 * @param numero_modalita A vector defining the number of modalities per variable.
 * @return LexMrpResult containing the MRP matrix and generated profiles.
 */
LexMrpResult LexMrp(const std::vector<std::uint64_t>& numero_modalita);
