/**
 * @file rwrapper.h
 * @brief R/C++ interface wrapper for the POSetic Data Analysis library.
 *
 * Provides the full glue layer between R and C++:
 * - @ref RProtectGuard     — RAII wrapper for PROTECT/UNPROTECT.
 * - @ref forward_exception_to_r / @ref forward_warning_to_r — exception forwarding.
 * - @ref pointerFinalizer  — generic R external-pointer finalizer.
 * - @ref makeDisplayMessage / @ref buildLEMatrix — shared utilities.
 * - @ref RConvert           — R → C++ type conversion helpers.
 * - @ref RCreate            — C++ → R object creation helpers.
 *
 * ### RAII pattern example
 * @code
 * {
 *     RProtectGuard guard;
 *     SEXP vec = guard.protect(Rf_allocVector(REALSXP, 10));
 *     SEXP mat = guard.protect(Rf_allocMatrix(INTSXP, 5, 5));
 *     // ... use vec and mat ...
 * } // UNPROTECT(2) called automatically
 * @endcode
 *
 * @author Alessandro Avellone
 * @version 2.1
 * @date 2025
 */

#pragma once


#include <cstdint>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>

#include "tensor.h"
#include "display_message.h"
#include "r_display.h"

// ===========================================================================
// RAII guard for R memory protection
// ===========================================================================

/**
 * @class RProtectGuard
 * @brief RAII wrapper for automatic R PROTECT / UNPROTECT management.
 *
 * Tracks how many R objects have been protected and calls `UNPROTECT` with
 * the correct count when the guard goes out of scope — even on exception paths.
 *
 * @note Non-copyable and non-movable to prevent double-unprotect bugs.
 */
class RProtectGuard {
private:
    int count_ = 0;  ///< Number of protected objects
    
public:
    RProtectGuard() = default;
    
    /// Calls `UNPROTECT(count_)` if any objects were protected.
    ~RProtectGuard() {
        if (count_ > 0) {
            UNPROTECT(count_);
        }
    }
    
    // Non-copyable, non-movable
    RProtectGuard(const RProtectGuard&) = delete;
    RProtectGuard& operator=(const RProtectGuard&) = delete;
    RProtectGuard(RProtectGuard&&) = delete;
    RProtectGuard& operator=(RProtectGuard&&) = delete;
    
    /**
     * @brief Protects @p sexp and tracks it for automatic unprotection.
     * @param sexp  R object to protect.
     * @return The same @p sexp (enables one-liner usage).
     */
    SEXP Protect(SEXP sexp) {
        ++count_;
        return PROTECT(sexp);
    }
    
    /** @brief Returns the number of currently protected objects. */
    int count() const noexcept { return count_; }
};

// ===========================================================================
// Exception / warning forwarding
// ===========================================================================

/**
 * @brief Forwards a C++ exception message to R as a fatal error.
 *
 * Calls `Rf_error()`, which **does not return** (it longjmps).
 * @param message  Human-readable error description.
 */
inline void forward_exception_to_r(const std::string& message) {
    Rf_error("%s", message.c_str());
}

/**
 * @brief Forwards a warning message to R via `Rf_warning()`.
 * @param message  Human-readable warning text.
 */
inline void forward_warning_to_r(const std::string& message) {
    Rf_warning("%s", message.c_str());
}

// ===========================================================================
// R-callable function declarations
// ===========================================================================

extern "C" {
    
    // Relation property checks
    SEXP isReflexive(SEXP elements_r, SEXP comparabilities_r);
    SEXP isSymmetric(SEXP comparabilities_r);
    SEXP isAntisymmetric(SEXP comparabilities_r);
    SEXP isTransitive(SEXP comparabilities_r);
    SEXP TransitiveClosure(SEXP comparabilities_r);
    SEXP ReflexiveClosure(SEXP elements_r, SEXP comparabilities_r);
    SEXP isPreorder(SEXP elements_r, SEXP comparabilities_r);
    SEXP isPartialOrder(SEXP elements_r, SEXP comparabilities_r);
    
    // POSet construction
    SEXP BuildPOSet(SEXP elements, SEXP comparabilities);
    SEXP BuildLinearPOSet(SEXP elements);
    SEXP BuildProductPOSet(SEXP posets_r);
    SEXP BuildLexicographicProductPOSet(SEXP posets_r);
    SEXP BuildIntersectionPOSet(SEXP posets_r);
    SEXP BuildLinearSumPOSet(SEXP posets_r);
    SEXP BuildDisjointSumPOSet(SEXP posets_r);
    SEXP BuildLiftingPOSet(SEXP poset_r, SEXP new_element_r);
    SEXP BuildBucketPOSet(SEXP elements, SEXP comparabilities);
    SEXP BuildBinaryVariablePOSet(SEXP variables);
    SEXP BuildFencePOSet(SEXP elements_r, SEXP orientation_r);
    SEXP BuildCrownPOSet(SEXP elements_1_r, SEXP elements_2_r);
    SEXP BuildDualPOSet(SEXP poset_r);
    
    // POSet properties
    SEXP Elements(SEXP poset_r);
    SEXP IncidenceMatrix(SEXP poset_r);
    SEXP OrderRelation(SEXP poset_r);
    SEXP IsDominatedBy(SEXP poset_r, SEXP v1, SEXP v2);
    SEXP Dominates(SEXP poset_r, SEXP v1, SEXP v2);
    SEXP IsComparableWith(SEXP poset_r, SEXP v1, SEXP v2);
    SEXP IsIncomparableWith(SEXP poset_r, SEXP v1, SEXP v2);
    SEXP UpsetOf(SEXP poset_r, SEXP insieme);
    SEXP IsUpset(SEXP poset_r, SEXP insieme);
    SEXP IsDownset(SEXP poset_r, SEXP insieme);
    SEXP DownsetOf(SEXP poset_r, SEXP insieme);
    SEXP ComparabilitySetOf(SEXP poset_r, SEXP elemento);
    SEXP IncomparabilitySetOf(SEXP poset_r, SEXP elemento);
    SEXP Maximal(SEXP poset_r);
    SEXP Minimal(SEXP poset_r);
    SEXP IsMaximal(SEXP poset_r, SEXP elemento);
    SEXP IsMinimal(SEXP poset_r, SEXP elemento);
    SEXP Meet(SEXP poset_r, SEXP insieme_r);
    SEXP Join(SEXP poset_r, SEXP insieme_r);
    SEXP CoverRelation(SEXP poset_r);
    SEXP CoverMatrix(SEXP poset_r);
    SEXP Incomparabilities(SEXP poset_r);
    SEXP IsExtensionOf(SEXP poset_r_1, SEXP poset_r_2);
    
    // Linear extension generators
    SEXP BuildLEGenerator(SEXP poset_r);
    SEXP BuildBubleyDyerLEGenerator(SEXP poset_r, SEXP arg);
    SEXP LEGBubleyDyerGet(SEXP generator_r, SEXP from_start_r, SEXP quante_r, SEXP errore_r, SEXP output_ogni_in_sec_r);
    SEXP LEGGet(SEXP generator_r, SEXP from_start_r, SEXP quante_r, SEXP output_ogni_in_sec_r);
    
    // Evaluation functions
    SEXP ExactMRP(SEXP poset_r, SEXP output_ogni_in_sec_r);
    SEXP BuildBubleyDyerMRPGenerator(SEXP poset_r, SEXP seed_r);
    SEXP BubleyDyerMRP(SEXP bubley_dyer_r, SEXP quante_r, SEXP errore_r, SEXP output_ogni_in_sec_r);
    SEXP ExactEvaluation(SEXP poset_r, SEXP output_ogni_in_sec_r, SEXP functions_r);
    SEXP BuildBubleyDyerEvaluationGenerator(SEXP poset_r, SEXP seed_r, SEXP functions_r);
    SEXP BubleyDyerEvaluation(SEXP generator_r, SEXP quante_r, SEXP errore_r, SEXP output_ogni_in_sec_r);
    SEXP BruggemannLercheSorensenDominance(SEXP poset_r, SEXP type_r);
    
    // Separation functions
    SEXP FuzzySeparation(SEXP dominance_matrix_r, SEXP quale_r, SEXP tnorm_r, SEXP tconorm_r);
    SEXP FuzzyInBetweenness(SEXP dominance_matrix_r, SEXP tnorm_r, SEXP tconorm_r, SEXP quali_r);
    SEXP ExactSeparation(SEXP poset_r, SEXP output_ogni_in_sec_r, SEXP quali_r);
    SEXP RunLexSeparation(SEXP modalita_r);
    SEXP RunLexMRP(SEXP modalita_r);
    SEXP BuildBubleyDyerSeparationGenerator(SEXP poset_r, SEXP seed_r, SEXP quali_r);
    SEXP BubleyDyerSeparation(SEXP generator_r, SEXP quante_r, SEXP errore_r, SEXP output_ogni_in_sec_r);
    
    // Dimensionality reduction
    SEXP RunDimensionalityReduction(SEXP profile_r, SEXP weights_r, SEXP loss_r, SEXP lpom_strategy_r, SEXP output_ogni_in_sec_r, SEXP thread_percentage_r);
    SEXP RunBidimentionalPosetRepresentation(SEXP profile_r, SEXP weights_r, SEXP loss_r, SEXP lpom_strategy_r, SEXP variable_priority_r);

    // First Order Dominance Analisys
    SEXP FirstOrderDominanceAnalysis(SEXP poset_r, SEXP freq_matrix_r, SEXP subpopulation_count_r, SEXP metric_r, SEXP total_bins_r, SEXP count_r, SEXP seed_r, SEXP output_interval, SEXP sep_r, SEXP linear_extensions_r, SEXP n_threads_r);
    
}
