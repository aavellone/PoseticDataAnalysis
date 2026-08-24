/**
 * @file rwrapper_poset_properties.cpp
 * @brief R/C++ interface — POSet property query functions (HPC Optimized).
 *
 * @details Implements all R-callable (`extern "C"`) entry points for querying POSet properties:
 * - Basic structure: elements, incidence matrix, order relation, cover relation.
 * - Pairwise comparisons: dominates, is-dominated-by, is-comparable, is-incomparable.
 * - Set properties: upset, downset, comparability set, incomparability set.
 * - Extremal elements: maximal, minimal, meet, join.
 * - Structural queries: incomparabilities, extension check.
 *
 * Every public function follows the strict zero-overhead C++20 pattern:
 * 1. Unwrap R arguments safely via `RConvert` helpers.
 * 2. Query the underlying C++ POSet without redundant string copies.
 * 3. Convert results natively into R objects mapped continuously in memory.
 *
 * @author Alessandro Avellone
 * @version 3.0 (HPC C++20)
 * @date 2025
 */

#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <list>
#include <utility>
#include <stdexcept>

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>

#include "rwrapper.h"
#include "my_exception.h"
#include "poset_wrapper.h"
#include "poset.h"
#include "rwrapper_conversion.h"

// ===========================================================================
// R-callable functions (C-linkage)
// ===========================================================================

extern "C" {
    /**
     * @brief Retrieves the elements of the POSet.
     *
     * @param poset_r ExternalPtr to the POSetWrap object.
     * @return SEXP Character vector containing the names of the elements.
     */
    SEXP Elements(SEXP poset_r) {
        SEXP result_r = R_NilValue;
        RProtectGuard guard;

        try {
            const POSetWrap *poset_wrap = RConvert::ToPOSetWrap(poset_r);
            POSet *poset = poset_wrap->GetPOSet();

            std::size_t n = static_cast<std::size_t>(poset->size());
            result_r = guard.Protect(Rf_allocVector(STRSXP, n));
            for (std::size_t i = 0; i < n; ++i) {
                std::string_view name_view = poset->GetElementName(i);
                SET_STRING_ELT(
                               result_r,
                               static_cast<R_xlen_t>(i),
                               Rf_mkCharLenCE(name_view.data(), static_cast<int>(name_view.size()), CE_UTF8)
                               );
            }
            
        } catch (std::exception &ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }

        return result_r;
    }
    
    /**
     * @brief Retrieves the incidence matrix of the POSet.
     *
     * @param poset_r ExternalPtr to the POSetWrap object.
     * @return SEXP An integer matrix representing the directed graph.
     */
    SEXP IncidenceMatrix(SEXP poset_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;
        
        try {
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);
            POSet* poset = poset_wrap->GetPOSet();

            auto r = poset->IncidenceMatrix();

            const int nrow = static_cast<int>(r.Extent(0));
            const int ncol = static_cast<int>(r.Extent(1));
            
            // Allocate native R matrix
            result_r = guard.Protect(Rf_allocMatrix(INTSXP, nrow, ncol));
            for (int i = 0; i < nrow; ++i) {
                for (int j = 0; j < ncol; ++j) {
                    INTEGER(result_r)[i + nrow * j] = r(static_cast<std::uint64_t>(i), static_cast<std::uint64_t>(j));
                }
            }
            
            // Assign dimnames
            SEXP dimnames = guard.Protect(Rf_allocVector(VECSXP, 2));
            SEXP r_names = guard.Protect(Rf_allocVector(STRSXP, nrow));
            
            for (int i = 0; i < nrow; ++i) {
                std::string_view name_view = poset->GetElementName(static_cast<std::uint64_t>(i));
                
                SET_STRING_ELT(
                               r_names,
                               i,
                               Rf_mkCharLenCE(name_view.data(), static_cast<int>(name_view.size()), CE_UTF8)
                               );
            }
            
            SET_VECTOR_ELT(dimnames, 0, r_names);
            SET_VECTOR_ELT(dimnames, 1, r_names);
            Rf_setAttrib(result_r, R_DimNamesSymbol, dimnames);
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        return result_r;
    }
    
    /**
     * @brief Retrieves the strict order relations (edges) of the POSet.
     *
     * @param poset_r ExternalPtr to the POSetWrap object.
     * @return SEXP N x 2 character matrix defining the (u < v) relations.
     */
    SEXP OrderRelation(SEXP poset_r) {
        SEXP result_r = R_NilValue;
        RProtectGuard guard;
        
        try {
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);
            POSet* poset = poset_wrap->GetPOSet();

            auto relations = poset->OrderRelation();
            
            const int n = static_cast<int>(relations.size());
            result_r = guard.Protect(Rf_allocMatrix(STRSXP, n, 2));
            
            int pos = 0;
            for (const auto& pair : relations) {
                SET_STRING_ELT(result_r, pos, Rf_mkChar(pair.first.c_str()));
                SET_STRING_ELT(result_r, pos + n, Rf_mkChar(pair.second.c_str()));
                pos++;
            }
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    /**
     * @brief Retrieves the cover relations (Hasse diagram edges) of the POSet.
     *
     * @param poset_r ExternalPtr to the POSetWrap object.
     * @return SEXP N x 2 character matrix defining the direct connections.
     */
    SEXP CoverRelation(SEXP poset_r) {
        SEXP result_r = R_NilValue;
        RProtectGuard guard;
        
        try {
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);
            POSet* poset = poset_wrap->GetPOSet();

            auto relations = poset->CoverRelation();
            
            const int n = static_cast<int>(relations.size());
            result_r = guard.Protect(Rf_allocMatrix(STRSXP, n, 2));
            
            int pos = 0;
            for (const auto& p : relations) {
                std::string_view first_view = poset->GetElementName(p.first);
                std::string_view second_view = poset->GetElementName(p.second);
                
                SET_STRING_ELT(result_r, pos,
                               Rf_mkCharLenCE(first_view.data(), static_cast<int>(first_view.size()), CE_UTF8));
                
                SET_STRING_ELT(result_r, pos + n,
                               Rf_mkCharLenCE(second_view.data(), static_cast<int>(second_view.size()), CE_UTF8));
                
                pos++;
            }
            
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    /**
     * @brief Returns the cover matrix as a logical matrix.
     *
     * @details Entry `[i, j] = TRUE` if and only if element `i` is covered by element `j`.
     * This evaluates the Hasse diagram structure directly in a zero-copy manner.
     *
     * @param poset_r ExternalPtr to the POSetWrap object.
     * @return SEXP An R logical matrix with row and column names mapped to the POSet elements.
     */
    SEXP CoverMatrix(SEXP poset_r) {
        SEXP result_r = R_NilValue;
        RProtectGuard guard;
        
        try {
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);
            POSet* poset = poset_wrap->GetPOSet();
            auto r = poset->CoverMatrix();
            
            const int nrow = static_cast<int>(r.Extent(0));
            const int ncol = static_cast<int>(r.Extent(1));
            
            // Allocate native R matrix
            result_r = guard.Protect(Rf_allocMatrix(INTSXP, nrow, ncol));
            for (int i = 0; i < nrow; ++i) {
                for (int j = 0; j < ncol; ++j) {
                    INTEGER(result_r)[i + nrow * j] = r(static_cast<std::uint64_t>(i), static_cast<std::uint64_t>(j));
                }
            }
            
            // Assign dimnames
            SEXP dimnames = guard.Protect(Rf_allocVector(VECSXP, 2));
            SEXP r_names = guard.Protect(Rf_allocVector(STRSXP, nrow));
            
            for (int i = 0; i < nrow; ++i) {
                std::string_view name_view = poset->GetElementName(static_cast<std::size_t>(i));
                SET_STRING_ELT(
                               r_names,
                               static_cast<R_xlen_t>(i),
                               Rf_mkCharLenCE(name_view.data(), static_cast<int>(name_view.size()), CE_UTF8)
                               );
            }
            
            SET_VECTOR_ELT(dimnames, 0, r_names);
            SET_VECTOR_ELT(dimnames, 1, r_names);
            Rf_setAttrib(result_r, R_DimNamesSymbol, dimnames);
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    /**
     * @brief Performs a vectorized check to determine if elements in the first vector are dominated by elements in the second.
     *
     * @details Evaluates the partial order relation element-wise between two input character vectors.
     * For each pair `(v1[k], v2[k])`, it returns `TRUE` if `v1[k] <= v2[k]` (i.e., `v1` is dominated by or equal to `v2`).
     * It leverages internal element IDs (EID) for instantaneous O(1) adjacency matrix lookups.
     *
     * @param poset_r ExternalPtr to the POSetWrap object.
     * @param elements_r_1 Character vector of the first set of elements.
     * @param elements_r_2 Character vector of the second set of elements (must be the same length as v1_r).
     * @return SEXP An R logical vector of length N containing the element-wise dominance results.
     */
    SEXP IsDominatedBy(SEXP poset_r, SEXP elements_r_1, SEXP elements_r_2) {
        SEXP result_r = R_NilValue;
        RProtectGuard guard;
        try {
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);
            POSet* poset = poset_wrap->GetPOSet();
            
            const std::vector<std::string> elements_1 = RConvert::ToStringVector(elements_r_1);
            const std::vector<std::string> elements_2 = RConvert::ToStringVector(elements_r_2);
            
            const int n = static_cast<int>(elements_1.size());
            if (n != static_cast<int>(elements_2.size())) {
                throw MyException("IsDominatedBy: Input vectors must have the same length.");
            }
            
            result_r = guard.Protect(Rf_allocVector(LGLSXP, n));
            // Caching del puntatore per evitare l'overhead della macro LOGICAL()
            int* logical_ptr = LOGICAL(result_r);
            for (int k = 0; k < n; ++k) {
                const auto idx = static_cast<std::size_t>(k);
                
                const std::uint64_t n1 = poset->GetElementId(elements_1[idx]);
                const std::uint64_t n2 = poset->GetElementId(elements_2[idx]);
                
                logical_ptr[k] = (n1 == n2 || poset->GreaterThan(n2, n1)) ? TRUE : FALSE;
            }
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    // ---------------------------------------------------------------------------
    
    /**
     * @brief Checks whether each element of @p v1_r dominates the corresponding
     *        element of @p v2_r.
     *
     * Returns a logical vector where entry @c k is `TRUE` iff \f$v_1[k] \succ v_2[k]\f$.
     *
     * @param poset_r  ExternalPtrAddr of POSetWrap.
     * @param elements_r_1     R character vector of element names.
     * @param elements_r_2     R character vector of element names (same length as @p v1_r).
     * @return R logical vector.
     */
    SEXP Dominates(SEXP poset_r, SEXP elements_r_1, SEXP elements_r_2) {
        SEXP result_r = R_NilValue;
        RProtectGuard guard;
        
        try {
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);
            POSet* poset = poset_wrap->GetPOSet();
            
            const std::vector<std::string> elements_1 = RConvert::ToStringVector(elements_r_1);
            const std::vector<std::string> elements_2 = RConvert::ToStringVector(elements_r_2);
            
            const auto n = static_cast<int>(elements_1.size());
            result_r = guard.Protect(Rf_allocVector(LGLSXP, n));
            // Caching del puntatore per evitare l'overhead della macro LOGICAL()
            int* logical_ptr = LOGICAL(result_r);
            for (int k = 0; k < n; ++k) {
                const auto idx = static_cast<std::size_t>(k);
                
                const std::uint64_t n1 = poset->GetElementId(elements_1[idx]);
                const std::uint64_t n2 = poset->GetElementId(elements_2[idx]);
                
                logical_ptr[k] = (n1 == n2 || poset->GreaterThan(n1, n2)) ? TRUE : FALSE;
            }
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    /**
     * @brief Performs a vectorized check to determine if elements in the first vector are comparable to elements in the second.
     *
     * @details Evaluates the comparability relation element-wise between two input character vectors.
     * For each pair `(v1[k], v2[k])`, it returns `TRUE` if they are the same element, or if one strictly dominates the other
     * (`v1[k] > v2[k]` or `v2[k] > v1[k]`).
     * It leverages internal element IDs (EID) for instantaneous O(1) integer comparisons and matrix lookups.
     *
     * @param poset_r ExternalPtr to the POSetWrap object.
     * @param elements_r_1 Character vector of the first set of elements.
     * @param elements_r_2 Character vector of the second set of elements (must be the same length as v1_r).
     * @return SEXP An R logical vector of length N containing the element-wise comparability results.
     */
    SEXP IsComparableWith(SEXP poset_r, SEXP elements_r_1, SEXP elements_r_2) {
        SEXP result_r = R_NilValue;
        RProtectGuard guard;
        
        try {
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);
            POSet* poset = poset_wrap->GetPOSet();
            
            const std::vector<std::string> elements_1 = RConvert::ToStringVector(elements_r_1);
            const std::vector<std::string> elements_2 = RConvert::ToStringVector(elements_r_2);
            
            const auto n = static_cast<int>(elements_1.size());
            if (n != static_cast<int>(elements_2.size())) {
                throw std::invalid_argument("IsComparableWith: Input vectors must have the same length.");
            }
            
            
            result_r = guard.Protect(Rf_allocVector(LGLSXP, n));
            // Caching del puntatore per evitare l'overhead della macro LOGICAL()
            int* logical_ptr = LOGICAL(result_r);
            for (int k = 0; k < n; ++k) {
                const auto idx = static_cast<std::size_t>(k);
                
                const std::uint64_t n1 = poset->GetElementId(elements_1[idx]);
                const std::uint64_t n2 = poset->GetElementId(elements_2[idx]);
                
                logical_ptr[k] = (n1 == n2 || poset->GreaterThan(n1, n2) || poset->GreaterThan(n2, n1)) ? TRUE : FALSE;
            }
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    /**
     * @brief Performs a vectorized check to determine if elements in the first vector are incomparable to elements in the second.
     *
     * @details Evaluates the incomparability relation element-wise between two input character vectors.
     * For each pair `(v1[k], v2[k])`, it returns `TRUE` if they are distinct elements (`v1[k] != v2[k]`)
     * and neither dominates the other (`!(v1[k] > v2[k])` and `!(v2[k] > v1[k])`).
     * It leverages internal element IDs (EID) for instantaneous O(1) integer comparisons and matrix lookups.
     *
     * @param poset_r ExternalPtr to the POSetWrap object.
     * @param elements_r_1 Character vector of the first set of elements.
     * @param elements_r_2 Character vector of the second set of elements (must be the same length as v1_r).
     * @return SEXP An R logical vector of length N containing the element-wise incomparability results.
     */
    SEXP IsIncomparableWith(SEXP poset_r, SEXP elements_r_1, SEXP elements_r_2) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;
        
        try {
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);
            POSet* poset = poset_wrap->GetPOSet();
            
            const std::vector<std::string> elements_1 = RConvert::ToStringVector(elements_r_1);
            const std::vector<std::string> elements_2 = RConvert::ToStringVector(elements_r_2);

            const int n = static_cast<int>(elements_1.size());
            
            // HPC Safety: previene buffer overflow
            if (n != static_cast<int>(elements_2.size())) {
                throw std::invalid_argument("IsIncomparableWith: Input vectors must have the same length.");
            }
            
            result_r = guard.Protect(Rf_allocVector(LGLSXP, n));
            // HPC Hot-loop: Caching del puntatore per evitare l'overhead della macro LOGICAL() in R
            int* logical_ptr = LOGICAL(result_r);
            
            for (int k = 0; k < n; ++k) {
                const auto idx = static_cast<std::size_t>(k);
                
                const std::uint64_t n1 = poset->GetElementId(elements_1[idx]);
                const std::uint64_t n2 = poset->GetElementId(elements_2[idx]);
                
                /// HPC Logic: due elementi sono incomparabili se sono DIVERSI e NON c'è
                // relazione d'ordine in nessuno dei due versi.
                logical_ptr[k] = (n1 != n2 && !poset->GreaterThan(n1, n2) && !poset->GreaterThan(n2, n1)) ? TRUE : FALSE;
            }
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    /**
     * @brief Checks if a given subset of elements forms an upset in the POSet.
     *
     * @details An upset (or upward closed set) is a subset U such that if x is in U
     * and x <= y, then y is also strictly in U. This function takes a vector of element names
     * and verifies this structural property against the underlying POSet.
     *
     * @param poset_r ExternalPtr to the POSetWrap object.
     * @param elements_r Character vector specifying the subset of elements to check.
     * @return SEXP Logical TRUE if the subset is a valid upset, FALSE otherwise.
     */
    SEXP IsUpset(SEXP poset_r, SEXP elements_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;
        
        try {
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);
            POSet* poset = poset_wrap->GetPOSet();
            
            const R_len_t len = Rf_length(elements_r);
            std::vector<std::uint64_t> elements_id;
            elements_id.reserve(static_cast<std::size_t>(len));
            for (R_len_t i = 0; i < len; ++i) {
                const char* el_cstr = CHAR(STRING_ELT(elements_r, i));
                const std::uint64_t el_id = poset->GetElementId(el_cstr);
                elements_id.emplace_back(el_id);
            }
            const bool r = poset->IsUpSet(elements_id);
            result_r = RCreate::FromBool(guard, r);
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        return result_r;
    }
    
    // ---------------------------------------------------------------------------
    
    /**
     * @brief Computes the combined upset for a given subset of elements.
     *
     * @details The combined upset of a subset is the union of the individual upsets
     * of each element within that subset. This function queries the C++ core to
     * retrieve this aggregated set efficiently.
     *
     * @param poset_r ExternalPtr to the POSetWrap object.
     * @param elements_r Character vector specifying the subset of base elements.
     * @return SEXP A character vector containing the unique elements of the combined upset.
     */
    SEXP UpsetOf(SEXP poset_r, SEXP elements_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;
        
        try {
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);
            POSet* poset = poset_wrap->GetPOSet();

            const R_len_t len = Rf_length(elements_r);
            std::vector<std::uint64_t> elements_id;
            elements_id.reserve(static_cast<std::size_t>(len));
            for (R_len_t i = 0; i < len; ++i) {
                const char* el_cstr = CHAR(STRING_ELT(elements_r, i));
                const std::uint64_t el_id = poset->GetElementId(el_cstr);
                elements_id.emplace_back(el_id);
            }
            auto r = poset->UpSet(elements_id);
            result_r = guard.Protect(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(r.count())));
            R_xlen_t i = 0;
            for (const auto id : r) {
                std::string_view name_view = poset->GetElementName(static_cast<std::size_t>(id));
                SET_STRING_ELT(
                               result_r,
                               i++,
                               Rf_mkCharLenCE(name_view.data(), static_cast<int>(name_view.size()), CE_UTF8)
                               );
            }
            
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    /**
     * @brief Checks if a given subset of elements forms a downset in the POSet.
     *
     * @details A downset (also known as a downward closed set or ideal) is a subset D
     * such that if x is in D and y <= x, then y is also strictly in D. This function
     * takes a vector of element names and verifies this structural property against
     * the underlying POSet.
     *
     * @param poset_r ExternalPtr to the POSetWrap object.
     * @param elements_r Character vector specifying the subset of elements to check.
     * @return SEXP Logical TRUE if the subset is a valid downset, FALSE otherwise.
     */
    SEXP IsDownset(SEXP poset_r, SEXP elements_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;
        
        try {
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);
            POSet* poset = poset_wrap->GetPOSet();
            
            const R_len_t len = Rf_length(elements_r);
            std::vector<std::uint64_t> elements_id;
            elements_id.reserve(static_cast<std::size_t>(len));
            for (R_len_t i = 0; i < len; ++i) {
                const char* el_cstr = CHAR(STRING_ELT(elements_r, i));
                const std::uint64_t el_id = poset->GetElementId(el_cstr);
                elements_id.emplace_back(el_id);
            }
            
            const bool r = poset->IsDownSet(elements_id);
            
            result_r = RCreate::FromBool(guard, r);
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        return result_r;
    }
    
    /**
     * @brief Computes the combined downset for a given subset of elements.
     *
     * @details The combined downset of a subset is the union of the individual downsets
     * of each element within that subset. This function queries the C++ core to
     * retrieve this aggregated set efficiently.
     *
     * @param poset_r ExternalPtr to the POSetWrap object.
     * @param elements_r Character vector specifying the subset of base elements.
     * @return SEXP A character vector containing the unique elements of the combined downset.
     */
    SEXP DownsetOf(SEXP poset_r, SEXP elements_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;
        
        try {
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);
            POSet* poset = poset_wrap->GetPOSet();
            
            const R_len_t len = Rf_length(elements_r);
            
            std::vector<std::uint64_t> elements_id;
            elements_id.reserve(static_cast<std::size_t>(len));
            for (R_len_t i = 0; i < len; ++i) {
                const char* el_cstr = CHAR(STRING_ELT(elements_r, i));
                const std::uint64_t el_id = poset->GetElementId(el_cstr);
                elements_id.emplace_back(el_id);
            }
            auto r = poset->DownSet(elements_id);
            
            result_r = guard.Protect(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(r.count())));
            R_xlen_t i = 0;
            for (const auto id : r) {
                std::string_view name_view = poset->GetElementName(static_cast<std::size_t>(id));
                SET_STRING_ELT(
                               result_r,
                               i++,
                               Rf_mkCharLenCE(name_view.data(), static_cast<int>(name_view.size()), CE_UTF8)
                               );
            }
            
            
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    /**
     * @brief Retrieves the set of elements comparable to a given element.
     *
     * @details Two elements are comparable if one dominates the other (either x <= y or y <= x).
     * This function efficiently queries the underlying POSet using internal integer IDs (EID)
     * to locate all elements comparable to the specified target with minimal overhead.
     *
     * @param poset_r ExternalPtr to the POSetWrap object.
     * @param element_r Character scalar specifying the target element.
     * @return SEXP A character vector containing the names of all comparable elements.
     */
    SEXP ComparabilitySetOf(SEXP poset_r, SEXP element_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;
        
        try {
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);
            POSet* poset = poset_wrap->GetPOSet();
            
            const char* el_cstr = CHAR(STRING_ELT(element_r, 0));

            const auto r = poset->ComparabilitySetOf(poset->GetElementId(el_cstr));

            const int n      = static_cast<int>(r.count());
            result_r = guard.Protect(Rf_allocVector(STRSXP, n));
            int       pos    = 0;
            for (auto id : r) {
                std::string_view name_view = poset->GetElementName(static_cast<std::size_t>(id));
                
                SET_STRING_ELT(
                               result_r,
                               pos++,
                               Rf_mkCharLenCE(name_view.data(), static_cast<int>(name_view.size()), CE_UTF8)
                               );
            }
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    /**
     * @brief Retrieves the set of elements incomparable to a given element.
     *
     * @details Two elements are incomparable if neither dominates the other
     * (i.e., neither x <= y nor y <= x). This function efficiently queries the
     * underlying POSet using internal integer IDs. The conversion from internal
     * IDs to string names is mapped directly into natively allocated R memory,
     * ensuring zero overhead.
     *
     * @param poset_r ExternalPtr to the POSetWrap object.
     * @param element_r Character scalar specifying the target element.
     * @return SEXP A character vector containing the names of all incomparable elements.
     */
    SEXP IncomparabilitySetOf(SEXP poset_r, SEXP element_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;
        
        try {
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);
            POSet* poset = poset_wrap->GetPOSet();
            
            const char* el_cstr = CHAR(STRING_ELT(element_r, 0));
            
            const auto r = poset->IncomparabilitySetOf(poset->GetElementId(el_cstr));
            
            const int n      = static_cast<int>(r.count());
            result_r = guard.Protect(Rf_allocVector(STRSXP, n));
            int       pos    = 0;
            for (auto id : r) {
                std::string_view name_view = poset->GetElementName(static_cast<std::size_t>(id));
                
                SET_STRING_ELT(
                               result_r,
                               pos++,
                               Rf_mkCharLenCE(name_view.data(), static_cast<int>(name_view.size()), CE_UTF8)
                               );
            }
            
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    /**
     * @brief Retrieves all maximal elements within the POSet.
     *
     * @details In order theory, an element is considered maximal if it is not
     * strictly dominated by any other element in the partially ordered set.
     * This function queries the underlying C++ core to efficiently identify
     * all such elements and maps their names directly into natively allocated
     * R memory, ensuring zero overhead.
     *
     * @param poset_r ExternalPtr to the POSetWrap object.
     * @return SEXP A character vector containing the names of all maximal elements.
     */
    SEXP Maximal(SEXP poset_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;
        
        try {
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);
            POSet* poset = poset_wrap->GetPOSet();
            
            const auto r = poset->Maximals();
            
            const int n      = static_cast<int>(r.count());
            result_r = guard.Protect(Rf_allocVector(STRSXP, n));
            int       pos    = 0;
            for (auto id : r) {
                std::string_view name_view = poset->GetElementName(static_cast<std::size_t>(id));
                SET_STRING_ELT(
                               result_r,
                               pos++,
                               Rf_mkCharLenCE(name_view.data(), static_cast<int>(name_view.size()), CE_UTF8)
                               );
            }
            
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    /**
     * @brief Retrieves all minimal elements within the POSet.
     *
     * @details In order theory, an element is considered minimal if no other element
     * in the partially ordered set strictly precedes (or is strictly dominated by) it.
     * This function queries the underlying C++ core to efficiently identify all such
     * elements and maps their names directly into natively allocated R memory,
     * ensuring zero overhead.
     *
     * @param poset_r ExternalPtr to the POSetWrap object.
     * @return SEXP A character vector containing the names of all minimal elements.
     */
    SEXP Minimal(SEXP poset_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;
        
        try {
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);
            POSet* poset = poset_wrap->GetPOSet();
            
            const auto r = poset->Minimals();
            
            const int n      = static_cast<int>(r.count());
            result_r = guard.Protect(Rf_allocVector(STRSXP, n));
            int       pos    = 0;
            for (auto id : r) {
                std::string_view name_view = poset->GetElementName(static_cast<std::size_t>(id));
                SET_STRING_ELT(
                               result_r,
                               pos++,
                               Rf_mkCharLenCE(name_view.data(), static_cast<int>(name_view.size()), CE_UTF8)
                               );
            }
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    /**
     * @brief Checks if a given element is maximal within the POSet.
     *
     * @details In order theory, an element is considered maximal if there is no
     * other element in the partially ordered set that strictly dominates it.
     * This function queries the underlying C++ core to efficiently evaluate this
     * property for the specified element, ensuring zero overhead.
     *
     * @param poset_r ExternalPtr to the POSetWrap object.
     * @param element_r Character scalar specifying the target element to check.
     * @return SEXP Logical TRUE if the element is maximal, FALSE otherwise.
     */
    SEXP IsMaximal(SEXP poset_r, SEXP element_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;
        
        try {
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);
            POSet* poset = poset_wrap->GetPOSet();
            
            const char* el_cstr = CHAR(STRING_ELT(element_r, 0));
            
            const auto r = poset->IsMaximal(poset->GetElementId(el_cstr));
            
            result_r = RCreate::FromBool(guard, r);
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    /**
     * @brief Checks if a given element is minimal within the POSet.
     *
     * @details In order theory, an element is considered minimal if there is no
     * other element in the partially ordered set that strictly precedes it.
     * This function queries the underlying C++ core to efficiently evaluate this
     * property for the specified element, ensuring zero overhead.
     *
     * @param poset_r ExternalPtr to the POSetWrap object.
     * @param element_r Character scalar specifying the target element to check.
     * @return SEXP Logical TRUE if the element is minimal, FALSE otherwise.
     */
    SEXP IsMinimal(SEXP poset_r, SEXP element_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;
        
        try {
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);
            POSet* poset = poset_wrap->GetPOSet();
            
            const char* el_cstr = CHAR(STRING_ELT(element_r, 0));
            
            const auto r = poset->IsMinimal(poset->GetElementId(el_cstr));
            
            result_r = RCreate::FromBool(guard, r);
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    /**
     * @brief Computes the meet (greatest lower bound) of a given subset of elements.
     *
     * @details In order theory, the meet of a subset is the greatest element in the
     * partially ordered set that is less than or equal to all elements in the subset.
     * This function queries the underlying C++ core to efficiently calculate this
     * structural bound for the specified elements.
     *
     * @param poset_r ExternalPtr to the POSetWrap object.
     * @param elements_r Character vector specifying the subset of target elements.
     * @return SEXP A character scalar containing the name of the meet element
     * (or an appropriate empty/NULL value if the meet does not exist).
     */
    SEXP Meet(SEXP poset_r, SEXP elements_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;
        
        try {
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);
            POSet* poset = poset_wrap->GetPOSet();
            
            const R_len_t len = Rf_length(elements_r);
            std::vector<std::uint64_t> elements_id;
            elements_id.reserve(static_cast<std::size_t>(len));
            for (R_len_t i = 0; i < len; ++i) {
                const char* el_cstr = CHAR(STRING_ELT(elements_r, i));
                const std::uint64_t el_id = poset->GetElementId(el_cstr);
                elements_id.emplace_back(el_id);
            }
            auto r = poset->Meet(elements_id);
            
            result_r = guard.Protect(Rf_allocVector(STRSXP, 1));
            
            if (r.has_value()) {
                std::string_view name_view = poset->GetElementName(static_cast<std::size_t>(r.value()));
                SET_STRING_ELT(result_r, 0,
                               Rf_mkCharLenCE(name_view.data(), static_cast<int>(name_view.size()), CE_UTF8));
            } else {
                SET_STRING_ELT(result_r, 0, NA_STRING);
            }
            
            
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    /**
     * @brief Computes the join (least upper bound) of a given subset of elements.
     *
     * @details In order theory, the join of a subset is the least element in the
     * partially ordered set that is greater than or equal to all elements in the subset.
     * This function queries the underlying C++ core to efficiently calculate this
     * structural bound for the specified elements.
     *
     * @param poset_r ExternalPtr to the POSetWrap object.
     * @param elements_r Character vector specifying the subset of target elements.
     * @return SEXP A character scalar containing the name of the join element
     * (or NA_character_ if the join does not exist).
     */
    SEXP Join(SEXP poset_r, SEXP elements_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;
        
        try {
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);
            POSet* poset = poset_wrap->GetPOSet();
            
            const R_len_t len = Rf_length(elements_r);
            std::vector<std::uint64_t> elements_id;
            elements_id.reserve(static_cast<std::size_t>(len));
            for (R_len_t i = 0; i < len; ++i) {
                const char* el_cstr = CHAR(STRING_ELT(elements_r, i));
                const std::uint64_t el_id = poset->GetElementId(el_cstr);
                elements_id.emplace_back(el_id);
            }
            auto r = poset->Join(elements_id);
            
            result_r = guard.Protect(Rf_allocVector(STRSXP, 1));
            if (r.has_value()) {
                std::string_view name_view = poset->GetElementName(static_cast<std::size_t>(r.value()));
                SET_STRING_ELT(result_r, 0,
                               Rf_mkCharLenCE(name_view.data(), static_cast<int>(name_view.size()), CE_UTF8));
            } else {
                SET_STRING_ELT(result_r, 0, NA_STRING);
            }
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    /**
     * @brief Retrieves all pairs of incomparable elements within the POSet.
     *
     * @details In order theory, two elements are considered incomparable if neither
     * dominates the other (i.e., neither x <= y nor y <= x). This function queries
     * the underlying C++ core to efficiently identify all such mutually incomparable
     * pairs across the entire partially ordered set.
     *
     * @param poset_r ExternalPtr to the POSetWrap object.
     * @return SEXP An R data structure (e.g., a matrix or list) containing the
     * names of all incomparable element pairs.
     */
    SEXP Incomparabilities(SEXP poset_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;
        
        try {
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);
            POSet* poset = poset_wrap->GetPOSet();
            
            const auto r = poset->Incomparabilities();
            
            const int n = static_cast<int>(r.size());
            
            result_r = guard.Protect(Rf_allocMatrix(STRSXP, n, 2));
            
            int pos = 0;
            for (const auto& p : r) {
                std::string_view name1 = poset->GetElementName(p.first);
                std::string_view name2 = poset->GetElementName(p.second);
                
                SET_STRING_ELT(result_r, pos,
                               Rf_mkCharLenCE(name1.data(), static_cast<int>(name1.size()), CE_UTF8));
                
                SET_STRING_ELT(result_r, pos + n,
                               Rf_mkCharLenCE(name2.data(), static_cast<int>(name2.size()), CE_UTF8));
                
                pos++;
            }
            
            
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        return result_r;
    }
    
    /**
     * @brief Checks if one POSet is an extension of another.
     *
     * @details In order theory, a partially ordered set P1 is considered an extension
     * of another POSet P2 if they share the same base elements and every order relation
     * present in P2 is strictly preserved in P1 (i.e., if x <= y in P2, then x <= y in P1).
     * This function queries the underlying C++ core to efficiently evaluate this
     * structural relationship between the two specified POSets.
     *
     * @param poset_r_1 ExternalPtr to the POSetWrap object representing the potential extension.
     * @param poset_r_2 ExternalPtr to the POSetWrap object representing the base POSet.
     * @return SEXP Logical TRUE if the first POSet is a valid extension of the second, FALSE otherwise.
     */
    SEXP IsExtensionOf(SEXP poset_r_1, SEXP poset_r_2) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;
        
        try {
            const POSetWrap* poset_wrap_1 = RConvert::ToPOSetWrap(poset_r_1);
            POSet* poset_1 = poset_wrap_1->GetPOSet();

            const POSetWrap* poset_wrap_2 = RConvert::ToPOSetWrap(poset_r_2);
            POSet* poset_2 = poset_wrap_2->GetPOSet();
            
            
            const bool r = poset_1->IsExtensionOf(*poset_2);
            
            result_r = RCreate::FromBool(guard, r);
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
}
