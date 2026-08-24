/**
 * @file rwrapper_relations.cpp
 * @brief R/C++ interface for testing binary relation properties.
 *
 * @details This file implements the R-callable (extern "C") entry points.
 * It bridges R strings/matrices to high-performance C++ matrix structures
 * and computes mathematical properties (reflexivity, transitivity, etc.)
 * using algorithms optimized for cache locality.
 *
 * @author Alessandro Avellone
 * @version 4.0
 */

#include "generic_functions.h"
#include "tensor.h"
#include "rwrapper_conversion.h"
#include "rwrapper.h"

#include <cstdint>
#include <exception>
#include <list>
#include <string>
#include <utility>
#include <vector>

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>

extern "C" {
    
    // ===========================================================================
    // 2-PARAMETER FUNCTIONS (elements_r, comparabilities_r)
    // ===========================================================================
    
    /**
     * @brief Checks if a binary relation is reflexive.
     *
     * @param elements_r SEXP object representing the set of elements (character vector).
     * @param comparabilities_r SEXP object representing edges (Nx2 character matrix).
     * @return SEXP A logical scalar in R indicating if the relation is reflexive.
     */
    SEXP isReflexive(SEXP elements_r, SEXP comparabilities_r) {
        SEXP result_r = R_NilValue;
        RProtectGuard guard;
        try {
            std::vector<std::string> elements;
            auto adj_matrix = RConvert::ToAdjacencyMatrix(elements_r,
                                                          comparabilities_r,
                                                          elements);
            
            const bool is_property_met = generic::IsReflexive(elements.size(),
                                                              adj_matrix);
            result_r = RCreate::FromBool(guard, is_property_met);
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("C++ exception (unknown reason)");
        }
        return result_r;
    }
    
    /**
     * @brief Checks if a binary relation is a preorder (reflexive and transitive).
     *
     * @param elements_r SEXP object representing the set of elements.
     * @param comparabilities_r SEXP object representing edges.
     * @return SEXP A logical scalar in R indicating if the relation is a preorder.
     */
    SEXP isPreorder(SEXP elements_r, SEXP comparabilities_r) {
        SEXP result_r = R_NilValue;
        RProtectGuard guard;
        try {
            std::vector<std::string> elements;
            auto adj_matrix = RConvert::ToAdjacencyMatrix(elements_r,
                                                          comparabilities_r,
                                                          elements);
            
            const bool is_property_met = generic::IsPreorder(elements.size(),
                                                             adj_matrix);
            result_r = RCreate::FromBool(guard, is_property_met);
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("C++ exception (unknown reason)");
        }
        return result_r;
    }
    
    /**
     * @brief Checks if a binary relation is a partial order.
     *
     * @details A partial order requires the relation to be reflexive, antisymmetric,
     * and transitive.
     *
     * @param elements_r SEXP object representing the set of elements.
     * @param comparabilities_r SEXP object representing edges.
     * @return SEXP A logical scalar in R indicating if it is a partial order.
     */
    SEXP isPartialOrder(SEXP elements_r, SEXP comparabilities_r) {
        SEXP result_r = R_NilValue;
        RProtectGuard guard;
        try {
            std::vector<std::string> elements;
            auto adj_matrix = RConvert::ToAdjacencyMatrix(elements_r,
                                                          comparabilities_r,
                                                          elements);
            
            const bool is_property_met = generic::IsPartialOrder(elements.size(),
                                                                 adj_matrix);
            result_r = RCreate::FromBool(guard, is_property_met);
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("C++ exception (unknown reason)");
        }
        return result_r;
    }
    
    /**
     * @brief Computes the reflexive closure of a binary relation.
     *
     * @param elements_r SEXP object representing the set of elements.
     * @param comparabilities_r SEXP object representing edges.
     * @return SEXP A pair list (or suitable R structure) containing the new edges.
     */
    SEXP ReflexiveClosure(SEXP elements_r, SEXP comparabilities_r) {
        SEXP result_r = R_NilValue;
        RProtectGuard guard;
        try {
            std::vector<std::string> elements;
            auto adj_matrix = RConvert::ToAdjacencyMatrix(elements_r,
                                                          comparabilities_r,
                                                          elements);
            const std::uint64_t num_nodes = elements.size();
            
            generic::ReflexiveClosureInPlace(num_nodes, adj_matrix);
            
            std::list<std::pair<std::string, std::string>> edges;
            for (std::uint64_t i = 0; i < num_nodes; ++i) {
                for (std::uint64_t j = 0; j < num_nodes; ++j) {
                    if (adj_matrix(i, j) == 1) {
                        edges.emplace_back(elements[i], elements[j]);
                    }
                }
            }
            result_r = RCreate::FromStringPairList(guard, edges);
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("C++ exception (unknown reason)");
        }
        return result_r;
    }
    
    // ===========================================================================
    // 1-PARAMETER FUNCTIONS (comparabilities_r only)
    // ===========================================================================
    
    /**
     * @brief Checks if a binary relation is symmetric.
     *
     * @note This function infers the element set dynamically from the edge list,
     * as isolated nodes do not affect symmetry.
     *
     * @param comparabilities_r SEXP object representing edges (Nx2 character matrix).
     * @return SEXP A logical scalar in R indicating if the relation is symmetric.
     */
    SEXP isSymmetric(SEXP comparabilities_r) {
        SEXP result_r = R_NilValue;
        RProtectGuard guard;
        try {
            std::vector<std::string> elements;
            auto adj_matrix = RConvert::ToAdjacencyMatrixFromEdges(comparabilities_r,
                                                                   elements);
            
            const bool is_property_met = generic::IsSymmetric(elements.size(),
                                                              adj_matrix);
            result_r = RCreate::FromBool(guard, is_property_met);
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("C++ exception (unknown reason)");
        }
        return result_r;
    }
    
    /**
     * @brief Checks if a binary relation is antisymmetric.
     *
     * @param comparabilities_r SEXP object representing edges.
     * @return SEXP A logical scalar in R indicating if it is antisymmetric.
     */
    SEXP isAntisymmetric(SEXP comparabilities_r) {
        SEXP result_r = R_NilValue;
        RProtectGuard guard;
        try {
            std::vector<std::string> elements;
            auto adj_matrix = RConvert::ToAdjacencyMatrixFromEdges(comparabilities_r,
                                                                   elements);
            
            const bool is_property_met = generic::IsAntisymmetric(elements.size(),
                                                                  adj_matrix);
            result_r = RCreate::FromBool(guard, is_property_met);
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("C++ exception (unknown reason)");
        }
        return result_r;
    }
    
    /**
     * @brief Checks if a binary relation is transitive.
     *
     * @param comparabilities_r SEXP object representing edges.
     * @return SEXP A logical scalar in R indicating if the relation is transitive.
     */
    SEXP isTransitive(SEXP comparabilities_r) {
        SEXP result_r = R_NilValue;
        RProtectGuard guard;
        try {
            std::vector<std::string> elements;
            auto adj_matrix = RConvert::ToAdjacencyMatrixFromEdges(comparabilities_r,
                                                                   elements);
            
            const bool is_property_met = generic::IsTransitive(elements.size(),
                                                               adj_matrix);
            result_r = RCreate::FromBool(guard, is_property_met);
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("C++ exception (unknown reason)");
        }
        return result_r;
    }
    
    /**
     * @brief Computes the transitive closure of a binary relation.
     *
     * @param comparabilities_r SEXP object representing edges.
     * @return SEXP A pair list (or suitable R structure) containing the new edges.
     */
    SEXP TransitiveClosure(SEXP comparabilities_r) {
        SEXP result_r = R_NilValue;
        RProtectGuard guard;
        try {
            std::vector<std::string> elements;
            auto adj_matrix = RConvert::ToAdjacencyMatrixFromEdges(comparabilities_r, elements);
            const std::uint64_t num_nodes = elements.size();
            
            generic::TransitiveClosureInPlace(num_nodes, adj_matrix);
            
            std::list<std::pair<std::string, std::string>> edges;
            for (std::uint64_t i = 0; i < num_nodes; ++i) {
                for (std::uint64_t j = 0; j < num_nodes; ++j) {
                    if (adj_matrix(i, j) == 1) {
                        edges.emplace_back(elements[i], elements[j]);
                    }
                }
            }
            result_r = RCreate::FromStringPairList(guard, edges);
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("C++ exception (unknown reason)");
        }
        return result_r;
    }
    
}  // extern "C"
