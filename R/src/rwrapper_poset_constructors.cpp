/**
 * @file rwrapper_poset_constructors.cpp
 * @brief Interfaccia R/C++ per i costruttori della libreria PoseticDataAnalysis.
 * * @details Questo modulo implementa i punti di ingresso (entry-points) in stile C
 * (`extern "C"`) invocabili direttamente dalla sessione R tramite `.Call()`.
 * Il design pattern adottato prevede tre fasi per ogni funzione:
 * 1. **Estrazione (RConvert)**: Conversione sicura dei tipi R (`SEXP`) in strutture dati C++20
 * (`std::vector`, `std::string`, puntatori raw) senza copie inutili.
 * 2. **Elaborazione (Core HPC)**: Invocazione dei Factory Method di `POSetWrap` passando i
 * parametri per `const reference` o puntatore osservatore per massimizzare la locality della cache.
 * 3. **Restituzione (RCreate)**: Trasferimento della ownership dei memory-safe pointer C++
 * (`std::unique_ptr`) al Garbage Collector di R, esponendoli come `ExternalPtr`.
 * * @author Alessandro Avellone
 * @version 4.0
 */

#include "poset_wrapper.h"
#include "rwrapper.h"
#include "rwrapper_conversion.h" // Il nuovo header per le utility di conversione

#include <cstdint>
#include <exception>
#include <stdexcept>
#include <string>
#include <vector>
#include <utility>

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>



extern "C" {
    /**
     * @brief Costruisce un POSet generico a partire da un set di nodi e archi.
     * * @param elements_r SEXP di tipo `STRSXP`. Un vettore carattere R contenente i nomi univoci.
     * @param comparabilities_r SEXP matrice di stringhe, rappresentante le relazioni di copertura (a < b).
     * @return SEXP di tipo `EXTPTRSXP`. Puntatore esterno R che incapsula l'istanza `POSetWrap`.
     */
    SEXP BuildPOSet(SEXP elements_r, SEXP comparabilities_r) {
        RProtectGuard guard;
        try {
            auto elements = RConvert::ToStringVector(elements_r);
            auto comparabilities = RConvert::ToPairVector(comparabilities_r);
            
            auto result = POSetWrap::BuildPoset(elements, comparabilities);
            
            return RCreate::WrapExternalPtr(guard, std::move(result));
            
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("C++ exception (unknown reason) in BuildGenericPOSet");
        }
        return R_NilValue;
    }
    
    /**
     * @brief Costruisce un POSet lineare (una Catena o Total Order).
     * * @param elements_r SEXP di tipo `STRSXP`. Sequenza ordinata degli elementi.
     * @return SEXP di tipo `EXTPTRSXP` che incapsula la Catena Lineare.
     */
    SEXP BuildLinearPOSet(SEXP elements_r) {
        RProtectGuard guard;
        try {
            auto elements = RConvert::ToStringVector(elements_r);
            
            auto result = POSetWrap::BuildLinear(elements);
            
            return RCreate::WrapExternalPtr(guard, std::move(result));
            
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        }
        return R_NilValue;
    }
    
    /**
     * @brief Builds the Cartesian product of multiple POSets.
     *
     * @param posets_r  R list of ExternalPtrAddr to POSetWrap objects.
     * @return ExternalPtrAddr of the product POSetWrap.
     */
    SEXP BuildProductPOSet(SEXP posets_r) {
        RProtectGuard guard;
        try {
            auto posets = RConvert::ToPOSetWrapVector(posets_r);
            auto result = POSetWrap::BuildProduct(posets);
            return RCreate::WrapExternalPtr(guard, std::move(result));
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason) in BuildProductPOSet");
        }
        return R_NilValue;
    }
    
    /**
     * @brief Builds the lexicographic product of multiple POSets.
     *
     * @param posets_r  R list of ExternalPtrAddr to POSetWrap objects (order matters).
     * @return ExternalPtrAddr of the lexicographic-product POSetWrap.
     */
    SEXP BuildLexicographicProductPOSet(SEXP posets_r) {
        RProtectGuard guard;
        try {
            auto posets = RConvert::ToPOSetWrapVector(posets_r);
            auto result = POSetWrap::BuildLexicographicProduct(posets);
            return RCreate::WrapExternalPtr(guard, std::move(result));
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason) in BuildLexicographicProductPOSet");
        }
        return R_NilValue;
    }
    
    /**
     * @brief Builds the intersection of multiple POSets.
     *
     * @param posets_r  R list of ExternalPtrAddr to POSetWrap objects.
     * @return ExternalPtrAddr of the intersection POSetWrap.
     */
    SEXP BuildIntersectionPOSet(SEXP posets_r) {
        RProtectGuard guard;
        try {
            auto posets = RConvert::ToPOSetWrapVector(posets_r);
            auto result = POSetWrap::BuildIntersection(posets);
            return RCreate::WrapExternalPtr(guard, std::move(result));
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason) in BuildIntersectionPOSet");
        }
        return R_NilValue;
    }
    
    /**
     * @brief Builds the linear sum of multiple POSets.
     *
     * @param posets_r  R list of ExternalPtrAddr to POSetWrap objects.
     * @return ExternalPtrAddr of the linear-sum POSetWrap.
     */
    SEXP BuildLinearSumPOSet(SEXP posets_r) {
        RProtectGuard guard;
        try {
            auto posets = RConvert::ToPOSetWrapVector(posets_r);
            auto result = POSetWrap::BuildLinearSum(posets);
            return RCreate::WrapExternalPtr(guard, std::move(result));
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason) in BuildLinearSumPOSet");
        }
        return R_NilValue;
    }
    
    /**
     * @brief Builds the disjoint sum (ordinal sum with no cross-relations) of multiple POSets.
     *
     * @param posets_r  R list of ExternalPtrAddr to POSetWrap objects.
     * @return ExternalPtrAddr of the disjoint-sum POSetWrap.
     */
    SEXP BuildDisjointSumPOSet(SEXP posets_r) {
        RProtectGuard guard;
        try {
            auto posets = RConvert::ToPOSetWrapVector(posets_r);
            auto result = POSetWrap::BuildDisjointSum(posets);
            return RCreate::WrapExternalPtr(guard, std::move(result));
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason) in BuildDisjointSumPOSet");
        }
        return R_NilValue;
    }
    
    /**
     * @brief Builds a lifting POSet by adding a new minimum element.
     *
     * @param poset_r       ExternalPtrAddr of the base POSetWrap.
     * @param new_element_r R character scalar.
     * @return ExternalPtrAddr of the lifting POSetWrap.
     */
    SEXP BuildLiftingPOSet(SEXP poset_r, SEXP new_element_r) {
        RProtectGuard guard;
        try {
            const auto* base_poset = RConvert::ToPOSetWrap(poset_r);
            const std::string new_element = CHAR(STRING_ELT(new_element_r, 0));
            
            auto result = POSetWrap::BuildLiftingPOSet(base_poset, new_element);
            return RCreate::WrapExternalPtr(guard, std::move(result));
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason) in BuildLiftingPOSet");
        }
        return R_NilValue;
    }
    
    /**
     * @brief Builds a bucket POSet from elements and comparabilities.
     */
    SEXP BuildBucketPOSet(SEXP elements_r, SEXP comparabilities_r) {
        RProtectGuard guard;
        try {
            auto elements = RConvert::ToStringVector(elements_r);
            auto comparabilities = RConvert::ToPairVector(comparabilities_r);
            
            auto result = POSetWrap::BuildBucketPOSet(elements, comparabilities);
            return RCreate::WrapExternalPtr(guard, std::move(result));
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason) in BuildBucketPOSet");
        }
        return R_NilValue;
    }
    
    
    /**
     * @brief Builds a binary-variable POSet.
     */
    SEXP BuildBinaryVariablePOSet(SEXP variables_r) {
        RProtectGuard guard;
        try {
            auto variables = RConvert::ToStringVector(variables_r);
            
            auto result = POSetWrap::BuildBinaryVariablePOSet(variables);
            return RCreate::WrapExternalPtr(guard, std::move(result));
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason) in BuildBinaryVariablePOSet");
        }
        return R_NilValue;
    }
    
    // ---------------------------------------------------------------------------
    
    /**
     * @brief Builds a fence (zigzag) POSet.
     */
    SEXP BuildFencePOSet(SEXP elements_r, SEXP orientation_r) {
        RProtectGuard guard;
        try {
            auto elements = RConvert::ToStringVector(elements_r);
            // Sfruttiamo LOGICAL di R internamente o la futura RConvert::ToBool
            bool orientation = LOGICAL(orientation_r)[0];
            
            auto result = POSetWrap::BuildFencePOSet(elements, orientation);
            return RCreate::WrapExternalPtr(guard, std::move(result));
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason) in BuildFencePOSet");
        }
        return R_NilValue;
    }
    
    /**
     * @brief Builds a Crown POSet \f$S_n^0\f$.
     */
    SEXP BuildCrownPOSet(SEXP elements_1_r, SEXP elements_2_r) {
        RProtectGuard guard;
        try {
            auto elements_1 = RConvert::ToStringVector(elements_1_r);
            auto elements_2 = RConvert::ToStringVector(elements_2_r);
            const auto n = elements_1.size();
            
            // Alloco un vettore standard invece di fare uno std::shared_ptr overhead
            std::vector<std::string> all_elements;
            all_elements.reserve(2 * n);
            for (std::size_t k = 0; k < n; ++k) {
                all_elements.push_back(elements_1[k]);
                all_elements.push_back(elements_2[k]);
            }
            
            std::vector<std::pair<std::string, std::string>> comparabilities;
            comparabilities.reserve(n * (n - 1)); // Riserva per performance
            for (std::size_t k = 0; k < n; ++k) {
                for (std::size_t h = 0; h < n; ++h) {
                    if (k != h) {
                        comparabilities.emplace_back(elements_1[k], elements_2[h]);
                    }
                }
            }
            
            auto result = POSetWrap::BuildPoset(all_elements, comparabilities);
            return RCreate::WrapExternalPtr(guard, std::move(result));
            
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason) in BuildCrownPOSet");
        }
        return R_NilValue;
    }
    /**
     * @brief Builds the dual (order-reversed) POSet.
     */
    SEXP BuildDualPOSet(SEXP poset_r) {
        RProtectGuard guard;
        try {
            const auto* base_poset = RConvert::ToPOSetWrap(poset_r);
            
            auto result = POSetWrap::BuildDualPOSet(base_poset);
            return RCreate::WrapExternalPtr(guard, std::move(result));
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason) in BuildDualPOSet");
        }
        return R_NilValue;
    }
}
