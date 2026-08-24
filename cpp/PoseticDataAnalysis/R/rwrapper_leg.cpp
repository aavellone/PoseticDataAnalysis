/**
 * @file rwrapper_leg.cpp
 * @brief R/C++ interface — Linear Extension Generator (LEG) functions.
 *
 * Implements all R-callable (`extern "C"`) entry points for building and
 * querying Linear Extension Generators:
 * - Exact enumeration via Tree-of-Ideals (`BuildLEGenerator` / `LEGGet`).
 * - Approximate MCMC sampling via Bubley-Dyer (`BuildBubleyDyerLEGenerator` /
 * `LEGBubleyDyerGet`).
 *
 * @author Alessandro Avellone
 * @version 2.1
 * @date 2025
 */
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>

#include "rwrapper.h"
#include "my_exception.h"
#include "random.h"
#include "display_message.h"
#include "linear_generator_wrapper.h"
#include "poset_wrapper.h"
#include "rwrapper_conversion.h"
#include "r_display.h"
#include "poset.h"

// ***********************************************
// Linear Extension Generator Functions
// ***********************************************

extern "C" {
    /**
     * @brief Creates an exact Linear Extension Generator (Tree-of-Ideals).
     *
     * The returned opaque pointer is registered with an R finalizer so the
     * underlying C++ object is freed when the R object is garbage-collected.
     *
     * @param poset_r  ExternalPtrAddr of a POSetWrap.
     * @return ExternalPtrAddr wrapping a heap-allocated LinearGeneratorWrap.
     */
    SEXP BuildLEGenerator(SEXP poset_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;
        
        try {
            // 1. Estrazione del raw pointer
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);
            
            // 2. Costruzione dell'oggetto C++
            auto generator = LinearGeneratorWrap::BuildLEGenerator(poset_wrap);
            
            // 3. Creazione del puntatore esterno in R
            result_r = RCreate::WrapExternalPtr(guard, std::move(generator));
            
            // 4. Creazione dello "scudo" per il Garbage Collector
            R_SetExternalPtrProtected(result_r, poset_r);
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    /**
     * @brief Creates a Bubley-Dyer MCMC Linear Extension Generator.
     *
     * @details This function initializes an approximate linear extension generator
     * based on the Bubley-Dyer Markov Chain Monte Carlo (MCMC) algorithm.
     * The returned opaque pointer is safely registered with an R finalizer, ensuring
     * that the underlying C++ object is properly deallocated by R's garbage collector.
     *
     * @param poset_r ExternalPtr to the underlying POSetWrap object.
     * @param seed_r Integer scalar specifying the seed for the random number generator.
     * @return SEXP An ExternalPtr wrapping the heap-allocated LinearGeneratorWrap.
     */
    SEXP BuildBubleyDyerLEGenerator(SEXP poset_r, SEXP seed_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;
        
        try {
            // 1. Estrazione del raw pointer
            const POSetWrap* poset_wrap = RConvert::ToPOSetWrap(poset_r);
            
            // 2. Estrazione del seed
            std::optional<std::uint64_t> seed = RConvert::ToOptionalSeed(seed_r);
            if (!seed.has_value()) {
                seed = Random::GLOBAL.RndNextInt(0, std::numeric_limits<std::uint64_t>::max());
            }
            
            // 3. Costruzione dell'oggetto C++
            auto generator = LinearGeneratorWrap::BuildBubleyDyerGenerator(poset_wrap, seed.value());

            // 4. Creazione del puntatore esterno in R
            SEXP ptr_r = RCreate::WrapExternalPtr(guard, std::move(generator));

            // 5. Creazione dello "scudo" per il Garbage Collector
            R_SetExternalPtrProtected(ptr_r, poset_r);

            // 6. Il seme usato viene restituito a R come stringa: reimmetterlo
            // riproduce esattamente la stessa sequenza di estensioni lineari.
            result_r = RCreate::GeneratorWithSeed(guard, ptr_r, seed.value());
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    /**
     * @brief Draws linear extensions from a Bubley-Dyer generator.
     *
     * The number of extensions returned depends on the stopping criterion:
     * - If @p quante_r is non-empty: exactly `quante` extensions are generated.
     * - Otherwise the required count is computed as
     * \f$ N^4 \cdot (\ln N)^2 + N^3 \cdot \ln N \cdot \ln(\varepsilon^{-1}) \f$
     * where \f$N\f$ is the number of poset elements and
     * \f$\varepsilon\f$ = value of @p errore_r (0 if absent).
     *
     * Setting @p from_start_r to `TRUE` resets the chain to its initial state
     * (starting from the minimal elements).
     *
     * @param generator_r          ExternalPtrAddr of a LinearGeneratorWrap (Bubley-Dyer).
     * @param from_start_r         Logical scalar: reset the chain before sampling.
     * @param quante_r             Integer scalar (fixed sample size), or empty.
     * @param errore_r             Real scalar (target TV distance from uniform), or empty.
     * @param output_ogni_in_sec_r Integer scalar (progress interval in seconds), or empty.
     * @return An (LESize × k) R character matrix; each column is one linear extension.
     */
    SEXP LEGBubleyDyerGet(SEXP generator_r, SEXP from_start_r,
                          SEXP quante_r, SEXP errore_r,
                          SEXP output_ogni_in_sec_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;
        
        try {
            
            auto* generator = RConvert::ToWrapPtr<LinearGeneratorWrap>(
                generator_r, "Il generatore di estensioni lineari");
            
            const bool from_start = RConvert::ToBool(from_start_r);
            const std::optional<std::uint64_t> quante = RConvert::ToOptionalUInt(quante_r);
            const std::optional<double> errore = RConvert::ToOptionalDouble(errore_r);
            const std::optional<std::uint64_t> output_ogni = RConvert::ToOptionalUInt(output_ogni_in_sec_r);
            
            const auto les = generator->GetFromBubleyDyer(from_start, quante, errore, output_ogni);
             
            const int nrow = static_cast<int>(generator->LESize());
            const int ncol = static_cast<int>(les.size());
            result_r = guard.Protect(Rf_allocMatrix(STRSXP, nrow, ncol));
            
            int col = 0;
            for (const auto& le : les) {
                for (int row = 0; row < nrow; ++row)
                    SET_STRING_ELT(result_r, row + nrow * col,
                                   Rf_mkChar(le[static_cast<std::size_t>(row)].c_str()));
                ++col;
            }
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
    
    /**
     * @brief Draws linear extensions from an exact (Tree-of-Ideals) generator.
     *
     * - If @p quante_r is non-empty: up to `quante` not-yet-returned extensions
     * are generated (fewer are returned if the enumeration is exhausted).
     * - If @p quante_r is empty: all remaining extensions are returned.
     *
     * Setting @p from_start_r to `TRUE` resets the generator to the beginning.
     *
     * @param generator_r          ExternalPtrAddr of a LinearGeneratorWrap (Tree-of-Ideals).
     * @param from_start_r         Logical scalar: restart enumeration from the beginning.
     * @param quante_r             Integer scalar (max extensions to return), or empty.
     * @param output_ogni_in_sec_r Integer scalar (progress interval in seconds), or empty.
     * @return An (LESize × k) R character matrix; each column is one linear extension.
     */
    SEXP LEGGet(SEXP generator_r, SEXP from_start_r,
                SEXP quante_r, SEXP output_ogni_in_sec_r) {
        RProtectGuard guard;
        SEXP result_r = R_NilValue;
        
        try {
            auto* generator = RConvert::ToWrapPtr<LinearGeneratorWrap>(
                generator_r, "Il generatore di estensioni lineari");
            
            const bool from_start = RConvert::ToBool(from_start_r);
            const std::optional<std::uint64_t> quante = RConvert::ToOptionalUInt(quante_r);
            const std::optional<std::uint64_t> output_ogni = RConvert::ToOptionalUInt(output_ogni_in_sec_r);

            auto les = generator->GetFromLE(from_start, quante, output_ogni);
            
            const int nrow = static_cast<int>(generator->LESize());
            const int ncol = static_cast<int>(les.size());
            result_r = guard.Protect(Rf_allocMatrix(STRSXP, nrow, ncol));
            
            int col = 0;
            for (const auto& le : les) {
                for (int row = 0; row < nrow; ++row)
                    SET_STRING_ELT(result_r, row + nrow * col,
                                   Rf_mkChar(le[static_cast<std::size_t>(row)].c_str()));
                ++col;
            }
        } catch (const std::exception& ex) {
            forward_exception_to_r(ex.what());
        } catch (...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        
        return result_r;
    }
}
