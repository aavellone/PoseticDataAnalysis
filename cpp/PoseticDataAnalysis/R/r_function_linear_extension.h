/**
 * @file r_function_linear_extension.h
 * @brief Decoupled interface between C++ linear extension evaluation and user-defined R functions.
 *
 * @details
 * This header defines @c FLERInterface, a concrete subclass of
 * @c FunctionLinearExtension that delegates the evaluation of each linear
 * extension to a user-supplied R function.
 *
 * It is heavily optimized for High-Performance Computing (HPC), leveraging
 * C++20 features. It introduces **R String Caching (CHARSXP)** to completely
 * avoid `Rf_mkChar` hash-table lookups during the hot evaluation loop.
 *
 * @par Requirements
 * C++20 or later (`-std=c++20`)
 */

#pragma once

#include <memory>
#include <string>
#include <vector>
#include <cstdint>

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif

#include <R.h>
#include <Rinternals.h>

// Inclusione corretta dei file del tuo GEM
#include "function_linear_extension.h"
#include "linear_extension.h"
#include "poset.h"
#include "rwrapper.h"

/**
 * @class FLERInterface
 * @brief Concrete implementation of FunctionLinearExtension evaluating via an R function.
 */
class FLERInterface final : public FunctionLinearExtension {
private:
    /**
     * @brief The R function object (SEXP) supplied by the user.
     */
    SEXP r_function_;
    
    /**
     * @brief HPC Optimization: Pre-allocated R character vector containing element names.
     * By caching the R strings (CHARSXP) and pinning them in memory, we avoid calling
     * Rf_mkChar during the heavily repeated evaluation loop.
     */
    SEXP r_element_names_;
    
    SEXP r_le_buffer_;
    SEXP r_fcall_prebuilt_;

    /**
     * @brief Names of the rows inferred from the R matrix.
     */
    std::vector<std::string> rows_name_;
    
    /**
     * @brief Names of the columns inferred from the R matrix.
     */
    std::vector<std::string> cols_name_;
    
    /**
     * @brief Helper to call the R function and perform basic validation safely.
     */
    SEXP CallRFunction(SEXP le_r, RProtectGuard& guard);
    
public:
    /**
     * @brief Constructs a FLERInterface, pre-allocates R strings, and probes the output.
     *
     * @param poset     A sample linear extension used to probe the output shape.
     * @param r_function    The user-defined R function.
     *
     * @throws std::runtime_error if the R function does not return a valid matrix.
     */
    FLERInterface(const POSet* poset,
                  SEXP r_function);

    /**
     * @brief Non clonabile: possiede oggetti R (SEXP) rilasciati nel distruttore;
     * una copia condivisa causerebbe un doppio rilascio al GC di R.
     * @throws MyException sempre.
     */
    [[nodiscard]] std::unique_ptr<FunctionLinearExtension> Clone() const override {
        throw MyException("FLERInterface::Clone: R-based functions cannot be cloned.");
    }

    /**
     * @brief Non thread-safe: invoca l'API di R, consentita solo sul thread principale.
     */
    [[nodiscard]] bool IsThreadSafe() const noexcept override { return false; }
    
    /**
     * @brief Destructor.
     * @details Releases the preserved R memory block (r_element_names_) to the
     * R Garbage Collector.
     */
    ~FLERInterface() override;
    
    /**
     * @brief Evaluates the R function on a given linear extension.
     *
     * @details
     * Converts @p x to an R character vector using O(1) pointer copies, calls
     * the stored R function, and populates the HPC Structure of Arrays (SoA).
     * Iterates in column-major order to match R's memory layout.
     *
     * @param x Const reference to the linear extension to evaluate (matches base class).
     */
    void operator()(const LinearExtension& x) noexcept override;
     
    std::string_view GetRowNameAt(std::size_t k) const override;
    
    std::string_view GetColNameAt(std::size_t k) const override;
    
};
