/**
 * @file rFunctionLinearExtension.cpp
 * @brief Implementation of FLERInterface.
 *
 * @details
 * This file provides the concrete implementation of @c FLERInterface, completely
 * decoupled from the POSet class to maximize single-responsibility and
 * testability. It maps C++ linear extensions to R character vectors in O(1) time
 * by caching CHARSXP pointers, and processes the resulting matrices into
 * cache-friendly, SoA (Structure of Arrays) flat vectors.
 */

#include "r_function_linear_extension.h"

#include <format>
#include <stdexcept>
#include <string>
#include <cstdint>

#include "my_exception.h"
#include "poset.h"


// ============================================================
//  FLERInterface – Implementation
// ============================================================

SEXP FLERInterface::CallRFunction(SEXP le_r, RProtectGuard& guard) {
    auto r_fcall = guard.Protect(Rf_lang2(r_function_, le_r));
    
    SEXP ans = R_NilValue;
    try {
        ans = guard.Protect(Rf_eval(r_fcall, R_GlobalEnv));
    } catch (const std::exception& e) {
        throw MyException(std::format("R evaluation error: {}", e.what()));
    } catch (...) {
        throw MyException("Unknown R evaluation error");
    }
    
    if (!Rf_isMatrix(ans)) [[unlikely]] {
        throw MyException("The provided R function did not return a matrix.");
    }
    return ans;
}

FLERInterface::FLERInterface(const POSet* poset,
                             SEXP r_function)
: FunctionLinearExtension(),
r_function_(r_function) {
    
    // Allocate R string array and preserve it from GC
    const R_xlen_t num_elements = static_cast<R_xlen_t>(poset->size());
    r_element_names_ = Rf_allocVector(STRSXP, num_elements);
    R_PreserveObject(r_element_names_);
    
    for (R_xlen_t i = 0; i < num_elements; ++i) {
        std::string_view name_view = poset->GetElementName(static_cast<std::size_t>(i));
        
        SET_STRING_ELT(
                       r_element_names_,
                       i,
                       Rf_mkCharLenCE(name_view.data(), static_cast<int>(name_view.size()), CE_UTF8)
                       );
    }
    
    
    
    // Pre-alloca il vettore di input che verrà passato a R
    r_le_buffer_ = Rf_allocVector(STRSXP, num_elements);
    R_PreserveObject(r_le_buffer_);
    
    // Pre-costruisce la chiamata R (la funzione punterà sempre a r_le_buffer_)
    r_fcall_prebuilt_ = Rf_lang2(r_function_, r_le_buffer_);
    R_PreserveObject(r_fcall_prebuilt_);
    
    // Evaluate the R function using the sample extension to determine
    // the matrix shape and the row/column names ahead of time.
    LinearExtension sample_le(poset->size());
    poset->FirstLE(sample_le);
    RProtectGuard guard;
    
    const R_xlen_t len = static_cast<R_xlen_t>(sample_le.size());
    for (R_xlen_t k = 0; k < len; ++k) {
        const std::uint64_t element_id = sample_le.GetVal(static_cast<std::uint64_t>(k));
        SET_STRING_ELT(r_le_buffer_, k, STRING_ELT(r_element_names_, static_cast<R_xlen_t>(element_id)));
    }
    
    // Valuta la chiamata pre-costruita
    SEXP ans = guard.Protect(Rf_eval(r_fcall_prebuilt_, R_GlobalEnv));
    if (!Rf_isMatrix(ans)) {
        throw MyException("The provided R function did not return a matrix.");
    }
    
    // R matrix dimensions are int; convert to std::uint64_t for STL/HPC use.
    auto ans_dim = Rf_getAttrib(ans, R_DimSymbol);
    const std::uint32_t ans_nrow = static_cast<std::uint32_t>(INTEGER(ans_dim)[0]);
    const std::uint32_t ans_ncol = static_cast<std::uint32_t>(INTEGER(ans_dim)[1]);
    
    shape_ = {ans_nrow, ans_ncol};
    rows_name_.resize(static_cast<std::size_t>(ans_nrow));
    cols_name_.resize(static_cast<std::size_t>(ans_ncol));
    
    auto ans_dim_names = Rf_getAttrib(ans, R_DimNamesSymbol);
    if (ans_dim_names == R_NilValue) {
        for (std::uint64_t k = 0; k < ans_nrow; ++k) {
            rows_name_[static_cast<std::size_t>(k)] = std::to_string(k + 1);
        }
        for (std::uint64_t k = 0; k < ans_ncol; ++k) {
            cols_name_[static_cast<std::size_t>(k)] = std::to_string(k + 1);
        }
    } else {
        auto ans_dim_names_row = VECTOR_ELT(ans_dim_names, 0);
        for (std::uint64_t k = 0; k < ans_nrow; ++k) {
            rows_name_[static_cast<std::size_t>(k)] =
            CHAR(STRING_ELT(ans_dim_names_row, static_cast<R_xlen_t>(k)));
        }
        
        auto ans_dim_names_col = VECTOR_ELT(ans_dim_names, 1);
        for (std::uint64_t k = 0; k < ans_ncol; ++k) {
            cols_name_[static_cast<std::size_t>(k)] =
            CHAR(STRING_ELT(ans_dim_names_col, static_cast<R_xlen_t>(k)));
        }
    }
    
    // Clear base class Structure of Arrays (SoA) before starting batch processing
    idx0_.clear();
    idx1_.clear();
    values_.clear();
}


FLERInterface::~FLERInterface() {
    if (r_element_names_ != R_NilValue) R_ReleaseObject(r_element_names_);
    if (r_le_buffer_ != R_NilValue)     R_ReleaseObject(r_le_buffer_);
    if (r_fcall_prebuilt_ != R_NilValue) R_ReleaseObject(r_fcall_prebuilt_);
}

void FLERInterface::operator()(const LinearExtension& x) noexcept {
    ++calls_;
    
    try {
        RProtectGuard guard;
        
        // Aggiorna IN-PLACE il buffer esistente (Zero Allocazioni in R!)
        const R_xlen_t len = static_cast<R_xlen_t>(x.size());
        for (R_xlen_t k = 0; k < len; ++k) {
            const std::uint64_t element_id = x.GetVal(static_cast<std::uint64_t>(k));
            SET_STRING_ELT(r_le_buffer_, k, STRING_ELT(r_element_names_, static_cast<R_xlen_t>(element_id)));
        }
        
        // 2. Valuta la chiamata pre-costruita
        // Siccome r_le_buffer_ è stato modificato, r_fcall_prebuilt_ userà automaticamente i nuovi dati
        SEXP ans = guard.Protect(Rf_eval(r_fcall_prebuilt_, R_GlobalEnv));
        
        
        const std::uint64_t ans_nrow = shape_[0];
        const std::uint64_t ans_ncol = shape_[1];
        
        // HPC Optimization: Pre-allocate exact memory space required
        const std::size_t total_elements = static_cast<std::size_t>(ans_nrow * ans_ncol);
        idx0_.clear();
        idx1_.clear();
        values_.clear();
        
        idx0_.reserve(total_elements);
        idx1_.reserve(total_elements);
        values_.reserve(total_elements);
        
        // Rebuild the SoA vectors iterating in column-major order to match R's memory layout
        if (Rf_isReal(ans)) {
            const double* data_ptr = REAL(ans);
            for (std::uint32_t col = 0; col < ans_ncol; ++col) {
                for (std::uint32_t row = 0; row < ans_nrow; ++row) {
                    const R_xlen_t idx = static_cast<R_xlen_t>(row + ans_nrow * col);
                    idx0_.push_back(row);
                    idx1_.push_back(col);
                    values_.push_back(data_ptr[idx]);
                }
            }
        } else if (Rf_isInteger(ans)) {
            const int* data_ptr = INTEGER(ans);
            for (std::uint32_t col = 0; col < ans_ncol; ++col) {
                for (std::uint32_t row = 0; row < ans_nrow; ++row) {
                    const R_xlen_t idx = static_cast<R_xlen_t>(row + ans_nrow * col);
                    idx0_.push_back(row);
                    idx1_.push_back(col);
                    values_.push_back(static_cast<double>(data_ptr[idx]));
                }
            }
        } else if (Rf_isLogical(ans)) {
            const int* data_ptr = LOGICAL(ans);
            for (std::uint32_t col = 0; col < ans_ncol; ++col) {
                for (std::uint32_t row = 0; row < ans_nrow; ++row) {
                    const R_xlen_t idx = static_cast<R_xlen_t>(row + ans_nrow * col);
                    idx0_.push_back(row);
                    idx1_.push_back(col);
                    values_.push_back(static_cast<double>(data_ptr[idx]));
                }
            }
        } else [[unlikely]] {
            throw MyException("The R function evaluated a non-numerical matrix type.");
        }
        
    } catch (const std::exception& e) {
        Rf_error("Exception in FLERInterface during evaluation: %s", e.what());
    } catch (...) {
        Rf_error("Unknown exception caught in FLERInterface.");
    }
}

std::string_view FLERInterface::GetRowNameAt(std::size_t k) const {
    if (k >= static_cast<std::size_t>(shape_[0])) throw MyException("Index out of bounds");
    return rows_name_[k];
}

std::string_view FLERInterface::GetColNameAt(std::size_t k) const {
    if (k >= static_cast<std::size_t>(shape_[1])) throw MyException("Index out of bounds");
    return cols_name_[k];
}

