/**
 * @file tensor.h
 * @brief Optimized generic N-dimensional array (Tensor) for HPC environments.
 *
 * @par Requirements
 * C++20 or later (`-std=c++20`).
 *
 * @note Naming follows the Google C++ Style Guide:
 * - Types (classes, structs, concepts): PascalCase
 * - Methods: PascalCase (except STL-compatible ones like at(), size(), data())
 * - Data members: snake_case_
 * - Local variables / parameters: snake_case
 * - Constants / constexpr variables: kCamelCase
 */

#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <utility>

#include "core_tags.h"
#include "my_exception.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <format>
#include <memory>
#include <format>
#include <iterator>
#include <stdexcept>
#include <string>

// ---------------------------------------------------------------------------
// Concepts
// ---------------------------------------------------------------------------

/**
 * @concept TensorElement
 * @brief Constrains the element type usable within the Tensor class.
 * @details The type must be default-initializable and copyable to be stored
 * in the contiguous memory buffer.
 */
template <typename T>
concept TensorElement = std::default_initializable<T> && std::copyable<T>;


// ---------------------------------------------------------------------------
// Main Class
// ---------------------------------------------------------------------------

/**
 * @class Tensor
 * @brief Generic N-dimensional tensor container optimized for HPC.
 * * @details This class provides a multidimensional view over a flat, row-major,
 * contiguous memory block. Dimension bounds and strides are stored on
 * the stack via std::array to avoid heap allocation overhead and to
 * allow the compiler to unroll index computations.
 * * @tparam T The type of elements stored in the tensor.
 * @tparam Rank The number of dimensions (must be strictly greater than 0).
 */
template <TensorElement T, std::size_t Rank>
class Tensor final {
private:
    std::array<std::uint64_t, Rank> shape_{};   ///< Extents of each dimension.
    std::array<std::uint64_t, Rank> strides_{}; ///< Memory steps for each dimension.
    std::uint64_t                   size_{0};   ///< Total number of elements.
    std::unique_ptr<T[]>            data_;      ///< Flat row-major storage buffer.
    
    /**
     * @brief Computes the total size and the strides for flat index calculation.
     * @details Strides are calculated in a row-major fashion (C-style contiguous).
     */
    void ComputeStridesAndSize() noexcept {
        static_assert(Rank > 0, "Tensor: Rank must be greater than 0.");
        
        size_ = 1;
        for (std::size_t i = 0; i < Rank; ++i) {
            size_ *= shape_[i];
        }
        
        if (size_ > 0) {
            strides_[Rank - 1] = 1;
            for (std::size_t i = Rank - 1; i > 0; --i) {
                strides_[i - 1] = strides_[i] * shape_[i];
            }
        }
    }
    
    /**
     * @brief Validates the provided multidimensional indices against the shape.
     * @tparam Idx Variadic template for index types.
     * @param indices The indices to check.
     * @throws MyException if any index is out of bounds.
     */
    template <std::integral... Idx>
    void CheckBounds(Idx... indices) const {
        static_assert(sizeof...(indices) == Rank, "Tensor: Incorrect number of indices provided.");
        const std::array<std::uint64_t, Rank> idx{ static_cast<std::uint64_t>(indices)... };
        for (std::size_t i = 0; i < Rank; ++i) {
            if (idx[i] >= shape_[i]) {
                throw MyException(std::format("Tensor: Index out of bounds at dimension {} ({} >= {}).",
                                              i, idx[i], shape_[i]));
            }
        }
    }
    
public:
    // -------------------------------------------------------------------------
    // Constructors & Destructor
    // -------------------------------------------------------------------------
    
    /**
     * @brief Default constructor. Creates an empty tensor.
     */
    Tensor() = default;
    
    /**
     * @brief HPC constructor: allocates memory WITHOUT initializing it.
     * @param shape An array defining the size of each dimension.
     * @param tag Tag to explicitly request uninitialized memory.
     */
    Tensor(const std::array<std::uint64_t, Rank>& shape, Uninitialized tag)
    : shape_(shape) {
        ComputeStridesAndSize();
        data_ = std::make_unique_for_overwrite<T[]>(static_cast<std::size_t>(size_));
    }
    
    /**
     * @brief Standard constructor: allocates and value-initializes memory.
     * @param shape An array defining the size of each dimension.
     * @param init The default value used to fill the tensor (defaults to T{}).
     */
    explicit Tensor(const std::array<std::uint64_t, Rank>& shape, const T& init = T{})
    : shape_(shape) {
        ComputeStridesAndSize();
        data_ = std::make_unique_for_overwrite<T[]>(static_cast<std::size_t>(size_));
        if (size_ > 0) {
            std::fill_n(data_.get(), static_cast<std::size_t>(size_), init);
        }
    }
    
    /**
     * @brief Copy constructor. Performs a deep copy of the underlying buffer.
     * @param other The tensor to copy from.
     */
    Tensor(const Tensor& other)
    : shape_(other.shape_),
    strides_(other.strides_),
    size_(other.size_),
    data_(std::make_unique_for_overwrite<T[]>(static_cast<std::size_t>(other.size_))) {
        if (size_ > 0 && other.data_) {
            std::copy_n(other.data_.get(), static_cast<std::size_t>(size_), data_.get());
        }
    }
    
    /**
     * @brief Move constructor. Takes ownership of the resources from another tensor.
     * @param other The tensor to move from.
     */
    Tensor(Tensor&& other) noexcept = default;
    
    /**
     * @brief Default destructor.
     */
    ~Tensor() = default;
    
    // -------------------------------------------------------------------------
    // Assignment Operators
    // -------------------------------------------------------------------------
    
    /**
     * @brief Copy assignment operator.
     * @param other The tensor to copy from.
     * @return Reference to the current tensor.
     * @throws MyException If the shapes of the two tensors do not match.
     */
    Tensor& operator=(const Tensor& other) {
        if (this == &other) return *this;
        if (shape_ != other.shape_) {
            throw MyException("Tensor::operator=: Shape mismatch during assignment.");
        }
        if (size_ > 0 && data_ && other.data_) {
            std::copy_n(other.data_.get(), static_cast<std::size_t>(size_), data_.get());
        }
        return *this;
    }
    
    /**
     * @brief Move assignment operator.
     * @param other The tensor to move from.
     * @return Reference to the current tensor.
     * @throws MyException If the shapes of the two tensors do not match (unless the current one is empty).
     */
    Tensor& operator=(Tensor&& other) noexcept {
        if (this == &other) return *this;
        
        shape_   = std::move(other.shape_);
        strides_ = std::move(other.strides_);
        size_    = other.size_;
        data_    = std::move(other.data_);
        
        other.size_ = 0;
        other.shape_.fill(0);
        other.strides_.fill(0);

        return *this;
    }
    
    // -------------------------------------------------------------------------
    // Element Access
    // -------------------------------------------------------------------------
    
    /**
     * @brief Unchecked element access using N-dimensional indices.
     * @note Loop over dimensions is unrolled at compile-time by the optimizer.
     * @tparam Idx Variadic template for index types.
     * @param indices A sequence of indices corresponding to each dimension.
     * @return Reference to the requested element.
     */
    template <std::integral... Idx>
    [[nodiscard]] T& operator()(Idx... indices) noexcept {
        static_assert(sizeof...(indices) == Rank, "Tensor: Incorrect number of indices.");
        const std::array<std::uint64_t, Rank> idx{ static_cast<std::uint64_t>(indices)... };
        std::uint64_t flat_idx = 0;
        
        for (std::size_t i = 0; i < Rank; ++i) {
            assert(idx[i] < shape_[i] && "Tensor: Index out of bounds in operator().");
            flat_idx += idx[i] * strides_[i];
        }
        return data_[flat_idx];
    }
    
    /**
     * @brief Unchecked constant element access using N-dimensional indices.
     * @tparam Idx Variadic template for index types.
     * @param indices A sequence of indices corresponding to each dimension.
     * @return Constant reference to the requested element.
     */
    template <std::integral... Idx>
    [[nodiscard]] const T& operator()(Idx... indices) const noexcept {
        static_assert(sizeof...(indices) == Rank, "Tensor: Incorrect number of indices.");
        const std::array<std::uint64_t, Rank> idx{ static_cast<std::uint64_t>(indices)... };
        std::uint64_t flat_idx = 0;
        
        for (std::size_t i = 0; i < Rank; ++i) {
            assert(idx[i] < shape_[i] && "Tensor: Index out of bounds in operator().");
            flat_idx += idx[i] * strides_[i];
        }
        return data_[flat_idx];
    }
    
    /**
     * @brief Checked element access with boundary validation.
     * @tparam Idx Variadic template for index types.
     * @param indices A sequence of indices corresponding to each dimension.
     * @return Reference to the requested element.
     * @throws MyException if any index is out of bounds.
     */
    template <std::integral... Idx>
    [[nodiscard]] T& at(Idx... indices) {
        CheckBounds(indices...);
        return operator()(indices...);
    }
    
    /**
     * @brief Checked constant element access with boundary validation.
     * @tparam Idx Variadic template for index types.
     * @param indices A sequence of indices corresponding to each dimension.
     * @return Constant reference to the requested element.
     * @throws MyException if any index is out of bounds.
     */
    template <std::integral... Idx>
    [[nodiscard]] const T& at(Idx... indices) const {
        CheckBounds(indices...);
        return operator()(indices...);
    }
    
    // -------------------------------------------------------------------------
    // Observers
    // -------------------------------------------------------------------------
    
    /**
     * @brief Retrieves the size of a specific dimension.
     * @param dim The dimension index (0-based).
     * @return The extent (number of elements) along the specified dimension.
     */
    [[nodiscard]] std::uint64_t Extent(std::size_t dim) const {
        assert(dim < Rank && "Tensor::Extent: Dimension index out of bounds.");
        return shape_[dim];
    }
    
    /**
     * @brief Retrieves the full shape of the tensor.
     * @return A constant reference to the array containing the extents.
     */
    [[nodiscard]] const std::array<std::uint64_t, Rank>& Shape() const noexcept {
        return shape_;
    }
    
    /**
     * @brief Retrieves the total number of elements in the tensor.
     * @return The overall size.
     */
    [[nodiscard]] std::uint64_t size() const noexcept {
        return size_;
    }
    
    /**
     * @brief Provides raw access to the underlying memory buffer.
     * @return Pointer to the beginning of the flat array.
     */
    [[nodiscard]] T* data() noexcept {
        return data_.get();
    }
    
    /**
     * @brief Provides raw access to the underlying constant memory buffer.
     * @return Constant pointer to the beginning of the flat array.
     */
    [[nodiscard]] const T* data() const noexcept {
        return data_.get();
    }
    
    // -------------------------------------------------------------------------
    // Operators
    // -------------------------------------------------------------------------
    
    /**
     * @brief Equality operator.
     * @param other The tensor to compare against.
     * @return True if shapes and all elements match, false otherwise.
     */
    [[nodiscard]] bool operator==(const Tensor& other) const noexcept {
        if (shape_ != other.shape_) return false;
        if (size_ == 0) return true;
        if (data_ && other.data_) {
            return std::equal(data_.get(), data_.get() + static_cast<std::size_t>(size_), other.data_.get());
        }
        return false;
    }
    
    // -------------------------------------------------------------------------
    // Utilities
    // -------------------------------------------------------------------------
    
    /**
     * @brief Serialises the tensor content into a string format.
     * @details If the Rank is 2, it formats the output as a 2D matrix.
     * Otherwise, it defaults to a flat representation.
     * @param delimiter The character used to separate values.
     * @return A string representation including shape and data.
     */
    [[nodiscard]] std::string ToString(char delimiter = ';') const {
        std::string result;
        result.reserve(50 + static_cast<std::size_t>(size_) * 15);
        
        // Intestazione Shape
        std::format_to(std::back_inserter(result), "Shape: [");
        for (std::size_t i = 0; i < Rank; ++i) {
            std::format_to(std::back_inserter(result), "{}{}", shape_[i], (i < Rank - 1 ? " x " : ""));
        }
        std::format_to(std::back_inserter(result), "]\nData:\n");
        
        if constexpr (Rank == 2) {
            const uint64_t rows = shape_[0];
            const uint64_t cols = shape_[1];
            
            for (uint64_t i = 0; i < rows; ++i) {
                for (uint64_t j = 0; j < cols; ++j) {
                    std::format_to(std::back_inserter(result), "{}", (*this)(i, j));
                    if (j < cols - 1) {
                        std::format_to(std::back_inserter(result), "{} ", delimiter);
                    }
                }
                result += '\n'; // Nuova riga alla fine di ogni riga della matrice
            }
        } else {
            // Fallback per altri Rank
            for (std::size_t i = 0; i < static_cast<std::size_t>(size_); ++i) {
                if (i > 0) {
                    std::format_to(std::back_inserter(result), "{} ", delimiter);
                }
                std::format_to(std::back_inserter(result), "{}", data_[i]);
            }
            std::format_to(std::back_inserter(result), "\n");
        }
        
        return result;
    }
    
    // -------------------------------------------------------------------------
    // Iterators (Flat memory access)
    // -------------------------------------------------------------------------
    
    /**
     * @brief Returns an iterator to the beginning of the flat memory buffer.
     */
    [[nodiscard]] T* begin() noexcept { return data_.get(); }
    
    /**
     * @brief Returns an iterator to the end of the flat memory buffer.
     */
    [[nodiscard]] T* end() noexcept { return data_.get() + size_; }
    
    /**
     * @brief Returns a constant iterator to the beginning of the flat memory buffer.
     */
    [[nodiscard]] const T* begin() const noexcept { return data_.get(); }
    
    /**
     * @brief Returns a constant iterator to the end of the flat memory buffer.
     */
    [[nodiscard]] const T* end() const noexcept { return data_.get() + size_; }
};
