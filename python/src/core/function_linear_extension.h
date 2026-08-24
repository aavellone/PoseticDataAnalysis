#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <string_view>
#include <cstddef>

/**
 * @file   function_linear_extension.h
 * @brief  Abstract base class for evaluating statistical functions over the set
 * of linear extensions of a partially ordered set.
 */

#include <cstdint>
#include <memory>
#include <vector>
#include <string>

class LinearExtension;

/**
 * @class  FunctionLinearExtension
 * @brief  Abstract base class for stateful functors evaluated over linear extensions.
 * * This class provides the base structure for computing statistics over linear extensions
 * in HPC environments. It maintains flat, contiguous vectors for indices and values.
 */
class FunctionLinearExtension {
protected:
    std::uint64_t calls_{0};
    std::vector<std::uint32_t> idx0_;
    std::vector<std::uint32_t> idx1_;
    std::vector<double> values_;
    std::vector<std::uint32_t> shape_;
    
public:
    /**
     * @brief Default constructor.
     */
    FunctionLinearExtension() = default;
    
    /**
     * @brief Virtual destructor.
     */
    virtual ~FunctionLinearExtension() = default;
    
    /**
     * @brief  Returns the first index (row/group) for the flattened result at position k.
     * @param  k  The flat index.
     * @return The pre-computed uint64_t index.
     */
    [[nodiscard]] inline std::uint64_t At0(std::size_t k) const noexcept {
        return idx0_[k];
    }
    
    /**
     * @brief  Returns the second index (col/bin) for the flattened result at position k.
     * @param  k  The flat index.
     * @return The pre-computed uint64_t index.
     */
    [[nodiscard]] inline std::uint64_t At1(std::size_t k) const noexcept {
        return idx1_[k];
    }
    
    /**
     * @brief  Returns the accumulated value at position k.
     * @param  k  The flat index.
     * @return The stored double value.
     */
    [[nodiscard]] inline double At2(std::size_t k) const noexcept {
        return values_[k];
    }
    
    /**
     * @brief  Returns the total number of flattened results stored.
     * @return The size of the values array.
     */
    [[nodiscard]] std::size_t ResSize() const noexcept {
        return values_.size();
    }
    
    /**
     * @brief  Retrieves the shape of the resulting matrix/tensor.
     * @return Reference to the shape vector.
     */
    [[nodiscard]] std::vector<std::uint32_t>& Shape() noexcept { return shape_; }
    
    /**
     * @brief  Creates an independent copy of the functor (internal state included).
     * @details Required by the multi-chain parallel evaluation: each worker thread
     * operates on its own clone to avoid data races on the internal buffers
     * (idx0_, idx1_, values_).
     * @return Owning pointer to the new copy.
     */
    [[nodiscard]] virtual std::unique_ptr<FunctionLinearExtension> Clone() const = 0;

    /**
     * @brief  Tells whether the functor may run outside the main thread.
     * @details Implementations that call the R API (not thread-safe) must
     * return false; they are rejected by the parallel evaluation engine.
     */
    [[nodiscard]] virtual bool IsThreadSafe() const noexcept { return true; }

    /**
     * @brief  Evaluates the statistic for a single linear extension.
     * @param  x  The linear extension to process.
     */
    virtual void operator()(const LinearExtension& x) noexcept = 0;
    
    /**
     * @brief  Retrieves the human-readable name for a given row index.
     * @param  k  The row index.
     * @return The name as a string.
     */
    [[nodiscard]] virtual std::string_view GetRowNameAt(std::size_t k) const = 0;
    
    /**
     * @brief  Retrieves the human-readable name for a given column index.
     * @param  k  The column index.
     * @return The name as a string.
     */
    [[nodiscard]] virtual std::string_view GetColNameAt(std::size_t k) const = 0;
};
