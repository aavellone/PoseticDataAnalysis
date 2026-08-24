#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <string_view>
#include <memory>
#include <cstddef>

/**
 * @file function_linear_extension_average_height.h
 * @brief Average height evaluator.
 */

#include "function_linear_extension.h"

class POSet;

/**
 * @class  FLEAverageHeight
 * @brief  Accumulates the 1-based rank (height) of each POSet element.
 */
class FLEAverageHeight final : public FunctionLinearExtension {
private:
    const POSet* poset_;
    
public:
    /**
     * @brief Constructs the average height evaluator.
     * @param poset Pointer to the POSet instance.
     */
    explicit FLEAverageHeight(const POSet* poset);

    [[nodiscard]] std::unique_ptr<FunctionLinearExtension> Clone() const override {
        return std::make_unique<FLEAverageHeight>(*this);
    }

    void operator()(const LinearExtension& x) noexcept override;
    [[nodiscard]] std::string_view GetRowNameAt(std::size_t k) const override;
    [[nodiscard]] std::string_view GetColNameAt(std::size_t k) const override;
};
