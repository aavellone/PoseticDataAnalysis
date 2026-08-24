#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <string_view>
#include <memory>
#include <cstddef>

/**
 * @file function_linear_extension_separation_asymmetric_lower.h
 * @brief Lower asymmetric separation evaluator.
 */

#include "function_linear_extension.h"

class POSet;

/**
 * @class  FLESeparationAsymmetricLower
 * @brief  Computes the lower asymmetric rank separation between element pairs.
 */
class FLESeparationAsymmetricLower final : public FunctionLinearExtension {
private:
    const POSet* poset_;
    
public:
    /**
     * @brief Constructs the lower asymmetric separation evaluator.
     * @param poset Pointer to the POSet instance.
     */
    explicit FLESeparationAsymmetricLower(const POSet* poset);

    [[nodiscard]] std::unique_ptr<FunctionLinearExtension> Clone() const override {
        return std::make_unique<FLESeparationAsymmetricLower>(*this);
    }

    void operator()(const LinearExtension& x) noexcept override;
    [[nodiscard]] std::string_view GetRowNameAt(std::size_t k) const override;
    [[nodiscard]] std::string_view GetColNameAt(std::size_t k) const override;
};
