#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <string_view>
#include <memory>
#include <cstddef>

/**
 * @file function_linear_extension_separation_asymmetric_upper.h
 * @brief Upper asymmetric separation evaluator.
 */

#include "function_linear_extension.h"

class POSet;

/**
 * @class  FLESeparationAsymmetricUpper
 * @brief  Computes the upper asymmetric rank separation between element pairs.
 */
class FLESeparationAsymmetricUpper final : public FunctionLinearExtension {
private:
    const POSet* poset_;
    
public:
    /**
     * @brief Constructs the upper asymmetric separation evaluator.
     * @param poset Pointer to the POSet instance.
     */
    explicit FLESeparationAsymmetricUpper(const POSet* poset);

    [[nodiscard]] std::unique_ptr<FunctionLinearExtension> Clone() const override {
        return std::make_unique<FLESeparationAsymmetricUpper>(*this);
    }

    void operator()(const LinearExtension& x) noexcept override;
    [[nodiscard]] std::string_view GetRowNameAt(std::size_t k) const override;
    [[nodiscard]] std::string_view GetColNameAt(std::size_t k) const override;
};
