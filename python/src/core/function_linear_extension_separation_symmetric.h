#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <string_view>
#include <memory>
#include <cstddef>

/**
 * @file function_linear_extension_separation_symmetric.h
 * @brief Symmetric separation evaluator.
 */

#include "function_linear_extension.h"

class POSet;

/**
 * @class  FLESeparationSymmetric
 * @brief  Computes the absolute rank separation between element pairs.
 */
class FLESeparationSymmetric final : public FunctionLinearExtension {
private:
    const POSet* poset_;
    
public:
    /**
     * @brief Constructs the symmetric separation evaluator.
     * @param poset Pointer to the POSet instance.
     */
    explicit FLESeparationSymmetric(const POSet* poset);

    [[nodiscard]] std::unique_ptr<FunctionLinearExtension> Clone() const override {
        return std::make_unique<FLESeparationSymmetric>(*this);
    }

    void operator()(const LinearExtension& x) noexcept override;
    [[nodiscard]] std::string_view GetRowNameAt(std::size_t k) const override;
    [[nodiscard]] std::string_view GetColNameAt(std::size_t k) const override;
};
