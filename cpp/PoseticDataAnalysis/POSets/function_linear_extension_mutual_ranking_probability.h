#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <string_view>
#include <memory>
#include <cstddef>

/**
 * @file function_linear_extension_mutual_ranking_probability.h
 * @brief Mutual Ranking Probability evaluator.
 */

#include "function_linear_extension.h"

class POSet;

/**
 * @class  FLEMutualRankingProbability
 * @brief  Estimates the Mutual Ranking Probability (MRP) matrix of a POSet.
 */
class FLEMutualRankingProbability final : public FunctionLinearExtension {
private:
    const POSet* poset_;
    
public:
    /**
     * @brief Constructs the evaluator for a specific POSet.
     * @param poset Pointer to the POSet instance.
     */
    explicit FLEMutualRankingProbability(const POSet* poset);

    [[nodiscard]] std::unique_ptr<FunctionLinearExtension> Clone() const override {
        return std::make_unique<FLEMutualRankingProbability>(*this);
    }

    void operator()(const LinearExtension& x) noexcept override;
    [[nodiscard]] std::string_view GetRowNameAt(std::size_t k) const override;
    [[nodiscard]] std::string_view GetColNameAt(std::size_t k) const override;
};
