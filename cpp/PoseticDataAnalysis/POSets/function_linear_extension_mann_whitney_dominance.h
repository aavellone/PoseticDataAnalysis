#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <vector>
#include <memory>
#include <cstdint>
#include <cstddef>
#include <limits>

/**
 * @file function_linear_extension_mann_whitney_dominance.h
 * @brief Mann-Whitney dominance evaluator.
 */

#include "tensor.h"
#include "function_linear_extension.h"
#include <span>
#include <string>
#include <string_view>

class POSet;

/**
 * @class  FLEMannWhitneyDominance
 * @brief  Accumulates a Mann–Whitney U-like dominance statistic over empirical distributions.
 */
class FLEMannWhitneyDominance final : public FunctionLinearExtension {
private:
    std::size_t poset_size_;
    const Tensor<double, 2>& frequency_matrix_;
    const Tensor<std::string_view, 1>& result_labels_;
    std::span<const std::uint32_t> observed_value_;
    std::size_t k_cols_;
    std::size_t n_obs_;
    std::size_t total_bins_;
    
    std::vector<std::uint32_t> obs_row_map_;
    std::vector<std::uint32_t> ordered_indices_;
    Tensor<double, 2> cumulative_matrix_;
    std::vector<double> running_sums_;
    std::vector<double> adjusted_col_g1_;
    
    mutable std::string row_name_scratchpad_;
    
    const std::uint32_t UNOBSERVED_FLAG = std::numeric_limits<std::uint32_t>::max();

    mutable std::string col_name_scratchpad_;
public:
    /**
     * @brief Constructs the Mann-Whitney dominance evaluator.
     * @param poset_size Size of the POSet.
     * @param frequency_matrix Flattened frequency matrix across subpopulations.
     * @param result_labels Labels for the columns/groups.
     * @param observed_value The observed element map.
     * @param total_bins Number of bins for histogram allocation.
     */
    explicit FLEMannWhitneyDominance(std::size_t poset_size,
                                     const Tensor<double, 2>& frequency_matrix,
                                     const Tensor<std::string_view, 1>& result_labels,
                                     std::span<const std::uint32_t> observed_value,
                                     std::size_t total_bins);

    // NB: il clone condivide (per riferimento) frequency_matrix_ e result_labels_,
    // che sono dati di input a sola lettura: sicuro anche tra thread.
    [[nodiscard]] std::unique_ptr<FunctionLinearExtension> Clone() const override {
        return std::make_unique<FLEMannWhitneyDominance>(*this);
    }

    void operator()(const LinearExtension& x) noexcept override;
    [[nodiscard]] std::string_view GetRowNameAt(std::size_t k) const override;
    [[nodiscard]] std::string_view GetColNameAt(std::size_t k) const override;
};
