#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <vector>
#include <memory>
#include <cstdint>
#include <cstddef>
#include <limits>

/**
 * @file function_linear_extension_mann_whitney_inferential_dominance.h
 * @brief Inferential Mann-Whitney dominance evaluator.
 */

#include "tensor.h"
#include "function_linear_extension.h"
#include <span>
#include <string>
#include <string_view>

class POSet;

/**
 * @class  FLEMannWhitneyInferentialDominance
 * @brief  Computes inferential statistics based on the Mann-Whitney U distribution.
 */
class FLEMannWhitneyInferentialDominance final : public FunctionLinearExtension {
private:
    std::size_t poset_size_;
    const Tensor<double, 2>& frequency_matrix_;
    const Tensor<std::string_view, 1>& result_labels_;
    const Tensor<double, 1>& subpopulation_count_;
    std::span<const std::uint32_t> observed_value_;
    std::size_t k_cols_;
    std::size_t n_obs_;
    std::size_t total_bins_;
    
    std::vector<std::uint32_t> obs_row_map_;
    std::vector<std::uint32_t> ordered_indices_;
    Tensor<double, 2> cumulative_matrix_;
    std::vector<double> running_sums_;
    std::vector<double> null_variance_;
    std::vector<double> adjusted_col_g1_;

    inline double standard_normal_cdf(double x) const noexcept;
    
    
    mutable std::string row_name_scratchpad_;
    mutable std::string col_name_scratchpad_;
    
    const std::uint32_t UNOBSERVED_FLAG = std::numeric_limits<std::uint32_t>::max();

public:
    /**
     * @brief Constructs the inferential Mann-Whitney dominance evaluator.
     * @param poset_size Size of the POSet.
     * @param frequency_matrix Flattened frequency matrix across subpopulations.
     * @param result_labels Labels for the columns/groups.
     * @param subpopulation_count Array representing the sample size per subpopulation.
     * @param observed_value The observed element map.
     * @param total_bins Number of bins for histogram allocation of p-values.
     */
    explicit FLEMannWhitneyInferentialDominance(std::size_t poset_size,
                                                const Tensor<double, 2>& frequency_matrix,
                                                const Tensor<std::string_view, 1>& result_labels,
                                                const Tensor<double, 1>& subpopulation_count,
                                                std::span<const std::uint32_t> observed_value,
                                                std::size_t total_bins);

    // NB: il clone condivide (per riferimento) frequency_matrix_, result_labels_ e
    // subpopulation_count_, dati di input a sola lettura: sicuro anche tra thread.
    [[nodiscard]] std::unique_ptr<FunctionLinearExtension> Clone() const override {
        return std::make_unique<FLEMannWhitneyInferentialDominance>(*this);
    }

    void operator()(const LinearExtension& x) noexcept override;
    [[nodiscard]] std::string_view GetRowNameAt(std::size_t k) const override;
    [[nodiscard]] std::string_view GetColNameAt(std::size_t k) const override;
};
