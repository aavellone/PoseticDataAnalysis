
// MSVC: explicit standard includes (not pulled in transitively).
#include <string_view>
#include <numeric>
#include <cstdint>
#include <cstddef>
#include "function_linear_extension_mann_whitney_inferential_dominance.h"
#include "linear_extension.h"
#include "my_exception.h"
#include "tensor.h"
#include "poset.h"
#include <cmath>
#include <numbers>
#include <algorithm>
#include <format>
#include <limits>

FLEMannWhitneyInferentialDominance::FLEMannWhitneyInferentialDominance(
                                                                       std::size_t poset_size,
                                                                       const Tensor<double, 2>& frequency_matrix,
                                                                       const Tensor<std::string_view, 1>& result_labels,
                                                                       const Tensor<double, 1>& subpopulation_count,
                                                                       std::span<const std::uint32_t> observed_value,
                                                                       std::size_t total_bins)
: FunctionLinearExtension(),
poset_size_(poset_size),
frequency_matrix_(frequency_matrix),
result_labels_(result_labels),
subpopulation_count_(subpopulation_count),
observed_value_(observed_value),
k_cols_(result_labels.size()),
n_obs_(observed_value.size()),
total_bins_(total_bins),
cumulative_matrix_({n_obs_, k_cols_}, 0.0)
{
    const std::size_t n = poset_size_;
    
    // Controlliamo che la matrice abbia le dimensioni attese
    if (frequency_matrix_.Shape()[0] != static_cast<std::uint64_t>(n_obs_) ||
        frequency_matrix_.Shape()[1] != static_cast<std::uint64_t>(k_cols_)) {
        throw MyException("Frequency matrix shape mismatch");
    }
    if (subpopulation_count_.size() != k_cols_) {
        throw MyException("Subpopulation count size mismatch");
    }
    
    // Inizializzazione della mappa osservazioni inversa
    obs_row_map_.assign(n, UNOBSERVED_FLAG);
    for (std::size_t i = 0; i < n_obs_; ++i) {
        obs_row_map_[static_cast<std::size_t>(observed_value_[i])] = static_cast<std::uint32_t>(i);
    }
    
    const std::size_t rows = k_cols_ * k_cols_;
    const bool has_bins = (total_bins_ > 0);
    const std::size_t res_size = has_bins ? rows * 2 : rows;
    
    idx0_.resize(res_size);
    idx1_.resize(res_size);
    values_.assign(res_size, 0.0);
    adjusted_col_g1_.assign(n_obs_, 0.0);

    if (has_bins) {
        for (std::size_t i = 0; i < rows; ++i) {
            idx0_[i * 2] = static_cast<std::uint32_t>(i);
            idx0_[i * 2 + 1] = static_cast<std::uint32_t>(i);
        }
    } else {
        std::iota(idx0_.begin(), idx0_.begin() + rows, 0ull);
    }
    
    ordered_indices_.assign(n_obs_, 0);
    running_sums_.assign(k_cols_, 0.0);
    null_variance_.assign(rows, 0.0);
    
    Tensor<double, 2> abs_freq({n_obs_, k_cols_}, 0.0);
    
    for (std::size_t obs = 0; obs < n_obs_; ++obs) {
        for (std::size_t gruppo = 0; gruppo < k_cols_; ++gruppo) {
            abs_freq(obs, gruppo) = frequency_matrix_(obs, gruppo) * subpopulation_count_(gruppo);
        }
    }

    for (std::size_t g1 = 0; g1 < k_cols_; ++g1) {
        const double sg1 = subpopulation_count_(g1);
        for (std::size_t g2 = g1; g2 < k_cols_; ++g2) {
            const double sg2 = subpopulation_count_(g2);
            
            const double a1 = 1.0 / (12.0 * sg1 * sg2);
            const double N = sg1 + sg2;
            const double inv_N_poly = 1.0 / (N * (N + 1.0));
            
            double somma = 0.0;
            for (std::uint64_t k = 0; k < n_obs_; ++k) {
                const double tk = abs_freq(k, g1) + abs_freq(k, g2);
                
                somma += tk * (tk * tk - 1.0);
            }
            
            const double val = a1 * (N + 1.0 - (somma * inv_N_poly));
            
            const std::size_t flat_idx = g1 * k_cols_ + g2;
            const std::size_t flat_idx_inv = g2 * k_cols_ + g1;
            
            null_variance_[flat_idx] = val;
            null_variance_[flat_idx_inv] = val;
        }
    }
    
    shape_ = {static_cast<std::uint32_t>(rows), static_cast<std::uint32_t>(total_bins_ + 1)};
}

inline double FLEMannWhitneyInferentialDominance::standard_normal_cdf(double x) const noexcept {
    return 0.5 * std::erfc(-x / std::numbers::sqrt2);
}

void FLEMannWhitneyInferentialDominance::operator()(const LinearExtension& x) noexcept {
    ++calls_;
    const std::size_t n = poset_size_;
    
    std::size_t p_obs = 0;
    for (std::size_t p = 0; p < n; ++p) {
        const std::uint64_t elem = x.GetVal(static_cast<std::uint64_t>(p));
        const std::uint64_t row_idx = obs_row_map_[static_cast<std::size_t>(elem)];
        if (row_idx != UNOBSERVED_FLAG) ordered_indices_[p_obs++] = static_cast<std::uint32_t>(row_idx);
    }
    
    std::fill(running_sums_.begin(), running_sums_.end(), 0.0);
    for (std::size_t p = 0; p < n_obs_; ++p) {
        for (std::size_t j = 0; j < k_cols_; ++j) {
            running_sums_[j] += frequency_matrix_(ordered_indices_[p], j);
            cumulative_matrix_(ordered_indices_[p], j) = running_sums_[j];
        }
    }
    
    // Calcoliamo U e P-value
    const bool has_bins = (total_bins_ > 0);
    
    if (has_bins) {
        std::size_t i = 0;
        for (std::size_t g1 = 0; g1 < k_cols_; ++g1) {
            for (std::size_t row = 0; row < n_obs_; ++row) {
                adjusted_col_g1_[row] = cumulative_matrix_(row, g1) - 0.5 * frequency_matrix_(row, g1);
            }
            for (std::size_t g2 = 0; g2 < k_cols_; ++g2) {
                double u_stat = 0.0;
                for (std::size_t row = 0; row < n_obs_; ++row) {
                    u_stat += frequency_matrix_(row, g2) * adjusted_col_g1_[row];
                }
                
                const std::size_t flat_idx = static_cast<std::size_t>(idx0_[i]);

                // Z-score ed Estrazione p-value
                double z = (u_stat - 0.5) / std::sqrt(null_variance_[flat_idx]);
                double p_val = 1.0 - standard_normal_cdf(z);
                
                std::size_t bin = std::min(static_cast<std::size_t>(p_val * total_bins_), total_bins_ - 1);
                idx1_[i] = static_cast<std::uint32_t>(bin);
                values_[i] = 1.0;
                
                idx1_[i+1] = static_cast<std::uint32_t>(total_bins_);
                values_[i+1] = p_val;
                
                i += 2;
            }
        }
    } else {
        std::size_t i = 0;
        for (std::size_t g1 = 0; g1 < k_cols_; ++g1) {
            for (std::size_t row = 0; row < n_obs_; ++row) {
                adjusted_col_g1_[row] = cumulative_matrix_(row, g1) - 0.5 * frequency_matrix_(row, g1);
            }
            for (std::size_t g2 = 0; g2 < k_cols_; ++g2) {
                double u_stat = 0.0;
                for (std::size_t row = 0; row < n_obs_; ++row) {
                    u_stat += frequency_matrix_(row, g2) * adjusted_col_g1_[row];
                }
                
                // Z-score ed Estrazione p-value
                const std::size_t flat_idx = static_cast<std::size_t>(idx0_[i]);

                double z = (u_stat - 0.5) / std::sqrt(null_variance_[flat_idx]);
                double p_val = 1.0 - standard_normal_cdf(z);
                
                idx1_[i] = 0;
                values_[i] = p_val;
                
                i++;
            }
        }
    }
    
}

std::string_view FLEMannWhitneyInferentialDominance::GetRowNameAt(std::size_t k) const {
    if (k >= static_cast<std::size_t>(shape_[0])) throw MyException("Index out of bounds");
    const std::size_t g1 = k / k_cols_;
    const std::size_t g2 = k % k_cols_;
    
    row_name_scratchpad_.clear();
    std::format_to(std::back_inserter(row_name_scratchpad_), "{} <= {}", result_labels_(g1), result_labels_(g2));
    return row_name_scratchpad_;
    
}

std::string_view FLEMannWhitneyInferentialDominance::GetColNameAt(std::size_t k) const {
    if (k >= static_cast<std::size_t>(shape_[1])) throw MyException("Index out of bounds");
    if (k == total_bins_) return "MWID";
    
    col_name_scratchpad_.clear();
    std::format_to(std::back_inserter(col_name_scratchpad_), "Bin_{}", k);
    return col_name_scratchpad_;
}
