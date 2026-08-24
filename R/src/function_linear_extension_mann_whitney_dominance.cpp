
// MSVC: explicit standard includes (not pulled in transitively).
#include <string_view>
#include <numeric>
#include <cstdint>
#include <cstddef>
#include "function_linear_extension_mann_whitney_dominance.h"
#include "linear_extension.h"
#include "my_exception.h"
#include "tensor.h"
#include "poset.h"
#include <limits>
#include <algorithm>
#include <format>

FLEMannWhitneyDominance::FLEMannWhitneyDominance(
                                                 std::size_t poset_size,
                                                 const Tensor<double, 2>& frequency_matrix,
                                                 const Tensor<std::string_view, 1>& result_labels,
                                                 std::span<const std::uint32_t> observed_value,
                                                 std::size_t total_bins)
: FunctionLinearExtension(),
poset_size_(poset_size),
frequency_matrix_(frequency_matrix),
result_labels_(result_labels),
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
    
    shape_ = {static_cast<std::uint32_t>(rows), static_cast<std::uint32_t>(total_bins_ + 1)};
}

void FLEMannWhitneyDominance::operator()(const LinearExtension& x) noexcept {
    ++calls_;
    const std::size_t n = poset_size_;

    // Estraiamo e filtriamo gli indici ordinati
    std::size_t p_obs = 0;
    for (std::size_t p = 0; p < n; ++p) {
        const std::uint64_t elem = x.GetVal(static_cast<std::uint64_t>(p));
        const std::uint64_t row_idx = obs_row_map_[static_cast<std::size_t>(elem)];
        if (row_idx != UNOBSERVED_FLAG) {
            ordered_indices_[p_obs++] = static_cast<std::uint32_t>(row_idx);
        }
    }
    
    // Matrice cumulativa delle frequenze
    std::fill(running_sums_.begin(), running_sums_.end(), 0.0);
    for (std::size_t p = 0; p < n_obs_; ++p) {
        for (std::size_t j = 0; j < k_cols_; ++j) {
            running_sums_[j] += frequency_matrix_(ordered_indices_[p], j);
            cumulative_matrix_(ordered_indices_[p], j) = running_sums_[j];
        }
    }
    
    const bool has_bins = (total_bins_ > 0);
    
    // Calcoliamo le stime U di Mann-Whitney
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
                std::size_t bin = std::min(static_cast<std::size_t>(u_stat * total_bins_), total_bins_ - 1);
                
                idx1_[i] = static_cast<std::uint32_t>(bin);
                values_[i] = 1.0;
                
                idx1_[i + 1] = static_cast<std::uint32_t>(total_bins_);
                values_[i + 1] = u_stat;
                
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
                
                idx1_[i] = 0;
                values_[i] = u_stat;
                i++;
            }
        }
    }
}

std::string_view FLEMannWhitneyDominance::GetRowNameAt(std::size_t k) const {
    if (k >= static_cast<std::size_t>(shape_[0])) throw MyException("Index out of bounds");
    const std::size_t g1 = k / k_cols_;
    const std::size_t g2 = k % k_cols_;
    
    row_name_scratchpad_.clear();
    std::format_to(std::back_inserter(row_name_scratchpad_), "{} <= {}", result_labels_(g1), result_labels_(g2));
    return row_name_scratchpad_;
}

std::string_view FLEMannWhitneyDominance::GetColNameAt(std::size_t k) const {
    if (k >= static_cast<std::size_t>(shape_[1])) throw MyException("Index out of bounds");
    if (k == total_bins_) return "MWD";
    
    col_name_scratchpad_.clear();
    std::format_to(std::back_inserter(col_name_scratchpad_), "Bin_{}", k);
    return col_name_scratchpad_;
    
}
