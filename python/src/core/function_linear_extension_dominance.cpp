
// MSVC: explicit standard includes (not pulled in transitively).
#include <string_view>
#include <cstdint>
#include <cstddef>
#include "function_linear_extension_dominance.h"
#include "linear_extension.h"
#include "my_exception.h"
#include "tensor.h"
#include "poset.h"
#include <limits>
#include <algorithm>

FLEDominance::FLEDominance(std::size_t poset_size,
                           const Tensor<double, 2>& frequency_matrix,
                           const Tensor<std::string_view, 1>& result_labels,
                           std::span<const std::uint32_t> observed_value)
: FunctionLinearExtension(),
poset_size_(poset_size),
frequency_matrix_(frequency_matrix),
result_labels_(result_labels),
observed_value_(observed_value),
k_cols_(result_labels.size()),
n_obs_(observed_value.size()),
cumulative_matrix_({n_obs_, k_cols_}, 0.0)
{
    const std::size_t n = poset_size_;
    
    // Controlliamo che la matrice abbia le dimensioni attese
    if (frequency_matrix_.Shape()[0] != static_cast<std::uint64_t>(n_obs_) ||
        frequency_matrix_.Shape()[1] != static_cast<std::uint64_t>(k_cols_)) {
        throw MyException("Frequency matrix shape mismatch");
    }
    
    // Inizializziamo la mappa delle osservazioni
    obs_row_map_.assign(n, UNOBSERVED_FLAG);
    for (std::size_t i = 0; i < n_obs_; ++i) {
        obs_row_map_[static_cast<std::size_t>(observed_value_[i])] = static_cast<std::uint32_t>(i);
    }
    
    const std::size_t res_size = k_cols_ * k_cols_;
    idx0_.reserve(res_size);
    idx1_.reserve(res_size);
    values_.assign(res_size, 0.0);
    
    for (std::size_t gruppo1 = 0; gruppo1 < k_cols_; ++gruppo1) {
        for (std::size_t gruppo2 = 0; gruppo2 < k_cols_; ++gruppo2) {
            idx0_.push_back(static_cast<std::uint32_t>(gruppo1));
            idx1_.push_back(static_cast<std::uint32_t>(gruppo2));
        }
    }
    
    // Allocazione buffer di lavoro
    ordered_indices_.assign(n_obs_, 0);
    running_sums_.assign(k_cols_, 0.0);
    
    shape_ = {static_cast<std::uint32_t>(k_cols_), static_cast<std::uint32_t>(k_cols_)};
}

void FLEDominance::operator()(const LinearExtension& x) noexcept {
    ++calls_;
    const std::size_t n = poset_size_;
    
    // 1. Troviamo l'ordine degli indici osservati rispetto a questa estensione lineare
    std::size_t p_obs = 0;
    for (std::size_t p = 0; p < n; ++p) {
        const std::uint64_t elem = x.GetVal(static_cast<std::uint64_t>(p));
        const std::uint64_t row_idx = obs_row_map_[static_cast<std::size_t>(elem)];
        
        if (row_idx != UNOBSERVED_FLAG) {
            ordered_indices_[p_obs++] = static_cast<std::uint32_t>(row_idx);
        }
    }
    
    // 2. Costruiamo la matrice cumulativa calcolando le somme correnti per colonna
    std::fill(running_sums_.begin(), running_sums_.end(), 0.0);
    for (std::size_t p = 0; p < n_obs_; ++p) {
        for (std::size_t j = 0; j < k_cols_; ++j) {
            running_sums_[j] += frequency_matrix_(ordered_indices_[p], j);
            cumulative_matrix_(ordered_indices_[p], j) = running_sums_[j];
        }
    }
    
    // 3. Verifichiamo la dominanza stocastica su tutte le coppie di variabili
    const std::size_t total_size = values_.size();
    for (std::size_t idx = 0; idx < total_size; ++idx) {
        const std::size_t gruppo1 = static_cast<std::size_t>(idx0_[idx]);
        const std::size_t gruppo2 = static_cast<std::size_t>(idx1_[idx]);
        
        if (gruppo1 == gruppo2) {
            values_[idx] = 1.0;
            continue;
        }
        
        bool always_greater_or_equal = true;
        for (std::size_t i = 0; i < n_obs_; ++i) {
            if (cumulative_matrix_(i, gruppo1) < cumulative_matrix_(i, gruppo2)) {
                always_greater_or_equal = false;
                break;
            }
        }
        values_[idx] = always_greater_or_equal ? 1.0 : 0.0;
    }
}

std::string_view FLEDominance::GetRowNameAt(std::size_t k) const {
    if (k >= static_cast<std::size_t>(shape_[0])) throw MyException("Index out of bounds");
    return result_labels_(k);
}

std::string_view FLEDominance::GetColNameAt(std::size_t k) const {
    if (k >= static_cast<std::size_t>(shape_[1])) throw MyException("Index out of bounds");
    return result_labels_(k);
}
