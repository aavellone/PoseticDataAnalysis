
// MSVC: explicit standard includes (not pulled in transitively).
#include <string_view>
#include <cstdint>
#include <cstddef>
#include "function_linear_extension_mutual_ranking_probability.h"
#include "linear_extension.h"
#include "my_exception.h"
#include "poset.h"

FLEMutualRankingProbability::FLEMutualRankingProbability(const POSet* poset)
: FunctionLinearExtension(), poset_(poset)
{
    // Otteniamo la dimensione del poset usando size_t per la logica locale
    const std::size_t n = static_cast<std::size_t>(poset_->size());
    
    // Pre-allochiamo la memoria per evitare riallocazioni
    idx0_.reserve(n * n);
    idx1_.reserve(n * n);
    values_.assign(n * n, 0.0);
    
    // Costruiamo la griglia degli indici
    for (std::size_t first = 0; first < n; ++first) {
        for (std::size_t second = 0; second < n; ++second) {
            // Salviamo a 64-bit come richiesto per l'HPC / compatibilità SIMD
            idx0_.push_back(static_cast<std::uint32_t>(first));
            idx1_.push_back(static_cast<std::uint32_t>(second));
        }
    }
    shape_ = {static_cast<std::uint32_t>(n), static_cast<std::uint32_t>(n)};
}

void FLEMutualRankingProbability::operator()(const LinearExtension& x) noexcept {
    ++calls_;
    const std::size_t total_size = values_.size();
    
    // Ciclo hot-path: iteriamo strettamente con std::size_t
    for (std::size_t k = 0; k < total_size; ++k) {
        const std::uint64_t first_pos = x.GetPos(idx0_[k]);
        const std::uint64_t second_pos = x.GetPos(idx1_[k]);
        
        // Assegniamo 1.0 se il primo precede il secondo, altrimenti 0.0
        values_[k] = static_cast<double>(first_pos <= second_pos);
    }
}

std::string_view FLEMutualRankingProbability::GetRowNameAt(std::size_t k) const {
    if (k >= static_cast<std::size_t>(shape_[0])) throw MyException("Index out of bounds");
    return poset_->GetElementName(k);
}

std::string_view FLEMutualRankingProbability::GetColNameAt(std::size_t k) const {
    if (k >= static_cast<std::size_t>(shape_[1])) throw MyException("Index out of bounds");
    return poset_->GetElementName(k);
}
