
// MSVC: explicit standard includes (not pulled in transitively).
#include <string_view>
#include <cstdint>
#include <cstddef>
#include "function_linear_extension_separation_asymmetric_lower.h"
#include "linear_extension.h"
#include "my_exception.h"
#include "poset.h"

FLESeparationAsymmetricLower::FLESeparationAsymmetricLower(const POSet* poset)
: FunctionLinearExtension(), poset_(poset)
{
    const std::size_t n = static_cast<std::size_t>(poset_->size());
    idx0_.reserve(n * n);
    idx1_.reserve(n * n);
    values_.assign(n * n, 0.0);
    
    // Generiamo le coppie di indici per ogni elemento
    for (std::size_t first = 0; first < n; ++first) {
        for (std::size_t second = 0; second < n; ++second) {
            idx0_.push_back(static_cast<std::uint32_t>(first));
            idx1_.push_back(static_cast<std::uint32_t>(second));
        }
    }
    shape_ = {static_cast<std::uint32_t>(n), static_cast<std::uint32_t>(n)};
}

void FLESeparationAsymmetricLower::operator()(const LinearExtension& x) noexcept {
    ++calls_;
    const std::size_t n = values_.size();
    
    for (std::size_t k = 0; k < n; ++k) {
        const std::uint64_t first_pos = x.GetPos(idx0_[k]);
        const std::uint64_t second_pos = x.GetPos(idx1_[k]);
        
        // Calcoliamo la distanza solo se il secondo elemento segue il primo
        values_[k] = (second_pos > first_pos) ? static_cast<double>(second_pos - first_pos) : 0.0;
    }
}

std::string_view FLESeparationAsymmetricLower::GetRowNameAt(std::size_t k) const {
    if (k >= static_cast<std::size_t>(shape_[0])) throw MyException("Index out of bounds");
    return poset_->GetElementName(k);
}

std::string_view FLESeparationAsymmetricLower::GetColNameAt(std::size_t k) const {
    if (k >= static_cast<std::size_t>(shape_[1])) throw MyException("Index out of bounds");
    return poset_->GetElementName(k);
}
