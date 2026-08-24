
// MSVC: explicit standard includes (not pulled in transitively).
#include <string_view>
#include <cstdint>
#include <cstddef>
#include "function_linear_extension_separation_symmetric.h"
#include "linear_extension.h"
#include "my_exception.h"
#include "poset.h"

FLESeparationSymmetric::FLESeparationSymmetric(const POSet* poset)
: FunctionLinearExtension(), poset_(poset)
{
    const std::size_t n = static_cast<std::size_t>(poset_->size());
    idx0_.reserve(n * n);
    idx1_.reserve(n * n);
    values_.assign(n * n, 0.0);
    
    for (std::size_t first = 0; first < n; ++first) {
        for (std::size_t second = 0; second < n; ++second) {
            idx0_.push_back(static_cast<std::uint32_t>(first));
            idx1_.push_back(static_cast<std::uint32_t>(second));
        }
    }
    shape_ = {static_cast<std::uint32_t>(n), static_cast<std::uint32_t>(n)};
}

void FLESeparationSymmetric::operator()(const LinearExtension& x) noexcept {
    ++calls_;
    const std::size_t n = values_.size();
    
    for (std::size_t k = 0; k < n; ++k) {
        const std::uint64_t first_pos = x.GetPos(idx0_[k]);
        const std::uint64_t second_pos = x.GetPos(idx1_[k]);
        
        // Calcoliamo la distanza in valore assoluto
        const std::uint64_t diff = (first_pos > second_pos) ? (first_pos - second_pos) : (second_pos - first_pos);
        values_[k] = static_cast<double>(diff);
    }
}

std::string_view FLESeparationSymmetric::GetRowNameAt(std::size_t k) const {
    if (k >= static_cast<std::size_t>(shape_[0])) throw MyException("Index out of bounds");
    return poset_->GetElementName(k);
}

std::string_view FLESeparationSymmetric::GetColNameAt(std::size_t k) const {
    if (k >= static_cast<std::size_t>(shape_[1])) throw MyException("Index out of bounds");
    return poset_->GetElementName(k);
}
