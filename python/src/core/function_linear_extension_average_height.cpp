
// MSVC: explicit standard includes (not pulled in transitively).
#include <string_view>
#include <cstdint>
#include <cstddef>
#include "function_linear_extension_average_height.h"
#include "linear_extension.h"
#include "my_exception.h"
#include "poset.h"

FLEAverageHeight::FLEAverageHeight(const POSet* poset)
: FunctionLinearExtension(), poset_(poset)
{
    const std::size_t n = static_cast<std::size_t>(poset_->size());
    idx0_.reserve(n);
    idx1_.assign(n, 0); // La colonna è sempre 0
    values_.assign(n, 0.0);
    
    for(std::size_t i = 0; i < n; ++i ) {
        idx0_.push_back(static_cast<std::uint32_t>(i));
    }
    shape_ = {static_cast<std::uint32_t>(n), 1};
}

void FLEAverageHeight::operator()(const LinearExtension& x) noexcept {
    ++calls_;
    const std::size_t n = values_.size();
    
    for (std::size_t k = 0; k < n; ++k) {
        values_[k] = static_cast<double>(x.GetPos(idx0_[k]) + 1);
    }
}

std::string_view FLEAverageHeight::GetRowNameAt(std::size_t k) const {
    if (k >= static_cast<std::size_t>(shape_[0])) throw MyException("Index out of bounds");
    return poset_->GetElementName(k);
}

std::string_view FLEAverageHeight::GetColNameAt(std::size_t k) const {
    if (k >= static_cast<std::size_t>(shape_[1])) throw MyException("Index out of bounds");
    return "AverageHeight";
}
