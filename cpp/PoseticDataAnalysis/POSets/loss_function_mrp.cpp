
// MSVC: explicit standard includes (not pulled in transitively).
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>
/**
 * @file loss_function_mrp.cpp
 * @brief Implementazione dei metodi comuni per l'interfaccia LossFunctionMRP.
 */

#include "loss_function_mrp.h"

LossFunctionMRP::LossFunctionMRP(const std::vector<Tensor<double, 2>>& p,
                                 const std::unordered_map<std::uint64_t, double>& w)
: mrps_ref_(&p) {
    // Appiattimento della mappa per azzerare il cache-miss rate durante le iterazioni
    fast_weights_.reserve(w.size());
    for (const auto& kv : w) {
        fast_weights_.push_back(kv);
    }
}

std::string LossFunctionMRP::to_string() const {
    return "Calls: " + std::to_string(calls_);
}
