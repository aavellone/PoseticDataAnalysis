
// MSVC: explicit standard includes (not pulled in transitively).
#include <vector>
#include <unordered_map>
#include <cstdint>
/**
 * @file loss_function_mrp_js.cpp
 * @brief Implementazione per le funzioni di perdita basate sulla Teoria dell'Informazione.
 */

#include "loss_function_mrp_js.h"
#include "loss_function_mrp_hpc_math.h"
#include <cmath>

JS0LossFunctionMRP::JS0LossFunctionMRP(const std::vector<Tensor<double, 2>>& p,
                                       const std::unordered_map<std::uint64_t, double>& w,
                                       std::uint64_t base)
: LossFunctionMRP(p, w), inv_log_base_(1.0 / std::log(static_cast<double>(base))) {}

double JS0LossFunctionMRP::operator()(const Tensor<double, 2>& x) {
    ++calls_;
    const auto n = x.Extent(0);
    const auto k = mrps_ref_->size();
    if (k == 0 || n < 2) return 0.0;
    
    double somma = 0.0;
    for (const auto& mc : *mrps_ref_) {
        double somma_js = 0.0;
        
        for (std::uint64_t r = 0; r < n; ++r) {
            for (std::uint64_t s = 0; s < r; ++s) {
                somma_js += hpc_math::JSDivergence2(mc(r, s), x(r, s), inv_log_base_);
            }
        }
        somma += somma_js;
    }
    return (2.0 / (static_cast<double>(k * n * (n - 1)))) * somma;
}


JS1LossFunctionMRP::JS1LossFunctionMRP(const std::vector<Tensor<double, 2>>& p,
                                       const std::unordered_map<std::uint64_t, double>& w,
                                       std::uint64_t base)
: LossFunctionMRP(p, w), inv_log_base_(1.0 / std::log(static_cast<double>(base))) {}

double JS1LossFunctionMRP::operator()(const Tensor<double, 2>& x) {
    ++calls_;
    const auto n = x.Extent(0);
    double somma_numeratore = 0.0;
    double somma_denominatore = 0.0;
    
    for (const auto& mc : *mrps_ref_) {
        double somma_js = 0.0;
        double somma_uno_meno_entropia = 0.0;
        
        for (std::uint64_t r = 0; r < n; ++r) {
            for (std::uint64_t s = 0; s < r; ++s) {
                double p_val = mc(r, s);
                double q_val = x(r, s);
                
                double uno_meno_h = 1.0 - hpc_math::Entropy2(p_val, inv_log_base_);
                somma_js += hpc_math::JSDivergence2(p_val, q_val, inv_log_base_) * uno_meno_h;
                somma_uno_meno_entropia += uno_meno_h;
            }
        }
        somma_numeratore += somma_js;
        somma_denominatore += somma_uno_meno_entropia;
    }
    return (somma_denominatore > 1e-15) ? (somma_numeratore / somma_denominatore) : 0.0;
}


JS2LossFunctionMRP::JS2LossFunctionMRP(const std::vector<Tensor<double, 2>>& p,
                                       const std::unordered_map<std::uint64_t, double>& w,
                                       std::uint64_t base)
: LossFunctionMRP(p, w), inv_log_base_(1.0 / std::log(static_cast<double>(base))) {}

double JS2LossFunctionMRP::operator()(const Tensor<double, 2>& x) {
    ++calls_;
    const auto n = x.Extent(0);
    if (n < 2) return 0.0;
    
    double somma_numeratore = 0.0;
    double somma_denominatore = 0.0;
    const double norm_factor = 2.0 / static_cast<double>(n * (n - 1));
    
    for (const auto& mc : *mrps_ref_) {
        double somma_h = 0.0;
        double somma_js = 0.0;
        
        for (std::uint64_t r = 0; r < n; ++r) {
            for (std::uint64_t s = 0; s < r; ++s) {
                double p_val = mc(r, s);
                double q_val = x(r, s);
                
                somma_js += hpc_math::JSDivergence2(p_val, q_val, inv_log_base_);
                somma_h += hpc_math::Entropy2(p_val, inv_log_base_);
            }
        }
        double h_medio = norm_factor * somma_h;
        somma_numeratore += (1.0 - h_medio) * somma_js;
        somma_denominatore += (1.0 - h_medio);
    }
    return (somma_denominatore > 1e-15) ? (somma_numeratore / somma_denominatore) * norm_factor : 0.0;
}
