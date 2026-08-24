
// MSVC: explicit standard includes (not pulled in transitively).
#include <vector>
#include <utility>
#include <cstdint>
/**
 * @file loss_function_mrp_l1.cpp
 * @brief Implementazione delle logiche di calcolo delle loss lineari e V2.
 */

#include "loss_function_mrp_l1.h"
#include <cmath>

double LBMRP::operator()(const Tensor<double, 2>& x) {
    ++calls_;
    double risultato = 0.0;
    
    for (const auto& mc : *mrps_ref_) {
        double denominatore = 0.0;
        double numeratore = 0.0;
        
        for (const auto& [riga, peso_r] : fast_weights_) {
            for (const auto& [colonna, peso_c] : fast_weights_) {
                double p_prod = peso_r * peso_c;
                numeratore += p_prod * std::fabs(x(riga, colonna) - mc(riga, colonna));
                denominatore += p_prod * mc(riga, colonna);
            }
        }
        if (denominatore > 1e-15) risultato += numeratore / denominatore;
    }
    return risultato;
}

double L1LossFunctionMRP::operator()(const Tensor<double, 2>& x) {
    ++calls_;
    double risultato = 0.0;
    const double n_rows = static_cast<double>(x.Extent(0));
    const double denominatore = (n_rows * (n_rows + 1.0)) / 2.0;
    
    for (const auto& mc : *mrps_ref_) {
        double numeratore = 0.0;
        
        for (const auto& [riga, peso_r] : fast_weights_) {
            for (const auto& [colonna, peso_c] : fast_weights_) {
                if (riga < colonna) {
                    numeratore += peso_r * std::fabs(x(riga, colonna) - mc(riga, colonna)) * peso_c;
                }
            }
        }
        risultato += (2.0 * numeratore) / denominatore;
    }
    
    return mrps_ref_->empty() ? 0.0 : risultato / static_cast<double>(mrps_ref_->size());
}

// ---------------- LBMRP2 (V2 Methods) ----------------

void LBMRP2::Resolve(std::vector<ResolvedIdx>& out,
                     const std::vector<std::vector<std::uint_fast64_t>>& funcs) const {
    out.clear();
    out.reserve(fast_weights_.size());
    for (const auto& [idx, w] : fast_weights_) {
        std::uint64_t curr = idx;
        for (const auto& f : funcs) {
            curr = f[curr];
        }
        out.push_back({idx, curr, w});
    }
}

double LBMRP2::operator()(const Tensor<double, 2>& x) {
    ++calls_;
    // LBMRP standard fallback
    double risultato = 0.0;
    for (const auto& mc : *mrps_ref_) {
        double denominatore = 0.0;
        double numeratore = 0.0;
        for (const auto& [riga, peso_r] : fast_weights_) {
            for (const auto& [colonna, peso_c] : fast_weights_) {
                double p_prod = peso_r * peso_c;
                numeratore += p_prod * std::fabs(x(riga, colonna) - mc(riga, colonna));
                denominatore += p_prod * mc(riga, colonna);
            }
        }
        if (denominatore > 1e-15) risultato += numeratore / denominatore;
    }
    return risultato;
}

double LBMRP2::operator()(const Tensor<double, 2>& x,
                          const std::vector<std::vector<std::uint_fast64_t>>& rows_functions,
                          const std::vector<std::vector<std::uint_fast64_t>>& cols_functions) {
    ++calls_;
    std::vector<ResolvedIdx> r_idx, c_idx;
    Resolve(r_idx, rows_functions);
    Resolve(c_idx, cols_functions);
    
    double risultato = 0.0;
    for (const auto& mc : *mrps_ref_) {
        double denominatore = 0.0;
        double numeratore = 0.0;
        
        for (const auto& r : r_idx) {
            for (const auto& c : c_idx) {
                double p_prod = r.w * c.w;
                numeratore += p_prod * std::fabs(x(r.res, c.res) - mc(r.orig, c.orig));
                denominatore += p_prod * mc(r.orig, c.orig);
            }
        }
        //std::cout << "\n" << numeratore << " / " << denominatore << " --> ";

        if (denominatore > 1e-15) risultato += numeratore / denominatore;
    }
    return risultato;
}

void LBMRP2::operator()(const Tensor<double, 2>& x,
                        const std::vector<std::vector<std::uint_fast64_t>>& rows_functions,
                        const std::vector<std::vector<std::uint_fast64_t>>& cols_functions,
                        std::vector<std::pair<std::uint_fast64_t, double>>& result) {
    ++calls_;
    result.clear();
    std::vector<ResolvedIdx> r_idx, c_idx;
    Resolve(r_idx, rows_functions);
    Resolve(c_idx, cols_functions);
    
    for (const auto& mc : *mrps_ref_) {
        for (const auto& c : c_idx) {
            double denominatore = 0.0;
            double numeratore = 0.0;
            
            for (const auto& r : r_idx) {
                numeratore += r.w * std::fabs(x(r.res, c.res) - mc(r.orig, c.orig));
                denominatore += r.w * mc(r.orig, c.orig);
            }
            result.emplace_back(c.orig, (denominatore > 1e-15 ? numeratore / denominatore : 0.0));
        }
    }
}
