
// MSVC: explicit standard includes (not pulled in transitively).
#include <string>
#include <string_view>
#include <vector>
#include <array>
#include <memory>
#include <utility>
#include <cstdint>
#include <cstddef>
#include <format>
#include "linear_poset.h"
#include "product_linear_poset.h"

#include <stdexcept>

// ***********************************************
// Build
// ***********************************************

std::unique_ptr<ProductLinearPOSet> ProductLinearPOSet::Build(const POSet& p1, const POSet& p2) {
    // 1. Validazione RTTI rigorosa (accesso istantaneo tramite reference)
    if (!dynamic_cast<const LinearPOSet*>(&p1) || dynamic_cast<const ProductLinearPOSet*>(&p1)) {
        throw std::runtime_error("ProductLinearPOSet error: wrong poset (p1 is not strictly a LinearPOSet)");
    }
    if (!dynamic_cast<const LinearPOSet*>(&p2) || dynamic_cast<const ProductLinearPOSet*>(&p2)) {
        throw std::runtime_error("ProductLinearPOSet error: wrong poset (p2 is not strictly a LinearPOSet)");
    }
    
    std::unique_ptr<ProductLinearPOSet> result(new ProductLinearPOSet());
    
    const std::uint64_t p1_size = p1.size();
    const std::uint64_t p2_size = p2.size();
    const std::uint64_t elements_size = p1_size * p2_size;
    
    // 2. Architettura RTTI per la storia (usando observer pointers: const POSet*)
    if (const auto* pp1 = dynamic_cast<const FromPOSets*>(&p1)) {
        result->posets_.insert(result->posets_.end(), pp1->posets_.begin(), pp1->posets_.end());
    } else {
        result->posets_.push_back(&p1);
    }
    
    if (const auto* pp2 = dynamic_cast<const FromPOSets*>(&p2)) {
        result->posets_.insert(result->posets_.end(), pp2->posets_.begin(), pp2->posets_.end());
    } else {
        result->posets_.push_back(&p2);
    }
    
    // 3. Pre-allocazione memoria per prestazioni estreme
    std::vector<std::string> elements_str;
    elements_str.reserve(elements_size);
    result->elements_to_original_.reserve(elements_size);
    
    struct PairIndices { std::uint64_t a1, a2; };
    std::vector<PairIndices> indices;
    indices.reserve(elements_size);
    
    const auto* pp1_composite = dynamic_cast<const FromPOSets*>(&p1);
    const auto* pp2_composite = dynamic_cast<const FromPOSets*>(&p2);
    
    for (std::uint64_t i = 0; i < p1_size; ++i) {
        std::string_view name1 = p1.GetElementName(i);
        
        for (std::uint64_t j = 0; j < p2_size; ++j) {
            std::string_view name2 = p2.GetElementName(j);
            elements_str.emplace_back(std::format("{}{}{}", name1, FromPOSets::kSep, name2));
            indices.push_back({i, j});
            
            std::vector<std::uint64_t> orig_eids;
            orig_eids.reserve(result->posets_.size());
            
            if (pp1_composite) {
                const auto& p1_hist = pp1_composite->elements_to_original_[i];
                orig_eids.insert(orig_eids.end(), p1_hist.begin(), p1_hist.end());
            } else {
                orig_eids.push_back(i);
            }
            
            if (pp2_composite) {
                const auto& p2_hist = pp2_composite->elements_to_original_[j];
                orig_eids.insert(orig_eids.end(), p2_hist.begin(), p2_hist.end());
            } else {
                orig_eids.push_back(j);
            }
            
            result->elements_to_original_.push_back(std::move(orig_eids));
        }
    }
    
    
    // 5. Comparabilità e prodotto logico
    std::vector<std::pair<std::string, std::string>> comparabilities;
    comparabilities.reserve(elements_size * 4); // Euristica dimensionale
    
    for (std::size_t k = 0; k < elements_size; ++k) {
        for (std::size_t h = k + 1; h < elements_size; ++h) {
            auto a1 = indices[k].a1; auto b1 = indices[h].a1;
            auto a2 = indices[k].a2; auto b2 = indices[h].a2;
            
            bool p1_leq = (a1 == b1) || p1.IsLessOrEqual(a1, b1);
            bool p2_leq = (a2 == b2) || p2.IsLessOrEqual(a2, b2);
            
            if (p1_leq && p2_leq) {
                comparabilities.emplace_back(elements_str[k], elements_str[h]);
            } else {
                bool p1_geq = (a1 == b1) || p1.IsLessOrEqual(b1, a1);
                bool p2_geq = (a2 == b2) || p2.IsLessOrEqual(b2, a2);
                if (p1_geq && p2_geq) {
                    comparabilities.emplace_back(elements_str[h], elements_str[k]);
                }
            }
        }
    }
    
    // 6. Configurazione Base Class
    result->FillBaseAttributes(elements_str, comparabilities, nullptr, false);
    
    return result;
}

// ***********************************************
// ComputeMRPLex
// ***********************************************

Tensor<double, 2> ProductLinearPOSet::ComputeMRPLex() const {
    const std::uint64_t n = this->size();
    Tensor<double, 2> matrice(std::array<std::uint64_t, 2>{n, n}, 0.0);
    
    for (std::uint64_t k = 0; k < n; ++k) {
        matrice(k, k) = 1.0;
    }
    
    const std::uint64_t k_size = posets_.size();
    
    for (std::uint64_t v1 = 0; v1 < n; ++v1) {
        for (std::uint64_t v2 = v1 + 1; v2 < n; ++v2) {
            std::uint64_t k1_1 = 0;
            std::uint64_t k1_2 = 0;
            std::uint64_t k2 = 0;
            
            for (std::uint64_t p = 0; p < k_size; ++p) {
                auto v1_p = elements_to_original_[v1][p];
                auto v2_p = elements_to_original_[v2][p];
                
                if (v1_p == v2_p) {
                    ++k2;
                } else if (posets_[p]->IsLessOrEqual(v1_p, v2_p)) {
                    ++k1_1;
                } else if (posets_[p]->IsLessOrEqual(v2_p, v1_p)) {
                    ++k1_2;
                }
            }
            
            std::uint64_t k2_fatt = k2;
            std::uint64_t k2_termine = k2 > 0 ? k2 - 1 : 0;
            std::uint64_t k_fatt = k_size - 1;
            std::uint64_t k_termine = k_size - 2;
            std::uint64_t somma = 1;
            
            for (std::uint64_t s = 1; s <= k2; ++s) {
                somma += k2_fatt * k_fatt;
                k2_fatt *= k2_termine;
                k_fatt *= k_termine;
                if (k2_termine > 0) --k2_termine;
                if (k_termine > 0) --k_termine;
            }
            
            matrice(v1, v2) = somma * (static_cast<double>(k1_1)) / static_cast<double>(k_size);
            matrice(v2, v1) = somma * (static_cast<double>(k1_2)) / static_cast<double>(k_size);
        }
    }
    
    return matrice;
}
