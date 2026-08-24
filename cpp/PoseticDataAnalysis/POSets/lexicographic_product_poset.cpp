
// MSVC: explicit standard includes (not pulled in transitively).
#include <string>
#include <string_view>
#include <vector>
#include <memory>
#include <utility>
#include <cstdint>
#include <cstddef>
#include <format>
#include "lexicographic_product_poset.h"

#include <stdexcept>

// ***********************************************
// Build - Generazione del Prodotto Lessicografico
// ***********************************************

std::unique_ptr<LexicographicProductPOSet> LexicographicProductPOSet::Build(const POSet& p1, const POSet& p2) {
    std::unique_ptr<LexicographicProductPOSet> result(new LexicographicProductPOSet());
    
    const std::uint64_t p1_size = p1.size();
    const std::uint64_t p2_size = p2.size();
    const std::uint64_t elements_size = p1_size * p2_size;
    
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
    
    
    
    
    std::vector<std::pair<std::string, std::string>> comparabilities;
    comparabilities.reserve(elements_size * 4); // euristica
    
    for (std::size_t k = 0; k < elements_size; ++k) {
        for (std::size_t h = k + 1; h < elements_size; ++h) {
            auto a1 = indices[k].a1;
            auto b1 = indices[h].a1;
            auto a2 = indices[k].a2;
            auto b2 = indices[h].a2;
            
            // Logica Lessicografica: (a1, a2) <= (b1, b2) SSE a1 < b1 OPPURE (a1 == b1 E a2 <= b2)
            bool p1_less = (a1 != b1) && p1.IsLessOrEqual(a1, b1);
            bool p1_eq_p2_leq = (a1 == b1) && ((a2 == b2) || p2.IsLessOrEqual(a2, b2));
            
            if (p1_less || p1_eq_p2_leq) {
                comparabilities.emplace_back(elements_str[k], elements_str[h]);
            } else {
                // Controllo inverso (b <= a)
                bool p1_greater = (a1 != b1) && p1.IsLessOrEqual(b1, a1);
                bool p1_eq_p2_geq = (a1 == b1) && ((a2 == b2) || p2.IsLessOrEqual(b2, a2));
                
                if (p1_greater || p1_eq_p2_geq) {
                    comparabilities.emplace_back(elements_str[h], elements_str[k]);
                }
            }
        }
    }
    
    result->FillBaseAttributes(elements_str, comparabilities, nullptr, false);
    
    return result;
}
