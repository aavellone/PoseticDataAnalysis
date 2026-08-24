
// MSVC: explicit standard includes (not pulled in transitively).
#include <string>
#include <string_view>
#include <vector>
#include <memory>
#include <utility>
#include <cstdint>
#include <cstddef>
#include <format>
#include "product_poset.h"
#include "linear_extension_generator.h"

// ***********************************************
// Build - Generazione del Prodotto Cartesiano
// ***********************************************
std::unique_ptr<ProductPOSet> ProductPOSet::Build(const POSet& p1, const POSet& p2) {
    // Inizializzazione con ownership esclusiva
    std::unique_ptr<ProductPOSet> result(new ProductPOSet());
    
    // 1. Eredita la storia dei POSet compositi (usando puntatori raw non-owning const POSet*)
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
    
    const std::uint64_t elements_size = p1.size() * p2.size();
    
    // 2. Pre-allocazione per prestazioni estreme
    std::vector<std::string> elements_str;
    elements_str.reserve(elements_size);
    result->elements_to_original_.reserve(elements_size);
    
    struct Indices { std::uint64_t a1, a2; };
    std::vector<Indices> indices;
    indices.reserve(elements_size);
    
    const auto* pp1_composite = dynamic_cast<const FromPOSets*>(&p1);
    const auto* pp2_composite = dynamic_cast<const FromPOSets*>(&p2);
    
    for (std::uint64_t i = 0; i < p1.size(); ++i) {
        // 1. Usiamo string_view: name1 non viene allocata qui, è solo un puntatore
        std::string_view name1 = p1.GetElementName(i);
        
        // Ottimizzazione opzionale: pre-cacciamo il riferimento alla storia di p1
        const std::vector<std::uint64_t>* p1_hist_ptr = pp1_composite ?
        &pp1_composite->elements_to_original_[i] : nullptr;
        
        for (std::uint64_t j = 0; j < p2.size(); ++j) {
            // 2. Costruzione Nome Ottimizzata
            std::string_view name2 = p2.GetElementName(j);
            
            // Usiamo emplace_back con std::format o concatenazione manuale
            // per minimizzare le riallocazioni della stringa finale
            elements_str.emplace_back(std::format("{}{}{}", name1, FromPOSets::kSep, name2));
            
            indices.push_back({i, j});
            
            // 3. Risoluzione Storia (HPC Style)
            std::vector<std::uint64_t> orig_eids;
            orig_eids.reserve(result->posets_.size()); // Dimensione nota, allocazione singola
            
            // Aggiungiamo i componenti di p1
            if (p1_hist_ptr) {
                orig_eids.insert(orig_eids.end(), p1_hist_ptr->begin(), p1_hist_ptr->end());
            } else {
                orig_eids.push_back(i);
            }
            
            // Aggiungiamo i componenti di p2
            if (pp2_composite) {
                const auto& p2_hist = pp2_composite->elements_to_original_[j];
                orig_eids.insert(orig_eids.end(), p2_hist.begin(), p2_hist.end());
            } else {
                orig_eids.push_back(j);
            }
            
            // Move semantics: il vector orig_eids viene "rubato" e messo nel risultato
            result->elements_to_original_.push_back(std::move(orig_eids));
        }
    }
    
    // 4. Generazione Comparabilità
    std::vector<std::pair<std::string, std::string>> comparabilities;
    // Riserviamo una dimensione stimata ragionevole per evitare continue allocazioni
    comparabilities.reserve(elements_size * 4);
    
    for (std::size_t k = 0; k < elements_size; ++k) {
        for (std::size_t h = k + 1; h < elements_size; ++h) {
            auto a1 = indices[k].a1; auto a2 = indices[k].a2;
            auto b1 = indices[h].a1; auto b2 = indices[h].a2;
            
            bool p1_leq = (a1 == b1) || p1.IsLessOrEqual(a1, b1);
            bool p2_leq = (a2 == b2) || p2.IsLessOrEqual(a2, b2);
            
            if (p1_leq && p2_leq) {
                comparabilities.emplace_back(elements_str[k], elements_str[h]);
            } else {
                // Controllo inverso se non sono comparabili nella prima direzione
                bool p1_geq = (a1 == b1) || p1.IsLessOrEqual(b1, a1);
                bool p2_geq = (a2 == b2) || p2.IsLessOrEqual(b2, a2);
                
                if (p1_geq && p2_geq) {
                    comparabilities.emplace_back(elements_str[h], elements_str[k]);
                }
            }
        }
    }
    
    // 5. Inizializzazione della Classe Base (CSR)
    result->FillBaseAttributes(elements_str, comparabilities, nullptr, false);
    
    return result;
}
