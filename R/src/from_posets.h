/**
 * @file from_posets.h
 * @brief Classe base HPC per i POSet compositi (Prodotti, Somme, ecc.).
 *
 * @details Questa classe traccia la genealogia di un POSet derivato da operazioni
 * su altri POSet. Per massimizzare le prestazioni, memorizza solo puntatori "osservatori"
 * (non-owning raw pointers) agli oggetti di origine, azzerando l'overhead di
 * reference counting (zero-cost abstraction).
 *
 * @author Alessandro Avellone
 * @version 4.0 (HPC Edition - 100% unique_ptr & observer pointers)
 */

#pragma once

#include "poset.h"

#include <vector>
#include <string>
#include <memory>
#include <cstdint>

class FromPOSets : public POSet {
protected:
    // Costruttore accessibile alle classi derivate (ProductPOSet, ecc.)
    FromPOSets() : POSet() {}
    
public:
    static constexpr char kSep = '_';
    
    // HPC Zero-Cost Abstraction: Nessuna ownership. Memorizza solo i riferimenti
    // costanti ai POSet genitrici per tracciare la gerarchia in O(1).
    std::vector<const POSet*> posets_;
    
    // Mappa la cronologia degli identificatori originali
    std::vector<std::vector<std::uint64_t>> elements_to_original_;
    
public:
    /**
     * @brief Factory Method generico per la base dei POSet compositi.
     * @param p1 Primo POSet genitore.
     * @param p2 Secondo POSet genitore.
     * @return unique_ptr al POSet composito di base.
     */
    [[nodiscard]] static std::unique_ptr<FromPOSets> Build(const POSet& p1, const POSet& p2) {
        std::unique_ptr<FromPOSets> result(new FromPOSets());
        
        // Risoluzione della storia del primo genitore
        if (const auto* pp1 = dynamic_cast<const FromPOSets*>(&p1)) {
            result->posets_.insert(result->posets_.end(), pp1->posets_.begin(), pp1->posets_.end());
        } else {
            result->posets_.push_back(&p1);
        }
        
        // Risoluzione della storia del secondo genitore
        if (const auto* pp2 = dynamic_cast<const FromPOSets*>(&p2)) {
            result->posets_.insert(result->posets_.end(), pp2->posets_.begin(), pp2->posets_.end());
        } else {
            result->posets_.push_back(&p2);
        }
        
        return result;
    }
};
