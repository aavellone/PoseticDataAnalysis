#pragma once

#include "from_posets.h"

#include <vector>
#include <string>
#include <memory>
#include <cstdint>

/**
 * @class LexicographicProductPOSet
 * @brief Implementazione HPC del Prodotto Lessicografico tra due POSet.
 * @details Nel prodotto lessicografico (A * B), un elemento (a1, a2) precede
 * o è uguale a (b1, b2) se e solo se a1 < b1 (in A), oppure se a1 = b1 e a2 <= b2 (in B).
 */
class LexicographicProductPOSet : public FromPOSets {
private:
    LexicographicProductPOSet() : FromPOSets() {}
    
public:
    /**
     * @brief Costruisce il prodotto lessicografico.
     * @details Evita il round-trip delle stringhe ereditando la storia degli
     * identificatori matematicamente in tempo O(1).
     */
    [[nodiscard]] static std::unique_ptr<LexicographicProductPOSet> Build(const POSet& p1, const POSet& p2);
};
