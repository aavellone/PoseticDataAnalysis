#pragma once

#include "from_posets.h"

#include <vector>
#include <string>
#include <memory>
#include <cstdint>

/**
 * @class ProductPOSet
 * @brief Implementazione HPC del Prodotto Diretto (componente per componente) tra due POSet.
 * @details Nel prodotto diretto (A x B), un elemento (a1, a2) precede
 * o è uguale a (b1, b2) se e solo se a1 <= b1 (in A) e a2 <= b2 (in B).
 */
class ProductPOSet : public FromPOSets {
private:
    ProductPOSet() : FromPOSets() {}

public:
    /**
     * @brief Costruisce il prodotto diretto.
     * @details Evita il round-trip delle stringhe ereditando la storia degli
     * identificatori matematicamente in tempo O(1).
     */
    [[nodiscard]] static std::unique_ptr<ProductPOSet> Build(const POSet& p1,
                                                             const POSet& p2);
};
