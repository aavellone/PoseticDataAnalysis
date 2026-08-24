#pragma once

#include "from_posets.h"
#include "tensor.h"

#include <vector>
#include <string>
#include <memory>
#include <cstdint>

/**
 * @class ProductLinearPOSet
 * @brief Implementazione HPC del Prodotto Cartesiano per POSet lineari.
 */
class ProductLinearPOSet : public FromPOSets {
private:
    ProductLinearPOSet() : FromPOSets() {}
    
public:
    [[nodiscard]] static std::unique_ptr<ProductLinearPOSet> Build(const POSet& p1, const POSet& p2);
    /**
     * @brief Calcola la matrice MRP lessicografica.
     */
    [[nodiscard]] Tensor<double, 2> ComputeMRPLex() const;
    
};

