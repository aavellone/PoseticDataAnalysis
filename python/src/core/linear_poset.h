#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <cstdint>

#include "tensor.h"
#include "poset.h"

#include <vector>
#include <string>
#include <memory>
#include <utility>

class LinearPOSet : public POSet {
private:
    // Costruttore privato
    LinearPOSet() : POSet() {}
    
    // Metodo helper statico per costruire la catena
    static std::vector<std::pair<std::string, std::string>> BuildComparability(
             const std::vector<std::string>& elements);
    
public:
    // Factory method che restituisce un unique_ptr
    static std::unique_ptr<LinearPOSet> Build(const std::vector<std::string>& elements);
    
    // Override del metodo della classe base utilizzando Matrice
    [[nodiscard]] Tensor<std::uint8_t, 2> IncidenceMatrix() const override;
};
