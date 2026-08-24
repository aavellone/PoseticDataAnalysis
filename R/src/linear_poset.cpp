
// MSVC: explicit standard includes (not pulled in transitively).
#include <string>
#include <vector>
#include <memory>
#include <utility>
#include <cstdint>
#include "linear_poset.h"
#include "tensor.h"

// ***********************************************
// Factory Method
// ***********************************************
std::unique_ptr<LinearPOSet> LinearPOSet::Build(const std::vector<std::string>& elements) {
    // Usiamo 'new' perché il costruttore è privato e make_unique fallirebbe
    std::unique_ptr<LinearPOSet> r(new LinearPOSet());
    
    // Genera gli archi della catena lineare
    auto comparabilities = BuildComparability(elements);
    
    // Passa nullptr al parametro Score* (come previsto da FillBaseAttributes)
    // e 'true' per la chiusura transitiva
    r->FillBaseAttributes(elements, comparabilities, nullptr, true);
    
    return r;
}

// ***********************************************
// Generazione archi (catena A -> B -> C)
// ***********************************************
std::vector<std::pair<std::string, std::string>>
LinearPOSet::BuildComparability(const std::vector<std::string>& elements) {
    std::vector<std::pair<std::string, std::string>> comparabilities;
    
    if (elements.size() < 2) {
        return comparabilities;
    }
    
    // Riserviamo la memoria esatta (N - 1 archi per N nodi in una linea)
    comparabilities.reserve(elements.size() - 1);
    
    for (std::uint64_t k = 0; k < elements.size() - 1; ++k) {
        comparabilities.emplace_back(elements[k], elements[k + 1]);
    }
    
    return comparabilities;
}

// ***********************************************
// Matrice di Incidenza
// ***********************************************
Tensor<std::uint8_t, 2> LinearPOSet::IncidenceMatrix() const {
    const std::uint64_t n = size();
    
    // Inizializza una matrice NxN con tutti zeri
    Tensor<std::uint8_t, 2> result({n, n}, 0);
    
    // Essendo un ordine lineare totale, l'elemento i è < di tutti gli elementi j > i.
    // Riempiamo la parte strettamente triangolare superiore con 1.
    for (std::uint64_t row = 0; row < n; ++row) {
        for (std::uint64_t col = row + 1; col < n; ++col) {
            result(row, col) = 1;
        }
    }
    
    return result;
}
