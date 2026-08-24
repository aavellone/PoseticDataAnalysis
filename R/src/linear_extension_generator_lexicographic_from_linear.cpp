
// MSVC: explicit standard includes (not pulled in transitively).
#include <vector>
#include <utility>
#include <cstdint>
#include <functional>
/**
 * @file leg_lexicographic_from_linear.cc
 * @brief Implementation of the lexicographical generator in HPC-oriented C++20.
 */

#include "linear_extension_generator_lexicographic_from_linear.h"

#include <algorithm>
#include <numeric>
#include <ranges>
#include <stdexcept>

#include "my_exception.h"

std::uint64_t LEGLexicographicFromLinear::CalculateTotalSize(
    const std::vector<std::uint64_t> &sizes) noexcept {
    return std::reduce(sizes.begin(), sizes.end(), 1ULL,
        std::multiplies<unsigned long long>());
}

LEGLexicographicFromLinear::LEGLexicographicFromLinear(
    std::vector<std::uint64_t> group_sizes)
    : LinearExtensionGenerator(CalculateTotalSize(group_sizes)),
      group_sizes_(std::move(group_sizes))
{
    
    total_elements_ = CalculateTotalSize(group_sizes_);
    const std::uint64_t num_groups = group_sizes_.size();

    if (num_groups == 0) {
        throw MyException("Inconsistent group sizes.");
    }

    priority_.resize(num_groups);
    strides_.resize(num_groups);
          
    // Calcolo degli strides per il row-major order del prodotto cartesiano
    std::uint64_t current_stride = 1;
    for (std::int64_t i = num_groups - 1; i >= 0; --i) {
        strides_[i] = current_stride;
        current_stride *= group_sizes_[i];
    }

    // Pre-calcoliamo l'offset totale di una dimensione per aggiornare l'ID in O(1)
    coords_.resize(num_groups);
    wrap_offsets_.resize(num_groups);
    for (std::uint64_t i = 0; i < num_groups; ++i) {
        wrap_offsets_[i] = (group_sizes_[i] - 1) * strides_[i];
    }
    
    // Pre-calcolo del numero totale di permutazioni lessicografiche (n_groups!)
    total_permutations_ = 1;
    for (std::uint64_t i = 2; i <= num_groups; ++i) {
        total_permutations_ *= i;
    }
}

void LEGLexicographicFromLinear::Start(std::uint64_t /*unused_id*/) {
    // Inizializza il vettore delle priorità come sequenza identità (0, 1, 2...)
    std::iota(priority_.begin(), priority_.end(), 0ULL);

    // Genera la prima estensione lineare
    RefreshData();

    this->started_ = true;
    this->current_number_le_ = 1;
}

void LEGLexicographicFromLinear::Next() {
    if (!this->started_) {
        throw MyException("Generator not started.");
    }

    // Genera la permutazione lessicografica successiva
    auto result = std::ranges::next_permutation(priority_);

    if (!result.found) {
        throw MyException("No further lexicographical extensions.");
    }

    RefreshData();
    this->current_number_le_++;
}

bool LEGLexicographicFromLinear::HasNext() noexcept {
    if (!this->started_)
        return false;

    // Se il vettore priority_ non è ordinato in modo decrescente, esiste una
    // permutazione successiva
    return !std::ranges::is_sorted(priority_, std::greater{});
}

std::uint64_t LEGLexicographicFromLinear::NumberOfLe() const noexcept {
    return total_permutations_;
}
/*
void LEGLexicographicFromLinear::RefreshData() noexcept {
    const std::uint64_t num_groups = group_sizes_.size();
    
    // Vettore per tenere traccia delle coordinate correnti nella griglia multidimensionale
    std::vector<std::uint64_t> coords(num_groups, 0);
    
    for (std::uint64_t step = 0; step < total_elements_; ++step) {
        // 1. Calcola l'ID reale (appiattito) partendo dalle coordinate e dagli strides
        std::uint64_t element_id = 0;
        for (std::uint64_t i = 0; i < num_groups; ++i) {
            element_id += coords[i] * strides_[i];
        }
        
        // 2. Salva l'elemento nella sequenza dell'estensione lineare
        this->current_linear_extension_.Set(step, element_id);
        
        // 3. Incrementa le coordinate comportandosi come un contatore "multi-base"
        // L'ordine di incremento è dettato dalla priorità (l'ultimo elemento di priority_ è il loop più interno)
        for (auto it = priority_.rbegin(); it != priority_.rend(); ++it) {
            const std::uint64_t dim = *it;
            coords[dim]++;
            
            // Se non c'è riporto per questa dimensione, fermati
            if (coords[dim] < group_sizes_[dim]) {
                break;
            }
            
            // Altrimenti azzera e propaga il riporto alla dimensione successiva in ordine di priorità
            coords[dim] = 0;
        }
    }
}*/


void LEGLexicographicFromLinear::RefreshData() noexcept {
    
    // Reset veloce dei buffer senza riallocazione
    std::fill(coords_.begin(), coords_.end(), 0ULL);
    std::uint64_t current_id = 0;
    
    for (std::uint64_t step = 0; step < total_elements_; ++step) {
        // Scrittura diretta dell'ID corrente
        this->current_linear_extension_.Set(step, current_id);
        
        // Aggiornamento incrementale delle coordinate (come un contatore meccanico)
        // La priorità decide quale dimensione "gira" più velocemente
        for (auto it = priority_.rbegin(); it != priority_.rend(); ++it) {
            const std::uint64_t dim = *it;
            
            if (++coords_[dim] < group_sizes_[dim]) {
                // Caso comune: incrementiamo solo questa dimensione e l'ID
                current_id += strides_[dim];
                break;
            } else {
                // Caso riporto: azzeriamo questa coordinata e aggiorniamo l'ID
                current_id -= wrap_offsets_[dim];
                coords_[dim] = 0;
                // Il loop prosegue alla prossima dimensione in ordine di priorità
            }
        }
    }
}
