
// MSVC: explicit standard includes (not pulled in transitively).
#include <algorithm>
#include <cstddef>
/**
 * @file linear_extension_generator_tree_of_ideals.cpp
 * @brief Implementation of the Tree of Ideals linear extension generator.
 */


#include "linear_extension_generator_tree_of_ideals.h"

#include <cstdint>
#include <string>
#include "lattice_of_ideals.h"
#include "my_exception.h"


LEGTreeOfIdeals::LEGTreeOfIdeals(std::uint64_t size, const LatticeOfIdeals& lattice)
: LinearExtensionGenerator(size),
lattice_of_ideals_(lattice) {
    
    lattice_of_ideals_crossing_.assign(size, 0);
    more_crossing_.assign(size, 0);
    out_result_buffer_.assign(size, 0);
}

void LEGTreeOfIdeals::Start(std::uint64_t /* unused_id */) {
    // Inizializza il vettore di attraversamento tutto a 0
    std::fill(lattice_of_ideals_crossing_.begin(), lattice_of_ideals_crossing_.end(), 0);
    
    // Ottiene la prima estensione lineare dal path iniziale
    lattice_of_ideals_.GetFromPath(lattice_of_ideals_crossing_,
                                   more_crossing_,
                                   out_result_buffer_);
    
    // riversa nel buffer protetto ottimizzato
    for (std::size_t k = 0; k < out_result_buffer_.size(); ++k) {
        current_linear_extension_.Set(k, out_result_buffer_[k]);
    }
    
    started_ = true;
    current_number_le_ = 1;
}

void LEGTreeOfIdeals::Next() {
    if (!started_) [[unlikely]] {
        throw MyException("LEGTreeOfIdeals error: Start() non eseguito.");
    }
    
    std::int64_t i = static_cast<std::int64_t>(more_crossing_.size()) - 1;
    
    // Risale l'albero finché trova bivi esauriti (more_crossing_ == 0)
    while (i >= 0 && more_crossing_[i] == 0) {
        lattice_of_ideals_crossing_[i] = 0; // Reset dei branch esplorati
        i--;
    }
    
    // Se siamo risaliti oltre la radice, l'esplorazione esaustiva è finita
    if (i < 0) [[unlikely]] {
        throw MyException("LEGTreeOfIdeals error: No more linear extensions available.");
    }
    
    // 2. Avanza di un passo nel branch individuato
    lattice_of_ideals_crossing_[i]++;
    
    // 3. Genera la nuova estensione calcolando il percorso
    lattice_of_ideals_.GetFromPath(lattice_of_ideals_crossing_,
                                   more_crossing_,
                                   out_result_buffer_);
    
    // 4. Copia sul buffer finale in O(1)
    for (std::size_t k = 0; k < out_result_buffer_.size(); ++k) {
        current_linear_extension_.Set(k, out_result_buffer_[k]);
    }
    
    current_number_le_++;
}

bool LEGTreeOfIdeals::HasNext() {
    if (!started_) return false;
    
    // Controlla se esiste almeno una diramazione non ancora esplorata
    for (std::uint8_t m : more_crossing_) {
        if (m != 0) return true;
    }
    
    return false;
}
