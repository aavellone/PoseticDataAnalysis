
// MSVC: explicit standard includes (not pulled in transitively).
#include <memory>
#include <optional>
#include <utility>
#include <cstdint>
/**
 * @file linear_extension_generator_bubley_dyer.cpp
 * @brief Implementazione del generatore di estensioni lineari Bubley-Dyer (MCMC).
 */

#include "linear_extension_generator_bubley_dyer.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

#include "linear_extension_generator.h"
#include "my_exception.h"
#include "poset.h"
#include "random.h"

LEGBubleyDyer::LEGBubleyDyer(POSet* poset, std::unique_ptr<Random> rnd)
: LinearExtensionGenerator(poset->size()),
poset_(poset),
rnd_(std::move(rnd)),
initial_le_(poset->size()) {
    rnd_->InitBubleyDyer(poset->size());
    poset->FirstLE(initial_le_);
}

LEGBubleyDyer::LEGBubleyDyer(POSet* poset, std::unique_ptr<Random> rnd, LinearExtension initial_le)
: LinearExtensionGenerator(poset->size()),
poset_(poset),
rnd_(std::move(rnd)),
initial_le_(std::move(initial_le)) {
    rnd_->InitBubleyDyer(poset->size());

    // Una LE iniziale non valida romperebbe silenziosamente l'invariante della
    // catena MCMC (che preserva la validità passo dopo passo): validazione
    // completa una tantum, costo trascurabile rispetto alla generazione.
    const std::uint64_t n = poset->size();
    if (initial_le_.size() != n) {
        throw MyException("LEGBubleyDyer: initial LE size does not match poset size.");
    }
    for (std::uint64_t pos = 0; pos < n; ++pos) {
        const std::uint64_t el = initial_le_.GetVal(pos);
        if (el >= n || initial_le_.GetPos(el) != pos) {
            throw MyException("LEGBubleyDyer: initial LE is not a coherent permutation.");
        }
    }
    for (std::uint64_t i = 0; i + 1 < n; ++i) {
        for (std::uint64_t j = i + 1; j < n; ++j) {
            if (poset_->GreaterThan(initial_le_.GetVal(i), initial_le_.GetVal(j))) {
                throw MyException("LEGBubleyDyer: initial LE violates the partial order.");
            }
        }
    }
}

void LEGBubleyDyer::Start(std::uint64_t quante) {
    rnd_->Restart();
    max_number_le_ = quante;
    
    for(std::uint64_t i = 0; i < n_elements_; ++i) {
        current_linear_extension_.Set(i, initial_le_.GetVal(i));
    }
    
    // Nelle classi HPC C++ queste variabili dovrebbero essere rinominate
    // con snake_case e underscore finale nel tuo .h (es: to_update_)
    to_update_ = false;
    is_switched_ = false;
    position_to_update_ = 0;
    
    started_ = true;
    current_number_le_ = 1;
}

void LEGBubleyDyer::Next() {
    if (!started_) {
        throw MyException("LEGBubleyDyer error: not started yet!");
    }
    
    ++current_number_le_;

    to_update_ = rnd_->RndUpdate();
    is_switched_ = false;

    // Con 0 o 1 elementi esiste una sola estensione lineare: nessuno swap
    // possibile. La guardia evita anche la lettura fuori dai limiti in
    // GetVal(position_to_update_ + 1) (con n==1, RndPosition() darebbe 0 e si
    // accederebbe alla posizione 1, inesistente).
    if (n_elements_ < 2) {
        return;
    }

    if (to_update_ != 0) {
        position_to_update_ = rnd_->RndPosition();
        
        std::uint64_t e1 = current_linear_extension_.GetVal(position_to_update_);
        std::uint64_t e2 = current_linear_extension_.GetVal(position_to_update_ + 1);
        
        if (!poset_->GreaterThan(e2, e1)) {
            current_linear_extension_.Set(position_to_update_, e2);
            current_linear_extension_.Set(position_to_update_ + 1, e1);
            is_switched_ = true;
        }
    }
}

std::uint64_t LEGBubleyDyer::EvaluateNumberOfIteration(double error) const {
    std::uint64_t nelementi = current_linear_extension_.size();
    
    double add1 = std::pow(nelementi, 4) * std::pow(std::log(nelementi), 2);
    double add2 = std::pow(nelementi, 3) * std::log(nelementi) * std::log(1.0 / error);
    
    return static_cast<std::uint64_t>(std::max(add1 + add2, 1.0));
}

bool LEGBubleyDyer::UpdateCounters(std::optional<std::uint64_t> quante, std::optional<double> error) {
    if (quante.has_value() && quante.value() != 0) {
        max_number_le_ += quante.value();
        return true;
    } else if (!quante.has_value()) {
        auto n = this->EvaluateNumberOfIteration(error.value());
        if (max_number_le_ < n) {
            max_number_le_ = n;
            return true;
        }
    }
    return false;
}
