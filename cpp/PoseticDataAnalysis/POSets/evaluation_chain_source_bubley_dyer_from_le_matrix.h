/**
 * @file evaluation_chain_source_bubley_dyer_from_le_matrix.h
 * @brief Sorgente lazy di catene di valutazione basate su LEGBubleyDyer.
 *
 * @par Requirements
 *   C++20 or later (`-std=c++20`)
 *
 * @par Naming conventions (Google C++ Style Guide)
 *   - Types / classes:              PascalCase
 *   - Methods:                      PascalCase
 *   - Data members:                 snake_case_  (trailing underscore)
 *   - Local variables / parameters: snake_case
 */

#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <cstddef>

#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

#include "function_linear_extension.h"
#include "linear_extension.h"
#include "linear_extension_generator_bubley_dyer.h"
#include "my_exception.h"
#include "poset.h"
#include "random.h"
#include "tensor.h"

/**
 * @class ECSBubleyDyerFromLEMatrix
 * @brief Sorgente di EvaluationChainTask: una catena Bubley-Dyer per punto iniziale.
 *
 * @details Implementa il concept EvaluationChainSource. Riceve nel costruttore
 * i punti di partenza (una LinearExtension per catena, tipicamente le colonne
 * di una matrice fornita dall'utente), i seed e i prototipi di funzioni e
 * tensori. Ogni Next() costruisce LA SOLA catena successiva: generatore
 * LEGBubleyDyer con seed e LE iniziale propri (validata dal costruttore del
 * generatore e già Start()-ato), cloni delle funzioni e tensori locali
 * azzerati. La costruzione lazy, combinata con il merge incrementale di
 * POSet::evaluation_parallel, limita la memoria di picco a O(n_threads)
 * catene invece di O(K).
 *
 * @warning HasNext()/Next() NON sono thread-safe: evaluation_parallel li
 * invoca sotto lock. La classe non tocca MAI l'API di R: i punti iniziali
 * devono arrivare già convertiti in strutture C++ (l'estrazione dai SEXP va
 * fatta sul main thread, l'API di R non è thread-safe).
 *
 * @note I vettori di prototipi (funzioni e tensori) sono riferiti, non
 * copiati: devono restare vivi per tutta la vita della sorgente.
 */
class ECSBubleyDyerFromLEMatrix final {
public:
    using GeneratorType = LEGBubleyDyer;

    /**
     * @brief Costruisce la sorgente e valida i parametri (fail-fast sul main thread).
     * @param poset             POSet target (non owning, non nullo).
     * @param initial_les       Punti di partenza, uno per catena (consumati via move).
     * @param seeds             Un seed per catena (stessa lunghezza di initial_les).
     * @param count_per_chain   LE da campionare per catena (> 0), passato a Start().
     * @param fle_prototypes    Funzioni prototipo, clonate una volta per catena.
     * @param result_prototypes Tensori di output: definiscono numero e shape dei
     *                          tensori locali (azzerati) di ogni catena.
     * @throws MyException se i parametri sono incoerenti o se una funzione non è
     *         thread-safe (es. funzioni definite in R) e non può girare in parallelo.
     */
    ECSBubleyDyerFromLEMatrix(POSet* poset,
                  std::vector<LinearExtension> initial_les,
                  std::vector<std::uint64_t> seeds,
                  std::uint64_t count_per_chain,
                  const std::vector<std::unique_ptr<FunctionLinearExtension>>& fle_prototypes,
                  const std::vector<std::unique_ptr<Tensor<double, 2>>>& result_prototypes)
    : poset_(poset),
      initial_les_(std::move(initial_les)),
      seeds_(std::move(seeds)),
      count_per_chain_(count_per_chain),
      fle_prototypes_(&fle_prototypes),
      result_prototypes_(&result_prototypes)
    {
        if (poset_ == nullptr) {
            throw MyException("ECSBubleyDyerFromLEMatrix: 'poset' cannot be null.");
        }
        if (initial_les_.empty()) {
            throw MyException("ECSBubleyDyerFromLEMatrix: at least one initial linear extension is required.");
        }
        if (seeds_.size() != initial_les_.size()) {
            throw MyException("ECSBubleyDyerFromLEMatrix: seeds/initial linear extensions size mismatch.");
        }
        if (count_per_chain_ == 0) {
            throw MyException("ECSBubleyDyerFromLEMatrix: count_per_chain must be > 0.");
        }
        if (fle_prototypes.size() != result_prototypes.size()) {
            throw MyException("ECSBubleyDyerFromLEMatrix: fles/results size mismatch.");
        }
        for (const auto& fle : fle_prototypes) {
            if (!fle->IsThreadSafe()) {
                throw MyException("ECSBubleyDyerFromLEMatrix: a function is not thread-safe "
                                  "(e.g. R-defined functions) and cannot run in parallel chains.");
            }
        }
    }

    /// @brief true se esiste una catena successiva da produrre.
    [[nodiscard]] bool HasNext() const noexcept {
        return next_ < initial_les_.size();
    }

    /**
     * @brief Costruisce e avvia la catena successiva.
     * @return Il task completo (generatore Start()-ato, cloni, tensori azzerati).
     * @throws MyException se la sorgente è esaurita; propaga le eccezioni di
     *         validazione della LE iniziale (costruttore di LEGBubleyDyer).
     */
    [[nodiscard]] EvaluationChainTask<LEGBubleyDyer> Next() {
        if (!HasNext()) {
            throw MyException("ECSBubleyDyerFromLEMatrix::Next: source exhausted.");
        }
        const std::size_t k = next_++;

        EvaluationChainTask<LEGBubleyDyer> task;
        task.leg = std::make_unique<LEGBubleyDyer>(
            poset_, std::make_unique<Random>(seeds_[k]), std::move(initial_les_[k]));
        task.leg->Start(count_per_chain_);

        task.fles.reserve(fle_prototypes_->size());
        for (const auto& fle : *fle_prototypes_) {
            task.fles.push_back(fle->Clone());
        }

        task.results.reserve(result_prototypes_->size());
        for (const auto& result : *result_prototypes_) {
            task.results.push_back(std::make_unique<Tensor<double, 2>>(result->Shape(), 0.0));
        }
        return task;
    }

    /// @brief Numero totale di catene che la sorgente produrrà.
    [[nodiscard]] std::uint64_t TotalChains() const noexcept {
        return initial_les_.size();
    }

private:
    POSet* poset_;
    std::vector<LinearExtension> initial_les_;
    std::vector<std::uint64_t> seeds_;
    std::uint64_t count_per_chain_;
    const std::vector<std::unique_ptr<FunctionLinearExtension>>* fle_prototypes_;
    const std::vector<std::unique_ptr<Tensor<double, 2>>>* result_prototypes_;
    std::size_t next_ = 0;
};

static_assert(EvaluationChainSource<ECSBubleyDyerFromLEMatrix>,
              "ECSBubleyDyerFromLEMatrix must satisfy the EvaluationChainSource concept.");
