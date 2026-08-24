/**
 * @file loss_function_mrp.h
 * @brief Interfacce base per le funzioni di perdita MRP (Mutual Ranking Probability).
 *
 * @details Questo file contiene le dichiarazioni per l'architettura base delle Loss Function.
 * Le strutture evitano la proprietà condivisa e utilizzano layout di memoria contigui
 * (`std::vector` per i pesi) per garantire un prefetching ottimale in ambienti multi-core.
 *
 * @par Riferimenti Bibliografici Generali
 * - De Loof, K., De Baets, B., & De Meyer, H. (2006). "Properties of mutual ranking probabilities
 * in partially ordered sets". in *European Journal of Operational Research*.
 *
 * @par Requisiti
 * Compilatore C++20 o superiore.
 */


#pragma once

#include <vector>
#include <unordered_map>
#include <string>
#include <cstdint>
#include <utility>
#include "tensor.h"

/**
 * @class LossFunctionMRP
 * @brief Classe base astratta per le funzioni di perdita su matrici MRP.
 */
class LossFunctionMRP {
protected:
    /// @brief Contatore delle chiamate di valutazione.
    std::uint64_t calls_{0};
    
    /// @brief Riferimento raw al vettore di matrici MRP di riferimento (evita l'overhead di shared_ptr).
    const std::vector<Tensor<double, 2>>* mrps_ref_;
    
    /// @brief Vettore contiguo dei pesi per iterazioni cache-friendly.
    std::vector<std::pair<std::uint64_t, double>> fast_weights_;
    
public:
    /**
     * @brief Costruttore base HPC. Internalizza i pesi per iterazioni O(1).
     * @param p Riferimento al vettore delle matrici MRP reali.
     * @param w Mappa dei pesi dei nodi (viene appiattita internamente).
     */
    LossFunctionMRP(const std::vector<Tensor<double, 2>>& p,
                    const std::unordered_map<std::uint64_t, double>& w);
    
    virtual ~LossFunctionMRP() = default;
    
    /**
     * @brief Fornisce statistiche interne sulle chiamate sotto forma di stringa.
     */
    [[nodiscard]] virtual std::string to_string() const;
    
    /**
     * @brief Operatore principale per il calcolo della Loss.
     * @param x Matrice MRP target da valutare.
     * @return Il valore scalare della funzione di perdita.
     */
    virtual double operator()(const Tensor<double, 2>& x) = 0;
};

/**
 * @class LossFunctionMRPV2
 * @brief Versione avanzata delle funzioni di perdita per la gestione di mappature dinamiche degli indici.
 * @details Elimina l'uso di `std::list` a favore di `std::vector` di vettori per il tracciamento
 * delle trasformazioni lineari degli indici (estensioni).
 */
class LossFunctionMRPV2 : public LossFunctionMRP {
public:
    using LossFunctionMRP::LossFunctionMRP;
    virtual ~LossFunctionMRPV2() = default;
    
    virtual double operator()(const Tensor<double, 2>& x) = 0;
    
    /**
     * @brief Valutazione ottimizzata con mappatura non-lineare degli indici.
     */
    virtual double operator()(const Tensor<double, 2>& x,
                              const std::vector<std::vector<std::uint_fast64_t>>& rows_functions,
                              const std::vector<std::vector<std::uint_fast64_t>>& cols_functions) = 0;
    
    /**
     * @brief Valutazione che restituisce vettori di risultati parziali per elemento.
     */
    virtual void operator()(const Tensor<double, 2>& x,
                            const std::vector<std::vector<std::uint_fast64_t>>& rows_functions,
                            const std::vector<std::vector<std::uint_fast64_t>>& cols_functions,
                            std::vector<std::pair<std::uint_fast64_t, double>>& result) = 0;
};
