// ================================================================
// File: bucket.h
// ================================================================

#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <cstddef>

#include "bit_set.h"

#include <string>
#include <stdexcept>
#include <cstdint>
#include <memory>
#include <format>
#include <iterator>

/**
 * @class Bucket
 * @brief Insieme a capacità fissa basato su BitSet per alte prestazioni HPC.
 * * @details La classe maschera un BitSet per fornire un'interfaccia compatibile
 * con la STL (iteratori, size, empty), ma operando a livello hardware (O(1)
 * o O(N/64) per le operazioni). È move-only per prevenire copie costose.
 */
class Bucket {
private:
    std::uint64_t max_size_; ///< Capacità massima del bucket (dominio degli elementi).
    BitSet dati_;            ///< Struttura dati sottostante (array di bit).
    
public:
    using iterator       = BitSet::iterator;       ///< Iteratore forward.
    using const_iterator = BitSet::const_iterator; ///< Iteratore forward costante.
    using value_type     = BitSet::word_t;         ///< Tipo di base degli elementi.
    
    /**
     * @brief Costruttore del Bucket.
     * @param max_size Valore massimo inseribile (esclusivo).
     */
    explicit Bucket(std::uint64_t max_size)
    : max_size_(max_size), dati_(max_size) {}
    
    Bucket(Bucket&&) noexcept = default;
    Bucket& operator=(Bucket&&) noexcept = default;
    Bucket(const Bucket&) = delete;
    Bucket& operator=(const Bucket&) = delete;
    
    /** @brief Restituisce l'iteratore all'inizio del bucket. */
    [[nodiscard]] iterator begin() { return dati_.begin(); }
    /** @brief Restituisce l'iteratore alla fine del bucket. */
    [[nodiscard]] iterator end()   { return dati_.end(); }
    /** @brief Restituisce l'iteratore costante all'inizio del bucket. */
    [[nodiscard]] const_iterator begin() const { return dati_.begin(); }
    /** @brief Restituisce l'iteratore costante alla fine del bucket. */
    [[nodiscard]] const_iterator end() const   { return dati_.end(); }
    
    /** @brief Restituisce il numero di elementi inseriti nel bucket (popcount). */
    [[nodiscard]] std::size_t size() const { return dati_.count(); }
    /** @brief Alias di size(). */
    [[nodiscard]] std::size_t count() const { return dati_.count(); }
    /** @brief Verifica se il bucket è vuoto. */
    [[nodiscard]] bool empty() const { return dati_.empty(); }
    
    /**
     * @brief Verifica la presenza di un elemento in O(1).
     * @param val L'elemento da cercare.
     * @return true se l'elemento è presente, false altrimenti.
     */
    [[nodiscard]] inline bool contains(std::uint64_t val) const { return dati_.test(val); }
    
    /**
     * @brief Inserisce un elemento nel bucket in O(1).
     * @param val L'elemento da inserire.
     */
    void insert(std::uint64_t val) { dati_.set(val); }
    
    /**
     * @brief Calcola la differenza insiemistica (*this \ b).
     * @param b Il bucket sottraendo.
     * @return Un nuovo Bucket contenente il risultato.
     * @throws std::invalid_argument Se i bucket hanno capacità diverse.
     */
    [[nodiscard]] Bucket set_difference(const Bucket& b) const {
        if (max_size_ != b.max_size_) {
            throw std::invalid_argument("Buckets must have the same capacity");
        }
        Bucket result(max_size_);
        result.dati_.SetDifference(this->dati_, b.dati_);
        return result;
    }
    
    /**
     * @brief Serializza il contenuto del bucket.
     * @return Una stringa formattata, es: "{1 3 4}"
     */
    [[nodiscard]] std::string to_string() const {
        std::string r = "{";
        bool first = true;
        for (value_type v : dati_) {
            if (!first) r += ' ';
            r += std::to_string(v);
            first = false;
        }
        r += '}';
        return r;
    }
}; // <-- Questa chiusura mancava nel tuo file!
