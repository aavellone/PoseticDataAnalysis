/**
 * @file BitSet.h
 * @brief Bitset ad alte prestazioni (HPC) con iteratore forward zero-allocation.
 *
 * @details Implementazione basata su array di word a 64 bit con storage row-major.
 * Le operazioni insiemistiche sfruttano istruzioni SIMD/hardware (popcount, countr_zero).
 *
 * @par Requirements
 * C++20 or later (`-std=c++20`)
 *
 * @par Naming conventions (Google C++ Style Guide)
 *   - Types / structs:              PascalCase
 *   - Methods:                      PascalCase, eccetto quelli che replicano
 *                                   l'interfaccia STL (set, test, flip, count,
 *                                   empty, capacity, begin, end, …) che restano snake_case
 *   - Data members:                 snake_case_
 *   - Local variables / parameters: snake_case
 *   - Constants (constexpr):        kCamelCase
 */

#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <string>
#include <utility>
#include <cstddef>

#include <algorithm>
#include <bit>
#include <concepts>
#include <cstdint>
#include <cstring>
#include <memory>
#include <vector>

/**
 * @class BitSet
 * @brief Insieme di bit a larghezza fissa ottimizzato per ambienti HPC.
 *
 * @details Lo storage è un array flat di word a 64 bit allocato sul heap.
 * L'accesso ai singoli bit avviene tramite shift/mask in O(1);
 * le operazioni insiemistiche (AND, OR, differenza) operano word-by-word in O(N/64).
 */
class BitSet {
public:
    /// @brief Tipo della singola parola di storage (64 bit).
    using word_t = uint64_t;
    
    /// @brief Numero di bit contenuti in una singola word.
    static constexpr uint64_t kBitsPerWord = 64;
    
private:
    uint64_t                  capacity_; ///< Numero di bit logici nel set.
    uint64_t                  nblock_;   ///< Numero di word_t nel backing store.
    std::unique_ptr<word_t[]> store_;    ///< Array flat di word; ownership esclusiva.
    
public:
    // =========================================================================
    // Costruttori e gestione memoria
    // =========================================================================
    
    /**
     * @brief Costruisce un BitSet con @p capacity bit, tutti inizializzati a 0.
     * @param capacity Numero di bit logici da allocare.
     */
    explicit BitSet(uint64_t capacity)
    : capacity_(capacity),
    nblock_((capacity + kBitsPerWord - 1) / kBitsPerWord),
    store_(std::make_unique<word_t[]>(nblock_)) {
        std::memset(store_.get(), 0, nblock_ * sizeof(word_t));
    }
    
    /**
     * @brief Copy constructor: copia profonda del backing store.
     * @param other BitSet sorgente.
     */
    BitSet(const BitSet& other)
    : capacity_(other.capacity_),
    nblock_(other.nblock_),
    store_(std::make_unique<word_t[]>(other.nblock_)) {
        std::memcpy(store_.get(), other.store_.get(), nblock_ * sizeof(word_t));
    }
    
    /**
     * @brief Copy assignment: riallocca il backing store solo se le dimensioni differiscono.
     * @param other BitSet sorgente.
     * @return Riferimento a *this.
     */
    BitSet& operator=(const BitSet& other) {
        if (this != &other) {
            if (nblock_ != other.nblock_) {
                store_ = std::make_unique<word_t[]>(other.nblock_);
            }
            capacity_ = other.capacity_;
            nblock_   = other.nblock_;
            std::memcpy(store_.get(), other.store_.get(), nblock_ * sizeof(word_t));
        }
        return *this;
    }
    
    /**
     * @brief Move constructor: trasferisce ownership del backing store in O(1).
     * @param other BitSet sorgente (lasciato in stato vuoto).
     */
    BitSet(BitSet&& other) noexcept
    : capacity_(other.capacity_),
    nblock_(other.nblock_),
    store_(std::move(other.store_)) {
        other.capacity_ = 0;
        other.nblock_   = 0;
    }
    
    /**
     * @brief Move assignment: trasferisce ownership del backing store in O(1).
     * @param other BitSet sorgente (lasciato in stato vuoto).
     * @return Riferimento a *this.
     */
    BitSet& operator=(BitSet&& other) noexcept {
        if (this != &other) {
            capacity_       = other.capacity_;
            nblock_         = other.nblock_;
            store_          = std::move(other.store_);
            other.capacity_ = 0;
            other.nblock_   = 0;
        }
        return *this;
    }
    
    // =========================================================================
    // Accesso ai singoli bit  (interfaccia STL-compatibile → snake_case)
    // =========================================================================
    
    /**
     * @brief Verifica se il bit @p bit è impostato a 1.
     * @param bit Indice del bit (0-based).
     * @return @c true se il bit è 1, @c false se è 0 o fuori range.
     */
    [[nodiscard]] inline bool test(uint64_t bit) const noexcept {
        if (bit >= capacity_) return false;
        return (store_[bit / kBitsPerWord] &
                (static_cast<word_t>(1) << (bit % kBitsPerWord))) != 0;
    }
    
    /**
     * @brief Restituisce una rappresentazione ASCII del BitSet.
     * @details I bit sono mostrati in ordine crescente di indice (0 → capacity-1).
     * @return Stringa di '0' e '1' di lunghezza pari a capacity().
     */
    [[nodiscard]] std::string to_string() const {
        std::string s;
        s.reserve(capacity_);
        for (uint64_t i = 0; i < capacity_; ++i) {
            s += test(i) ? '1' : '0';
        }
        return s;
    }
    
    // =========================================================================
    // Metodi HPC sui singoli bit
    // =========================================================================
    
    /**
     * @brief Imposta a 1 il bit @p b (STL-compatible: snake_case).
     * @param b Indice del bit. @b Nessun bounds-check: UB se @p b >= capacity().
     * @note  b >> 6 ≡ b / 64 ;  b & 63 ≡ b % 64.
     */
    inline void set(uint64_t b) noexcept {
        store_[b >> 6] |= (static_cast<word_t>(1) << (b & 63));
    }
    
    /**
     * @brief Azzera il bit @p b.
     * @param b Indice del bit. @b Nessun bounds-check: UB se @p b >= capacity().
     */
    inline void Unset(uint64_t b) noexcept {
        store_[b >> 6] &= ~(static_cast<word_t>(1) << (b & 63));
    }
    
    /**
     * @brief Azzera tutti i bit del BitSet in un'unica operazione (HPC bulk-clear).
     * @details Utilizza @c std::memset sul backing store: O(N/64).
     */
    inline void Unset() noexcept {
        std::memset(store_.get(), 0, nblock_ * sizeof(word_t));
    }
    
    /**
     * @brief Legge il bit @p b con bounds-check.
     * @param b Indice del bit (0-based).
     * @return @c true se il bit è 1, @c false se è 0 o fuori range.
     */
    [[nodiscard]] inline bool Get(uint64_t b) const noexcept {
        if (b >= capacity_) return false;
        return (store_[b >> 6] & (static_cast<word_t>(1) << (b & 63))) != 0;
    }
    
    /**
     * @brief Inverte tutti i bit dell'insieme (NOT logico) — STL-compatible.
     * @details NOT bit-a-bit su ogni word; i bit di padding oltre capacity_
     *          vengono riazzerati per prevenire falsi positivi.
     */
    inline void flip() noexcept {
        for (uint64_t i = 0; i < nblock_; ++i) {
            store_[i] = ~store_[i];
        }
        const uint64_t remainder = capacity_ % kBitsPerWord;
        if (remainder != 0) {
            const word_t mask  = (static_cast<word_t>(1) << remainder) - 1;
            store_[nblock_ - 1] &= mask;
        }
    }
    
    // =========================================================================
    // Operazioni insiemistiche HPC
    // =========================================================================
    
    /**
     * @brief Intersezione in-place: @c this = this ∩ other  (AND bit-a-bit).
     * @param other Operando destro.
     * @return Riferimento a *this.
     * @note I blocchi di @c this eccedenti la dimensione di @p other vengono azzerati.
     */
    inline BitSet& operator&=(const BitSet& other) noexcept {
        const uint64_t min_blocks = std::min(nblock_, other.nblock_);
        for (uint64_t i = 0; i < min_blocks; ++i) {
            store_[i] &= other.store_[i];
        }
        for (uint64_t i = min_blocks; i < nblock_; ++i) {
            store_[i] = 0;
        }
        return *this;
    }
    
    /**
     * @brief Unione in-place: @c this = this ∪ other  (OR bit-a-bit).
     * @param other Operando destro.
     * @return Riferimento a *this.
     * @note I blocchi di @c this eccedenti la dimensione di @p other restano invariati.
     */
    inline BitSet& operator|=(const BitSet& other) noexcept {
        const uint64_t min_blocks = std::min(nblock_, other.nblock_);
        for (uint64_t i = 0; i < min_blocks; ++i) {
            store_[i] |= other.store_[i];
        }
        return *this;
    }
    
    /**
     * @brief Verifica se @c this ⊆ @p other.
     * @details Implementa (A AND NOT B) == 0 word-by-word.
     * @param other Insieme potenzialmente contenente @c this.
     * @return @c true se ogni bit a 1 in @c this è a 1 anche in @p other.
     */
    [[nodiscard]] inline bool IsSubsetOf(const BitSet& other) const noexcept {
        const uint64_t min_blocks = std::min(nblock_, other.nblock_);
        for (uint64_t i = 0; i < min_blocks; ++i) {
            if ((store_[i] & ~other.store_[i]) != 0) return false;
        }
        for (uint64_t i = min_blocks; i < nblock_; ++i) {
            if (store_[i] != 0) return false;
        }
        return true;
    }
    
    /**
     * @brief Cardinalità della differenza insiemistica @c this \ @p other.
     * @details Calcola popcount(A AND NOT B) per ogni word: O(N/64).
     * @param other Insieme da sottrarre.
     * @return Numero di bit a 1 in @c this che non sono a 1 in @p other.
     */
    [[nodiscard]] inline std::uint64_t CountDifference(
                                                       const BitSet& other) const noexcept {
                                                           std::uint64_t cnt = 0;
                                                           for (std::uint64_t i = 0; i < nblock_; ++i) {
                                                               cnt += std::popcount(store_[i] & ~other.store_[i]);
                                                           }
                                                           return cnt;
                                                       }
    
    /**
     * @brief Cardinalità dell'insieme (numero di bit a 1) — STL-compatible.
     * @details Sfrutta l'istruzione hardware @c POPCNT: O(N/64).
     * @return Numero di bit impostati a 1.
     */
    [[nodiscard]] inline std::uint64_t count() const noexcept {
        std::uint64_t total = 0;
        for (uint64_t i = 0; i < nblock_; ++i) {
            total += static_cast<uint64_t>(std::popcount(store_[i]));
        }
        return total;
    }
    
    // =========================================================================
    // Iteratore
    // =========================================================================
    
    /**
     * @struct Iterator
     * @brief Forward iterator che enumera gli indici dei bit a 1.
     *
     * @details Avanza isolando il bit meno significativo acceso tramite
     * l'espressione @c (b & -b), evitando branch interni al loop principale.
     * Compatibile con i range-for del C++11 e con gli algoritmi STL forward.
     */
    struct Iterator {
        /// @name Alias STL richiesti dai concetti di iteratore
        ///@{
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = uint64_t;
        using pointer           = const uint64_t*;
        using reference         = uint64_t;
        ///@}
        
        /// @brief Costruisce l'iteratore end (stato sentinel).
        Iterator() = default;
        
        /**
         * @brief Costruttore interno usato da BitSet::begin().
         * @param p         Indice della word corrente nel backing store.
         * @param remaining Bit della word corrente non ancora visitati.
         * @param mask      Maschera con esattamente 1 bit acceso (prossimo da restituire).
         * @param nblock    Numero totale di word nel backing store.
         * @param store     Puntatore al backing store.
         */
        Iterator(std::uint64_t p, word_t remaining, word_t mask,
                 std::uint64_t nblock, const word_t* store)
        : storep_(p), remaining_(remaining), mask_(mask),
        nblock_(nblock), store_(store) {}
        
        /**
         * @brief Dereferenzia l'iteratore.
         * @return Indice (0-based) del bit correntemente puntato.
         */
        [[nodiscard]] reference operator*() const noexcept {
            return static_cast<uint64_t>(std::countr_zero(mask_)) +
            kBitsPerWord * storep_;
        }
        
        /**
         * @brief Pre-incremento: avanza al prossimo bit a 1.
         * @return Riferimento a *this dopo l'avanzamento.
         */
        Iterator& operator++() noexcept {
            if (remaining_) {
                mask_       = remaining_ & (~remaining_ + 1); // isola il LSB acceso
                remaining_ ^= mask_;                          // rimuove il LSB via XOR
            } else {
                AdvanceBlock();
            }
            return *this;
        }
        
        /**
         * @brief Post-incremento: restituisce copia prima dell'avanzamento.
         * @return Copia dell'iteratore prima dell'avanzamento.
         */
        Iterator operator++(int) noexcept {
            auto tmp = *this;
            ++(*this);
            return tmp;
        }
        
        /// @brief Uguaglianza tra iteratori.
        [[nodiscard]] bool operator==(const Iterator& o) const noexcept {
            return storep_ == o.storep_ && mask_ == o.mask_;
        }
        
        /// @brief Disuguaglianza tra iteratori.
        [[nodiscard]] bool operator!=(const Iterator& o) const noexcept {
            return !(*this == o);
        }
        
    private:
        /**
         * @brief Avanza all'inizio della prossima word non-zero.
         * @details Salta le word interamente a 0 senza ramo per bit; al termine
         *          imposta mask_ e remaining_ sulla nuova word, oppure porta
         *          l'iteratore nello stato sentinel (end).
         */
        void AdvanceBlock() noexcept {
            ++storep_;
            while (storep_ < nblock_ && store_[storep_] == 0) ++storep_;
            if (storep_ < nblock_) {
                const word_t b = store_[storep_];
                mask_          = b & (~b + 1);
                remaining_     = b ^ mask_;
            } else {
                *this = Iterator{}; // stato end
            }
        }
        
        std::uint64_t storep_    = 0;       ///< Indice della word corrente nel backing store.
        word_t        remaining_ = 0;       ///< Bit della word corrente non ancora visitati.
        word_t        mask_      = 0;       ///< Maschera con esattamente 1 bit acceso.
        std::uint64_t nblock_    = 0;       ///< Numero totale di word nel backing store.
        const word_t* store_     = nullptr; ///< Puntatore (non owning) al backing store.
    };
    
    /// @brief Alias STL per l'iteratore (i bit sono restituiti per valore).
    using iterator       = Iterator;
    /// @brief Alias STL per l'iteratore costante (identico a iterator).
    using const_iterator = Iterator;
    
    /**
     * @brief Restituisce un iteratore al primo bit a 1.
     * @return Iterator che punta al primo bit impostato, oppure end() se vuoto.
     */
    [[nodiscard]] iterator begin() const noexcept {
        std::uint64_t p = 0;
        while (p < nblock_ && store_[p] == 0) ++p;
        if (p >= nblock_) return {};
        const word_t b    = store_[p];
        const word_t mask = b & (~b + 1);
        return { p, b ^ mask, mask, nblock_, store_.get() };
    }
    
    /// @brief Restituisce il sentinel di fine iterazione.
    [[nodiscard]] iterator end() const noexcept { return {}; }
    
    /// @brief Equivalente a begin() (interfaccia STL const_iterator).
    [[nodiscard]] const_iterator cbegin() const noexcept { return begin(); }
    
    /// @brief Equivalente a end() (interfaccia STL const_iterator).
    [[nodiscard]] const_iterator cend() const noexcept { return end(); }
    
    // =========================================================================
    // Utilità
    // =========================================================================

    /**
     * @brief Imposta a 1 tutti i bit negli indici [0, @p n).
     * @param n Numero di bit da impostare (deve essere ≤ capacity()).
     */
    void FillUpTo(uint64_t n) noexcept {
        for (uint64_t i = 0; i < n; ++i) {
            set(i);
        }
    }
    
    /**
     * @brief Trova l'indice del primo bit che differisce tra @p a e @p b.
     * @details Calcola XOR word-by-word e usa @c countr_zero sul primo risultato non zero.
     * @param a Primo BitSet.
     * @param b Secondo BitSet.
     * @return Indice (0-based) del primo bit diverso, oppure
     *         @c static_cast<uint64_t>(-1) se i due insiemi sono identici
     *         nei blocchi confrontati.
     */
    static uint64_t FindFirstDifference(const BitSet& a,
                                        const BitSet& b) noexcept {
        const uint64_t min_blocks = std::min(a.nblock_, b.nblock_);
        for (uint64_t i = 0; i < min_blocks; ++i) {
            const word_t diff = a.store_[i] ^ b.store_[i];
            if (diff != 0) {
                const uint64_t bit_pos =
                static_cast<uint64_t>(std::countr_zero(diff));
                return (i * kBitsPerWord) + bit_pos;
            }
        }
        return static_cast<uint64_t>(-1);
    }
    
    /**
     * @brief Restituisce la capacità del BitSet — STL-compatible.
     * @return Numero di bit logici allocati.
     */
    [[nodiscard]] uint64_t capacity() const noexcept { return capacity_; }
    
    /**
     * @brief Converte il BitSet in un vettore degli indici dei bit a 1.
     * @details Mantenuto per retrocompatibilità; preferire l'iteratore nei
     *          cicli range-for per evitare l'allocazione del vettore.
     * @return @c std::vector<uint64_t> con gli indici (0-based) dei bit impostati.
     */
    [[nodiscard]] std::vector<uint64_t> ToVector() const {
        std::vector<uint64_t> res;
        for (uint64_t i = 0; i < nblock_; ++i) {
            word_t block = store_[i];
            while (block != 0) {
                const uint64_t bit_pos =
                static_cast<uint64_t>(std::countr_zero(block));
                res.push_back((i * kBitsPerWord) + bit_pos);
                block &= block - 1; // azzera il LSB acceso
            }
        }
        return res;
    }
    
    /**
     * @brief Verifica se il BitSet è vuoto — STL-compatible.
     * @return @c true se tutti i bit sono a 0, @c false altrimenti.
     */
    [[nodiscard]] inline bool empty() const noexcept {
        for (uint64_t i = 0; i < nblock_; ++i) {
            if (store_[i] != 0) return false;
        }
        return true;
    }
    
    /**
     * @brief Calcola la differenza insiemistica e salva il risultato in @c this.
     * @details Imposta @c this = @p a \ @p b, ovvero @c a AND NOT b, word-by-word.
     *          I blocchi di @c this eccedenti vengono azzerati.
     * @param a Minuendo.
     * @param b Sottraendo.
     */
    inline void SetDifference(const BitSet& a, const BitSet& b) noexcept {
        const uint64_t min_blocks = std::min({nblock_, a.nblock_, b.nblock_});
        for (uint64_t i = 0; i < min_blocks; ++i) {
            store_[i] = a.store_[i] & ~b.store_[i];
        }
        for (uint64_t i = min_blocks; i < nblock_; ++i) {
            store_[i] = 0;
        }
    }
    
    /**
     * @brief Differenza insiemistica in-place: @c this = this \ @p other.
     * @details Rimuove da @c this tutti i bit presenti in @p other: @c this AND NOT other.
     * @param other Insieme da sottrarre.
     */
    inline void DifferenceInPlace(const BitSet& other) noexcept {
        const uint64_t min_blocks = std::min(nblock_, other.nblock_);
        for (uint64_t i = 0; i < min_blocks; ++i) {
            store_[i] &= ~other.store_[i];
        }
    }
};
