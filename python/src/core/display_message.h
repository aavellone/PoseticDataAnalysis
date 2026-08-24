/**
 * @file display_message.h
 * @brief Astrazioni per il reporting del progresso di computazioni lunghe.
 *
 * @details Fornisce:
 *  - @ref DisplayMessage                            — interfaccia puramente virtuale.
 *  - @ref DisplayMessageNull                        — implementazione no-op.
 *  (il reporter su stdout DisplayMessageDimensionalityReductionCout vive in
 *   display_message_cout.h, fuori dal package R per la policy CRAN).
 *
 * @par Requirements
 *   C++20 or later (`-std=c++20`)
 *
 * @par Naming conventions (Google C++ Style Guide)
 *   - Types / classes:              PascalCase
 *   - Methods:                      PascalCase, eccetto quelli che replicano
 *                                   l'interfaccia STL (to_string(), …) che restano snake_case
 *   - Data members:                 snake_case_  (trailing underscore; NO prefisso _)
 *   - Local variables / parameters: snake_case
 *   - Constants (constexpr):        kCamelCase
 */

#pragma once

#include <atomic>
#include <chrono>
#include <cstdint>
#include <format>
#include <mutex>
#include <optional>
#include <string>

// =========================================================================
// DisplayMessage — interfaccia astratta
// =========================================================================

/**
 * @class DisplayMessage
 * @brief Interfaccia puramente virtuale per i reporter di progresso.
 *
 * @details Le sottoclassi concrete decidono *come* visualizzare il progresso
 * (stdout, file, GUI, …). La classe base possiede solo l'intervallo di output
 * e un mutex per la thread-safety.
 *
 * @note La classe non è copiabile perché possiede un @c std::mutex.
 */
class DisplayMessage {
protected:
    mutable std::mutex mutex_;          ///< Protegge tutte le operazioni di output.
    std::uint64_t      output_every_sec_; ///< Intervallo minimo in secondi tra due output.
    std::optional<std::uint64_t> total_; ///< Numero totale di elementi da elaborare.
public:
    /**
     * @brief Costruisce il reporter con l'intervallo di output desiderato.
     * @param output_every_sec Secondi minimi tra output consecutivi (default: 10).
     */
    explicit DisplayMessage(std::optional<std::uint64_t> total, std::uint64_t output_every_sec = 10)
    : output_every_sec_(output_every_sec), total_(total) {}
    
    /// @brief Distruttore virtuale.
    virtual ~DisplayMessage() = default;
    
    /// @brief Non copiabile (possiede un mutex).
    DisplayMessage(const DisplayMessage&)            = delete;
    /// @brief Non copy-assignable (possiede un mutex).
    DisplayMessage& operator=(const DisplayMessage&) = delete;
    
    /**
     * @brief Chiamata una volta prima dell'inizio della computazione.
     * @details Le sottoclassi devono inizializzare qui i timestamp o le risorse.
     */
    virtual void Start() = 0;
    
    /**
     * @brief Chiamata una volta al termine della computazione.
     * @details Le sottoclassi devono garantire che l'ultimo stato venga stampato.
     */
    virtual void Stop() = 0;
    
    /**
     * @brief Mostra il progresso corrente se l'intervallo configurato è trascorso.
     * @details Può essere chiamata ad alta frequenza senza overhead: l'implementazione
     *          decide internamente se è il momento di produrre output.
     */
    virtual void Display() = 0;
    
    /**
     * @brief Mostra immediatamente un messaggio arbitrario.
     * @param msg Messaggio da visualizzare.
     */
    virtual void Display(const std::string& msg) = 0;
    
    /**
     * @brief Restituisce l'intervallo di output configurato.
     * @return Secondi minimi tra output consecutivi.
     */
    [[nodiscard]] std::uint64_t OutputSeconds() const noexcept {
        return output_every_sec_;
    }
    
    /**
     * @brief Restituisce una stringa descrittiva dello stato corrente (opzionale).
     * @details L'implementazione di default restituisce una stringa vuota.
     *          Le sottoclassi possono fare override per fornire un riepilogo.
     * @return Stringa human-readable dello stato, oppure stringa vuota.
     */
    [[nodiscard]] virtual std::string to_string() const { return {}; }
    
    void UpdateTotal(std::uint64_t t) noexcept {
        total_ = t;
    }
};

// =========================================================================
// DisplayMessageNull — implementazione no-op
// =========================================================================

/**
 * @class DisplayMessageNull
 * @brief Implementazione di @ref DisplayMessage che non produce alcun output.
 *
 * @details Utile come reporter disabilitato o come valore di default.
 * Tutti i metodi virtuali sono no-op: il compilatore può inlinarli ed
 * eliminarli completamente in fase di ottimizzazione.
 */
class DisplayMessageNull final : public DisplayMessage {
public:
    DisplayMessageNull() : DisplayMessage(std::optional<std::uint64_t>(), 0) {}
    
    /// @copydoc DisplayMessage::Start()
    void Start()                     override {}
    
    /// @copydoc DisplayMessage::Stop()
    void Stop()                      override {}
    
    /// @copydoc DisplayMessage::Display()
    void Display()                   override {}
    
    /// @copydoc DisplayMessage::Display(const std::string&)
    void Display(const std::string&) override {}
};

// =========================================================================
// DisplayMessageAtomicTick — contatore atomico per worker thread
// =========================================================================

/**
 * @class DisplayMessageAtomicTick
 * @brief Reporter "muto" per i worker della valutazione parallela multi-catena.
 *
 * @details Ogni chiamata a Display() incrementa un contatore atomico condiviso
 * (una unità = una estensione lineare elaborata). Non produce MAI output:
 * la stampa del progresso spetta al thread principale, unico autorizzato a
 * chiamare l'API di R (Rprintf). Start(), Stop() e Display(msg) sono no-op.
 */
class DisplayMessageAtomicTick final : public DisplayMessage {
private:
    std::atomic<std::uint64_t>& ticks_; ///< Contatore condiviso delle LE elaborate.

public:
    explicit DisplayMessageAtomicTick(std::atomic<std::uint64_t>& ticks)
    : DisplayMessage(std::optional<std::uint64_t>(), 0), ticks_(ticks) {}

    void Start()                     override {}
    void Stop()                      override {}
    void Display()                   override {
        ticks_.fetch_add(1, std::memory_order_relaxed);
    }
    void Display(const std::string&) override {}
};

// Il reporter di progresso su std::cout
// (DisplayMessageDimensionalityReductionCout) è stato spostato in
// display_message_cout.h: scrive su stdout, vietato dalla policy CRAN, e non
// deve finire nel .so del package. Nel package si usano i reporter basati su
// Rprintf definiti in r_display.h.
