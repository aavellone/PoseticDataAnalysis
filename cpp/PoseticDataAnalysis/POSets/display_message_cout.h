/**
 * @file display_message_cout.h
 * @brief Reporter di progresso basato su std::cout (SOLO per uso standalone).
 *
 * @warning Questo header NON deve essere incluso da alcuna translation unit
 * compilata nel package R: scrive su stdout tramite std::cout, cosa vietata
 * dalla policy CRAN ("Compiled code should not ... write to stdout/stderr
 * instead of to the console"). Nel package usare i reporter che passano da
 * Rprintf (vedi r_display.h: DisplayMessageEvaluationExactR, ...). Questo file
 * resta disponibile per eseguibili di benchmark/standalone che non linkano R.
 *
 * @par Requirements
 *   C++20 or later (`-std=c++20`)
 */

#pragma once

#include <atomic>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <optional>
#include <string>

#include "display_message.h"

// =========================================================================
// DisplayMessageDimensionalityReductionCout
// =========================================================================

/**
 * @class DisplayMessageDimensionalityReductionCout
 * @brief Reporter di progresso per la riduzione dimensionale su @c std::cout.
 *
 * @details Stampa "X out of Y (Z%) linear extensions used." al massimo ogni
 * @c output_every_sec_ secondi, e sempre alla chiamata di Stop().
 *
 * Il contatore @c done_ è un riferimento a un @c std::atomic per evitare
 * data race durante la lettura del progresso da thread multipli.
 *
 * @par Thread safety
 *   - Display() e Stop() sono thread-safe (protetti da @c mutex_).
 *   - Display(const std::string&) è protetto da @c mutex_.
 */
class DisplayMessageDimensionalityReductionCout final : public DisplayMessage {
private:
    /// @brief Clock ad alta risoluzione usato per i timestamp.
    using Clock     = std::chrono::high_resolution_clock;
    /// @brief Tipo del time point associato a @ref Clock.
    using TimePoint = std::chrono::time_point<Clock>;

    const std::atomic<std::uint64_t>& done_;  ///< Riferimento al contatore di avanzamento (atomic).
    TimePoint                         start_; ///< Timestamp di inizio computazione.
    TimePoint                         last_;  ///< Timestamp dell'ultimo output prodotto.

public:
    /**
     * @brief Costruisce il reporter collegandolo al contatore atomico di progresso.
     * @param done             Riferimento al contatore atomico degli elementi completati.
     * @param total            Numero totale di elementi da elaborare.
     * @param output_every_sec Intervallo minimo in secondi tra output consecutivi (default: 10).
     */
    DisplayMessageDimensionalityReductionCout(
                                              const std::atomic<std::uint64_t>& done,
                                              std::optional<std::uint64_t> total,
                                              std::uint64_t output_every_sec)
    : DisplayMessage(total, output_every_sec),
    done_(done){}

    /**
     * @brief Inizializza i timestamp di inizio e dell'ultimo output.
     * @details Thread-safe: acquisisce @c mutex_ prima di scrivere i timestamp.
     */
    void Start() override {
        std::lock_guard lock(mutex_);
        start_ = last_ = Clock::now();
    }

    /**
     * @brief Forza la stampa dello stato finale al termine della computazione.
     * @details Chiama internamente Print() che acquisisce @c mutex_.
     */
    void Stop() override {
        Print();
    }

    /**
     * @brief Stampa il progresso corrente se l'intervallo configurato è trascorso.
     * @details Legge il tempo corrente prima di acquisire @c mutex_, poi
     *          confronta l'elapsed con @c output_every_sec_. La prima chiamata
     *          dopo Start() produce sempre output (last_ == start_).
     *          Thread-safe.
     */
    void Display() override {
        const auto now = Clock::now();
        std::lock_guard lock(mutex_);
        const auto elapsed =
        std::chrono::duration_cast<std::chrono::seconds>(now - last_).count();
        if (last_ == start_ || static_cast<std::uint64_t>(elapsed) >= output_every_sec_) {
            const std::uint64_t current_done = done_.load(std::memory_order_relaxed);
            if (total_.has_value()) {
                const std::uint64_t pct = total_.value() > 0 ? static_cast<std::uint64_t>(
                                                                                          (static_cast<double>(current_done) / total_.value()) * 100.0)
                : 0;
                std::cout << current_done << " out of " << total_.value()
                << " (" << pct << "%) linear extensions used.\n";
            } else {
                std::cout << current_done << " linear extensions used.\n";
            }
            last_ = now;
        }
    }

    /**
     * @brief Stampa immediatamente un messaggio arbitrario su @c std::cout.
     * @param msg Messaggio da visualizzare.
     * @details Thread-safe: acquisisce @c mutex_ per evitare output interleaved.
     */
    void Display(const std::string& msg) override {
        std::lock_guard lock(mutex_);
        std::cout << msg << '\n';
    }

private:
    /**
     * @brief Stampa lo stato corrente incondizionatamente.
     * @details Utilizzato da Stop() per garantire che l'ultimo stato venga sempre
     *          visualizzato indipendentemente dall'intervallo configurato.
     *          Acquisisce @c mutex_ internamente. Thread-safe.
     */
    void Print() {
        const auto now = Clock::now();
        std::lock_guard lock(mutex_);
        const std::uint64_t current_done = done_.load(std::memory_order_relaxed);
        if (total_.has_value()) {
            const std::uint64_t pct = total_.value() > 0 ? static_cast<std::uint64_t>(
                                                                                      (static_cast<double>(current_done) / total_.value()) * 100.0)
            : 0;
            std::cout << current_done << " out of " << total_.value() << " (" << pct << "%) linear extensions used.\n";
        } else {
            std::cout << current_done << "  linear extensions used.\n";
        }
        last_ = now;
    }
};
