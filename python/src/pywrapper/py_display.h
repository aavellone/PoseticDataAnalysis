/**
 * @file py_display.h
 * @brief Reporter di progresso per Python (analogo di r_display.h).
 *
 * @details Le computazioni lunghe (enumerazione esatta, campionamento MCMC)
 * riportano il progresso sul terminale. In R si usa Rprintf; qui si usa
 * PySys_WriteStdout, che si integra con sys.stdout di Python.
 *
 * IMPORTANTE: PySys_WriteStdout richiede che il GIL sia trattenuto. Tutti gli
 * entry-point che usano questi reporter chiamano il core in modo sincrono sul
 * thread chiamante SENZA rilasciare il GIL, quindi la condizione e' rispettata.
 * (Se in futuro si rilascera' il GIL attorno al calcolo, il reporter dovra'
 * riacquisirlo prima di stampare.)
 */

#pragma once

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <format>
#include <mutex>
#include <optional>
#include <string>
#include <utility>

#include "display_message.h"

namespace PyDisplay {

/// Stampa una riga su sys.stdout di Python (il '%s' rende sicuro il testo utente).
inline void Line(const std::string& msg) {
    PySys_WriteStdout("%s\n", msg.c_str());
}

/// Stampa un avviso su sys.stdout (analogo soft di forward_warning_to_r).
inline void Warning(const std::string& msg) {
    PySys_WriteStdout("Warning: %s\n", msg.c_str());
}

}  // namespace PyDisplay

/**
 * @class DisplayMessageStdout
 * @brief Reporter di progresso a intervallo su stdout Python.
 *
 * @details Copre sia il conteggio delle estensioni "valutate" sia quelle
 * "generate", scegliendo il verbo al costruttore. Se @c total e' presente
 * mostra anche la percentuale. Speculare a DisplayMessageEvaluationBDR /
 * DisplayMessageEvaluationExactR / DisplayMessageLEGetR di r_display.h.
 */
class DisplayMessageStdout final : public DisplayMessage {
  public:
    DisplayMessageStdout(std::uint64_t& count, std::optional<std::uint64_t> total,
                         std::uint64_t output_every_sec, std::string verb = "evaluated")
        : DisplayMessage(total, output_every_sec),
          count_(count),
          verb_(std::move(verb)),
          start_time_(Clock::now()),
          last_time_(start_time_),
          started_(false) {}

    void Start() override {
        start_time_ = Clock::now();
        last_time_ = start_time_;
        started_ = true;
    }

    void Stop() override {
        std::lock_guard<std::mutex> lock(mutex_);
        PrintLocked();
        started_ = false;
    }

    void Display() override {
        if (!started_) [[unlikely]] {
            return;
        }
        const auto elapsed = Clock::now() - last_time_;
        const auto elapsed_sec = static_cast<std::uint64_t>(
            std::chrono::duration_cast<std::chrono::seconds>(elapsed).count());
        if (elapsed_sec > output_every_sec_) [[unlikely]] {
            std::lock_guard<std::mutex> lock(mutex_);
            PrintLocked();
        }
    }

    void Display(const std::string& msg) override {
        std::lock_guard<std::mutex> lock(mutex_);
        PyDisplay::Line(msg);
    }

  private:
    using Clock = std::chrono::high_resolution_clock;
    using TimePoint = std::chrono::time_point<Clock>;

    std::uint64_t& count_;
    std::string verb_;
    TimePoint start_time_;
    TimePoint last_time_;
    bool started_;

    // Pre-condizione: il chiamante detiene mutex_.
    void PrintLocked() {
        std::string msg;
        if (total_.has_value() && total_.value() > 0) {
            const std::uint64_t pct = static_cast<std::uint64_t>(
                (static_cast<double>(count_) / static_cast<double>(total_.value())) * 100.0);
            msg = std::format("{} out of {} ({}%) linear extensions {}.", count_,
                              total_.value(), pct, verb_);
        } else {
            msg = std::format("{} linear extensions {}.", count_, verb_);
        }
        PyDisplay::Line(msg);
        last_time_ = Clock::now();
    }
};

/**
 * @class DisplayMessageAtomicStdout
 * @brief Reporter thread-safe SENZA GIL per computazioni multi-thread (dim-reduction).
 *
 * @details Legge un contatore std::atomic e stampa via std::fprintf(stdout) —
 * thread-safe a livello POSIX e indipendente dal GIL, cosi' puo' essere invocato
 * da qualunque worker mentre il chiamante ha rilasciato il GIL. Speculare a
 * DisplayMessageEvaluationBDRThreads di r_display.h.
 */
class DisplayMessageAtomicStdout final : public DisplayMessage {
  public:
    DisplayMessageAtomicStdout(const std::atomic<std::uint64_t>& count,
                               std::optional<std::uint64_t> total,
                               std::uint64_t output_every_sec)
        : DisplayMessage(total, output_every_sec),
          count_(count),
          start_time_(Clock::now()),
          last_time_(start_time_),
          started_(false) {}

    void Start() override {
        start_time_ = Clock::now();
        last_time_ = start_time_;
        started_ = true;
    }

    void Stop() override {
        std::lock_guard<std::mutex> lock(mutex_);
        PrintLocked();
        started_ = false;
    }

    void Display() override {
        if (!started_) [[unlikely]] {
            return;
        }
        const auto elapsed = Clock::now() - last_time_;
        const auto elapsed_sec = static_cast<std::uint64_t>(
            std::chrono::duration_cast<std::chrono::seconds>(elapsed).count());
        if (elapsed_sec > output_every_sec_) [[unlikely]] {
            std::lock_guard<std::mutex> lock(mutex_);
            PrintLocked();
        }
    }

    void Display(const std::string& msg) override {
        std::lock_guard<std::mutex> lock(mutex_);
        std::fprintf(stdout, "%s\n", msg.c_str());
        std::fflush(stdout);
    }

  private:
    using Clock = std::chrono::high_resolution_clock;
    using TimePoint = std::chrono::time_point<Clock>;

    const std::atomic<std::uint64_t>& count_;
    TimePoint start_time_;
    TimePoint last_time_;
    bool started_;

    void PrintLocked() {
        const std::uint64_t current = count_.load(std::memory_order_relaxed);
        const std::uint64_t total = total_.has_value() ? total_.value() : 0;
        const std::uint64_t pct = total > 0
            ? static_cast<std::uint64_t>((static_cast<double>(current) / static_cast<double>(total)) * 100.0)
            : 0;
        std::fprintf(stdout, "%llu out of %llu (%llu%%) linear extensions used.\n",
                     static_cast<unsigned long long>(current),
                     static_cast<unsigned long long>(total),
                     static_cast<unsigned long long>(pct));
        std::fflush(stdout);
        last_time_ = Clock::now();
    }
};
