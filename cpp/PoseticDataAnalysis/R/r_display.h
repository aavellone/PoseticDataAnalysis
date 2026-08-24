/**
 * @file r_display.h
 * @brief R-specific progress reporter implementations.
 */
#pragma once

#include <R_ext/Print.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <ctime>
#include <format>
#include <limits>
#include <mutex>
#include <string>

#include "display_message.h"

// ================================================================
// DisplayMessageEvaluationExactR
// ================================================================

class DisplayMessageEvaluationExactR final : public DisplayMessage {
public:
    /**
     * @brief Constructs the reporter.
     * @param count            Reference to the current extension count.
     * @param output_every_sec Minimum seconds between outputs (default: 10).
     */
    DisplayMessageEvaluationExactR(std::uint64_t& count, std::uint64_t output_every_sec)
        : DisplayMessage(std::optional<std::uint64_t>(), output_every_sec),
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
        std::lock_guard<std::mutex> lock(mutex_); // Uso del mutex_ dalla classe base
        PrintLocked();
        started_ = false;
    }
    
    void Display() override {
        if (!started_) [[unlikely]] return;
        
        const auto elapsed = Clock::now() - last_time_;
        const auto elapsed_sec = static_cast<std::uint64_t>(
                                                            std::chrono::duration_cast<std::chrono::seconds>(elapsed).count());
        
        // Uso di output_every_sec_ dalla classe base
        if (elapsed_sec > output_every_sec_) [[unlikely]] {
            std::lock_guard<std::mutex> lock(mutex_);
            PrintLocked();
        }
    }
    
    void Display(const std::string& msg) override {
        std::lock_guard<std::mutex> lock(mutex_);
        Rprintf("%s\n", msg.c_str());
    }
    
private:
    using Clock = std::chrono::high_resolution_clock;
    using TimePoint = std::chrono::time_point<Clock>;
    
    std::uint64_t& count_;
    TimePoint start_time_;
    TimePoint last_time_;
    bool started_;
    
    // Pre-condition: mutex_ must be held by the caller.
    void PrintLocked() {
        const auto now = Clock::now();
        
        const auto total_elapsed_sec = static_cast<std::uint64_t>(
                                                                  std::chrono::duration_cast<std::chrono::seconds>(now - start_time_).count()
                                                                  );
        
        // Calcolo di giorni, ore, minuti e secondi
        const std::uint64_t days = total_elapsed_sec / 86400;
        const std::uint64_t hours = (total_elapsed_sec % 86400) / 3600;
        const std::uint64_t minutes = (total_elapsed_sec % 3600) / 60;
        const std::uint64_t seconds = total_elapsed_sec % 60;
        
        std::string time_str;
        if (days > 0) {
            time_str = std::format("{}d {:02}h {:02}m {:02}s", days, hours, minutes, seconds);
        } else if (hours > 0) {
            time_str = std::format("{}h {:02}m {:02}s", hours, minutes, seconds);
        } else if (minutes > 0) {
            time_str = std::format("{}m {:02}s", minutes, seconds);
        } else {
            time_str = std::format("{}s", seconds);
        }
        
        const std::string msg = std::format("{} linear extensions evaluated in {}.",
                                            count_, time_str);
        
        Rprintf("%s\n", msg.c_str());
        
        last_time_ = now;
    }
};

// ================================================================
// DisplayMessageEvaluationBDRThreads
// ================================================================

class DisplayMessageEvaluationBDRThreads final : public DisplayMessage {
public:
    /**
     * @brief Constructs the reporter.
     * @param count            Reference to the atomic current extension count.
     * @param total            Optional total number of extensions.
     * @param output_every_sec Minimum seconds between outputs (default: 10).
     */
    DisplayMessageEvaluationBDRThreads(const std::atomic<std::uint64_t>& count,
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
        std::lock_guard<std::mutex> lock(mutex_); // Uso del mutex_ dalla classe base
        PrintLocked();
        started_ = false;
    }
    
    void Display() override {
        if (!started_) [[unlikely]] return;
        
        const auto elapsed = Clock::now() - last_time_;
        const auto elapsed_sec = static_cast<std::uint64_t>(
                                                            std::chrono::duration_cast<std::chrono::seconds>(elapsed).count());
        
        // Uso di output_every_sec_ dalla classe base
        if (elapsed_sec > output_every_sec_) [[unlikely]] {
            std::lock_guard<std::mutex> lock(mutex_);
            PrintLocked();
        }
    }
    
    void Display(const std::string& msg) override {
        std::lock_guard<std::mutex> lock(mutex_);
        Rprintf("%s\n", msg.c_str());
    }
    
private:
    using Clock = std::chrono::high_resolution_clock;
    using TimePoint = std::chrono::time_point<Clock>;
    
    // 1. Modificato per accettare una reference atomica in sola lettura
    const std::atomic<std::uint64_t>& count_;
    
    TimePoint start_time_;
    TimePoint last_time_;
    bool started_;
    
    // Pre-condition: mutex_ must be held by the caller.
    void PrintLocked() {
        // 2. Lettura ultra-veloce e thread-safe
        const std::uint64_t current_count = count_.load(std::memory_order_relaxed);
        
        const std::uint64_t pct = total_.value() > 0 ?
        static_cast<std::uint64_t>((static_cast<double>(current_count) /
                                    static_cast<double>(total_.value())) * 100.0) : 0;
        
        const std::string msg = std::format("{} out of {} ({}%) linear extensions evaluated.",
                                            current_count, total_.value(), pct);
        Rprintf("%s\n", msg.c_str());
        last_time_ = Clock::now();
    }
};

// ================================================================
// DisplayMessageEvaluationBDR
// ================================================================

class DisplayMessageEvaluationBDR final : public DisplayMessage {
public:
    /**
     * @brief Constructs the reporter.
     * @param count            Reference to the atomic current extension count.
     * @param total            Optional total number of extensions.
     * @param output_every_sec Minimum seconds between outputs (default: 10).
     */
    DisplayMessageEvaluationBDR(std::uint64_t& count,
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
        std::lock_guard<std::mutex> lock(mutex_); // Uso del mutex_ dalla classe base
        PrintLocked();
        started_ = false;
    }
    
    void Display() override {
        if (!started_) [[unlikely]] return;
        
        const auto elapsed = Clock::now() - last_time_;
        const auto elapsed_sec = static_cast<std::uint64_t>(
                                                            std::chrono::duration_cast<std::chrono::seconds>(elapsed).count());
        
        // Uso di output_every_sec_ dalla classe base
        if (elapsed_sec > output_every_sec_) [[unlikely]] {
            std::lock_guard<std::mutex> lock(mutex_);
            PrintLocked();
        }
    }
    
    void Display(const std::string& msg) override {
        std::lock_guard<std::mutex> lock(mutex_);
        Rprintf("%s\n", msg.c_str());
    }
    
private:
    using Clock = std::chrono::high_resolution_clock;
    using TimePoint = std::chrono::time_point<Clock>;
    
    // 1. Modificato per accettare una reference atomica in sola lettura
    const std::uint64_t& count_;
    
    TimePoint start_time_;
    TimePoint last_time_;
    bool started_;
    
    // Pre-condition: mutex_ must be held by the caller.
    void PrintLocked() {
        
        const std::uint64_t pct = total_.value() > 0 ?
        static_cast<std::uint64_t>((static_cast<double>(count_) /
                                    static_cast<double>(total_.value())) * 100.0) : 0;
        
        const std::string msg = std::format("{} out of {} ({}%) linear extensions evaluated.",
                                            count_, total_.value(), pct);
        Rprintf("%s\n", msg.c_str());
        last_time_ = Clock::now();
    }
};


// ================================================================
// DisplayMessageLEGetR
// ================================================================

class DisplayMessageLEGetR final : public DisplayMessage {
public:
    DisplayMessageLEGetR(std::uint64_t& count, std::optional<std::uint64_t> total, std::uint64_t output_every_sec)
    : DisplayMessage(total, output_every_sec),
    count_(count),
    start_tick_(std::clock()),
    start_time_(Clock::now()),
    last_time_(start_time_),
    show_expected_(true),
    started_(false) {}
    
    void Start() override {
        start_tick_ = std::clock();
        start_time_ = Clock::now();
        last_time_ = start_time_;
        show_expected_ = true;
        started_ = true;
    }
    
    void Stop() override {
        std::lock_guard<std::mutex> lock(mutex_);
        PrintLocked();
        started_ = false;
    }
    
    void Display() override;
    
    void Display(const std::string& msg) override {
        std::lock_guard<std::mutex> lock(mutex_);
        Rprintf("%s\n", msg.c_str());
    }
    
private:
    using Clock = std::chrono::high_resolution_clock;
    using TimePoint = std::chrono::time_point<Clock>;
    
    static constexpr std::uint64_t kExpectedTimeSec = 10;
    
    std::uint64_t& count_;
    std::clock_t start_tick_;
    TimePoint start_time_;
    TimePoint last_time_;
    bool show_expected_;
    bool started_;
    
    static std::string FormatDuration(std::uint64_t total_seconds) {
        const std::uint64_t hours = total_seconds / 3600;
        const std::uint64_t minutes = (total_seconds % 3600) / 60;
        const std::uint64_t seconds = total_seconds % 60;
        return std::format("{}h {}m {}s", hours, minutes, seconds);
    }
    
    // Pre-condition: mutex_ must be held by the caller.
    void PrintLocked() {
        std::string msg;
        if (!total_.has_value()) {
            msg = std::format("{} linear extensions generated.", count_);
        } else {
            const double pct = (static_cast<double>(count_) / static_cast<double>(total_.value())) * 100.0;
            msg = std::format("{} out of {} ({:.0f}%) linear extensions generated.",
                              count_, total_.value(), pct);
        }
        Rprintf("%s\n", msg.c_str());
        last_time_ = Clock::now();
    }
};

// ================================================================
// DisplayMessageDimensionalityReductionR
// ================================================================

class DisplayMessageDimensionalityReductionR final : public DisplayMessage {
public:
    DisplayMessageDimensionalityReductionR(std::uint64_t& done, std::optional<std::uint64_t> total, std::uint64_t output_every_sec)
    : DisplayMessage(total, output_every_sec),
    done_(done),
    started_(false) {}
    
    void Start() override {
        last_time_ = Clock::now();
        started_ = true;
    }
    
    void Stop() override {
        std::lock_guard<std::mutex> lock(mutex_);
        PrintLocked();
        started_ = false;
    }
    
    void Display() override {
        if (!started_) [[unlikely]] return;
        
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
        Rprintf("%s\n", msg.c_str());
    }
    
private:
    using Clock = std::chrono::high_resolution_clock;
    using TimePoint = std::chrono::time_point<Clock>;
    
    std::uint64_t& done_;
    TimePoint last_time_;
    bool started_;
    
    // Pre-condition: mutex_ must be held by the caller.
    void PrintLocked() {
        const std::uint64_t pct = total_ > 0 ? static_cast<std::uint64_t>((static_cast<double>(done_) /
                                                                           static_cast<double>(total_.value())) * 100.0)
        : 0;
        
        const std::string msg = std::format("{} out of {} ({}%) linear extensions used.",
                                            done_, total_.value(), pct);
        Rprintf("%s\n", msg.c_str());
        last_time_ = Clock::now();
    }
};
