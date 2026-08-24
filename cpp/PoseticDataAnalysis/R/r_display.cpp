/**
 * @file r_display.cpp
 * @brief R-specific progress reporter implementations.
 */

#include "r_display.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <ctime>
#include <format>
#include <limits>
#include <mutex>
#include <string>

#include <R_ext/Print.h>

void DisplayMessageLEGetR::Display() {
    if (!started_) [[unlikely]] return;
    
    const auto now = Clock::now();
    const auto wall_secs = static_cast<std::uint64_t>(
                                                      std::chrono::duration_cast<std::chrono::seconds>(now - start_time_).count());
    const auto elapsed_last = static_cast<std::uint64_t>(
                                                         std::chrono::duration_cast<std::chrono::seconds>(now - last_time_).count());
    
    // ---- Phase 1: expected-time block ----
    if (total_.has_value() && show_expected_ && wall_secs > kExpectedTimeSec) [[unlikely]] {
        const double cpu_secs = static_cast<double>(std::clock() - start_tick_) / CLOCKS_PER_SEC;
        const std::uint64_t safe_count = std::max(count_, std::uint64_t{1});
        
        // Cast total_ to double first to avoid integer overflow in large datasets
        const std::uint64_t eta_cpu = static_cast<std::uint64_t>(
                                                                 (static_cast<double>(total_.value()) / static_cast<double>(safe_count)) * cpu_secs);
        const std::uint64_t eta_wall = static_cast<std::uint64_t>(
                                                                  (static_cast<double>(total_.value()) / static_cast<double>(safe_count)) * wall_secs);
        
        const std::string msg_cpu = std::format("Expected cpu time: {}", FormatDuration(eta_cpu));
        const std::string msg_wall = std::format("Expected elapsed time: {}", FormatDuration(eta_wall));
        
        {
            // Uso diretto di mutex_ ereditato dalla base class
            std::lock_guard<std::mutex> lock(mutex_);
            Rprintf("%s\n%s\n", msg_cpu.c_str(), msg_wall.c_str());
        }
        show_expected_ = false;
    }
    
    // ---- Phase 2: periodic percentage output ----
    if (!show_expected_ && elapsed_last > output_every_sec_) [[unlikely]] {
        std::lock_guard<std::mutex> lock(mutex_);
        PrintLocked();
    }
}
