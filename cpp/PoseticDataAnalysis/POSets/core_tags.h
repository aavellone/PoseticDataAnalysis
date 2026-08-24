/**
 * @file memory_tags.h
 * @brief Global tags and policies for memory management in HPC contexts.
 */

#pragma once

/**
 * @struct Uninitialized
 * @brief Tag type used to disambiguate constructors for uninitialized allocation.
 */
struct Uninitialized {
    explicit Uninitialized() = default;
};

/**
 * @brief Global constant tag used to request uninitialized memory allocation.
 * @note Marked as inline (C++17) to prevent Multiple Definition / ODR violations
 * across multiple translation units.
 */
inline constexpr Uninitialized kUninitialized{};
