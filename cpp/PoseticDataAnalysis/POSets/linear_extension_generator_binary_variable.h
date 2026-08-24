/**
 * @file linear_extension_generator_binary_variable.h
 * @brief Header file for the Binary Variable linear extension generator.
 */

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "linear_extension.h"
#include "linear_extension_generator.h"

/**
 * @class LEGBinaryVariable
 * @brief Generator for binary variable profiles (e.g., boolean hypercubes).
 */
class LEGBinaryVariable final : public LinearExtensionGenerator {
public:
    explicit LEGBinaryVariable(std::uint64_t num_variables);
    ~LEGBinaryVariable() override = default;
    
    void Start(std::uint64_t unused_id = 0) override;
    void Next() override;
    
    [[nodiscard]] bool HasNext() noexcept override {
        return current_number_le_ < max_number_le_;
    }
    
    [[nodiscard]] std::uint64_t NumberOfLe() const noexcept override {
        return max_number_le_;
    }
    
private:
    [[nodiscard]] static std::uint64_t CalculateNumProfiles(std::uint64_t n_var) noexcept;
    [[nodiscard]] static std::uint64_t Factorial(std::uint64_t n) noexcept;
    void BuildLeFromPermutation() noexcept;
    
    // L'ordine di dichiarazione ORA corrisponde strettamente alla Initializer List
    LinearExtension le_;
    LinearExtension le_inv_;
    std::vector<std::uint64_t> permutation_;
    std::vector<std::uint64_t> permutation_inv_;
    std::vector<std::uint64_t> internal_permutation_;
    
    std::uint64_t num_variables_;
    std::uint64_t num_profiles_;
    std::uint64_t perm_var_1_;
    std::uint64_t perm_var_2_;
    std::uint64_t max_number_le_;
};
