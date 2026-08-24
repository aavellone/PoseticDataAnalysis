
// MSVC: explicit standard includes (not pulled in transitively).
#include <cstddef>
/**
 * @file linear_extension_generator_binary_variable.cpp
 * @brief Implementation of the Binary Variable linear extension generator.
 */

#include "linear_extension_generator_binary_variable.h"

#include <algorithm>
#include <cstdint>
#include <string>

#include "my_exception.h"

std::uint64_t LEGBinaryVariable::CalculateNumProfiles(std::uint64_t n_var) noexcept {
    return (1ULL << n_var);
}

std::uint64_t LEGBinaryVariable::Factorial(std::uint64_t n) noexcept {
    std::uint64_t res = 1;
    for (std::uint64_t i = 2; i <= n; ++i) {
        res *= i;
    }
    return res;
}

// CORREZIONE 1: L'ordine della Initializer List ora è perfettamente legale e
// pre-alloca direttamente i vettori per evitare riassegnazioni
LEGBinaryVariable::LEGBinaryVariable(std::uint64_t num_variables)
: LinearExtensionGenerator(CalculateNumProfiles(num_variables)),
le_(CalculateNumProfiles(num_variables)),
le_inv_(CalculateNumProfiles(num_variables)),
permutation_(num_variables, 0),
permutation_inv_(num_variables, 0),
internal_permutation_(num_variables > 1 ? num_variables - 2 : 0, 0),
num_variables_(num_variables),
num_profiles_(CalculateNumProfiles(num_variables)),
perm_var_1_(0),
perm_var_2_(0),
max_number_le_(Factorial(num_variables)) {
    
    if (num_variables_ < 2) {
        throw MyException("LEGBinaryVariable requires at least 2 variables.");
    }
}

void LEGBinaryVariable::BuildLeFromPermutation() noexcept {
    const std::uint64_t offset = num_variables_ - 1;
    
    for (std::uint64_t n = 0; n < le_.size(); ++n) {
        std::uint64_t value_lex = 0;
        std::uint64_t value_lex_inv = 0;
        
        for (std::uint64_t i = 0; i < num_variables_; ++i) {
            std::uint64_t old_val = (n >> (offset - i)) & 1ULL;
            value_lex |= (old_val << (offset - permutation_[i]));
            value_lex_inv |= (old_val << (offset - permutation_inv_[i]));
        }
        le_.Set(n, value_lex);
        le_inv_.Set(n, value_lex_inv);
    }
}

void LEGBinaryVariable::Start(std::uint64_t /*unused_id*/) {
    perm_var_1_ = 0;
    perm_var_2_ = 1;
    
    for (std::uint64_t i = 0, j = 0; i < num_variables_; ++i) {
        if (i != perm_var_1_ && i != perm_var_2_) {
            internal_permutation_[j++] = i;
        }
    }
    
    permutation_[0] = perm_var_1_;
    permutation_inv_[0] = perm_var_2_;
    
    const std::size_t perm_size = internal_permutation_.size();
    for (std::uint64_t p = 0; p < perm_size; p++) {
        permutation_[p + 1] = internal_permutation_[p];
        permutation_inv_[p + 1] = internal_permutation_[perm_size - p - 1];
    }
    
    permutation_[num_variables_ - 1] = perm_var_2_;
    permutation_inv_[num_variables_ - 1] = perm_var_1_;
    
    BuildLeFromPermutation();
    
    for (std::uint64_t k = 0; k < current_linear_extension_.size(); ++k) {
        current_linear_extension_.Set(k, le_.GetVal(k));
    }
    
    started_ = true;
    current_number_le_ = 1;
}

void LEGBinaryVariable::Next() {
    if (!started_) [[unlikely]] {
        throw MyException("LEGBinaryVariable error: Start() not called.");
    }
    
    if (current_number_le_ % 2 == 1) {
        for (std::uint64_t k = 0; k < current_linear_extension_.size(); ++k) {
            current_linear_extension_.Set(k, le_inv_.GetVal(k));
        }
    } else {
        if (!std::next_permutation(internal_permutation_.begin(), internal_permutation_.end())) {
            if (perm_var_2_ >= num_variables_ - 1) {
                
                // CORREZIONE 2: Ora l'eccezione viene sollevata correttamente appena
                // perm_var_1_ raggiunge l'ultimo valore valido (n - 2)
                if (perm_var_1_ >= num_variables_ - 2) {
                    throw MyException("LEGBinaryVariable error: Next() exhausted permutations.");
                }
                
                ++perm_var_1_;
                perm_var_2_ = perm_var_1_;
            }
            
            ++perm_var_2_;
            
            for (std::uint64_t i = 0, j = 0; i < num_variables_; ++i) {
                if (i != perm_var_1_ && i != perm_var_2_) {
                    internal_permutation_[j++] = i;
                }
            }
        }
        
        permutation_[0] = perm_var_1_;
        permutation_inv_[0] = perm_var_2_;
        
        const std::size_t perm_size = internal_permutation_.size();
        for (std::uint64_t p = 0; p < perm_size; p++) {
            permutation_[p + 1] = internal_permutation_[p];
            permutation_inv_[p + 1] = internal_permutation_[perm_size - p - 1];
        }
        
        permutation_[num_variables_ - 1] = perm_var_2_;
        permutation_inv_[num_variables_ - 1] = perm_var_1_;
        
        BuildLeFromPermutation();
        
        for (std::uint64_t k = 0; k < current_linear_extension_.size(); ++k) {
            current_linear_extension_.Set(k, le_.GetVal(k));
        }
    }
    
    ++current_number_le_;
}
