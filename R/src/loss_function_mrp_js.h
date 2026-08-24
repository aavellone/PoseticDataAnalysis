/**
 * @file loss_function_mrp_js.h
 * @brief Loss functions basate sulla Divergenza di Jensen-Shannon per matrici MRP.
 */

#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <vector>
#include <unordered_map>
#include <cstdint>

#include "loss_function_mrp.h"

/**
 * @class JS0LossFunctionMRP
 * @brief Funzione di perdita MRP non pesata basata su Jensen-Shannon.
 */
class JS0LossFunctionMRP final : public LossFunctionMRP {
private:
    double inv_log_base_;
    
public:
    JS0LossFunctionMRP(const std::vector<Tensor<double, 2>>& p,
                       const std::unordered_map<std::uint64_t, double>& w,
                       std::uint64_t base = 2);
    double operator()(const Tensor<double, 2>& x) override;
};

/**
 * @class JS1LossFunctionMRP
 * @brief Funzione JS normalizzata rispetto al complemento dell'entropia locale.
 */
class JS1LossFunctionMRP final : public LossFunctionMRP {
private:
    double inv_log_base_;
    
public:
    JS1LossFunctionMRP(const std::vector<Tensor<double, 2>>& p,
                       const std::unordered_map<std::uint64_t, double>& w,
                       std::uint64_t base = 2);
    double operator()(const Tensor<double, 2>& x) override;
};

/**
 * @class JS2LossFunctionMRP
 * @brief Funzione JS normalizzata rispetto all'entropia media globale della matrice.
 */
class JS2LossFunctionMRP final : public LossFunctionMRP {
private:
    double inv_log_base_;
    
public:
    JS2LossFunctionMRP(const std::vector<Tensor<double, 2>>& p,
                       const std::unordered_map<std::uint64_t, double>& w,
                       std::uint64_t base = 2);
    double operator()(const Tensor<double, 2>& x) override;
};
