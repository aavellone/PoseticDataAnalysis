/**
 * @file loss_function_mrp_l1.h
 * @brief Funzioni di perdita MRP basate su distanze lineari e Lower Bounds (L1).
 */

#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <vector>
#include <utility>
#include <cstdint>

#include "loss_function_mrp.h"

/**
 * @class LBMRP
 * @brief Funzione di Lower Bound per matrici MRP (media pesata L1 relativa).
 */
class LBMRP final : public LossFunctionMRP {
public:
    using LossFunctionMRP::LossFunctionMRP;
    double operator()(const Tensor<double, 2>& x) override;
};

/**
 * @class L1LossFunctionMRP
 * @brief Funzione di perdita basata sulla distanza di Manhattan (L1) strettamente superiore.
 * @cite Kemeny, J. G., & Snell, J. L. (1962). "Mathematical models in the social sciences".
 */
class L1LossFunctionMRP final : public LossFunctionMRP {
public:
    using LossFunctionMRP::LossFunctionMRP;
    double operator()(const Tensor<double, 2>& x) override;
};

/**
 * @class LBMRP2
 * @brief Estensione della LBMRP per operare su mappature multiple di indici (Pipeline V2).
 */
class LBMRP2 final : public LossFunctionMRPV2 {
private:
    /**
     * @struct ResolvedIdx
     * @brief Struttura piatta per la pre-risoluzione degli indici.
     */
    struct ResolvedIdx {
        std::uint64_t orig;
        std::uint64_t res;
        double w;
    };
    
    void Resolve(std::vector<ResolvedIdx>& out,
                 const std::vector<std::vector<std::uint_fast64_t>>& funcs) const;
    
public:
    using LossFunctionMRPV2::LossFunctionMRPV2;
    
    double operator()(const Tensor<double, 2>& x) override;
    
    double operator()(const Tensor<double, 2>& x,
                      const std::vector<std::vector<std::uint_fast64_t>>& rows_functions,
                      const std::vector<std::vector<std::uint_fast64_t>>& cols_functions) override;
    
    void operator()(const Tensor<double, 2>& x,
                    const std::vector<std::vector<std::uint_fast64_t>>& rows_functions,
                    const std::vector<std::vector<std::uint_fast64_t>>& cols_functions,
                    std::vector<std::pair<std::uint_fast64_t, double>>& result) override;
};
