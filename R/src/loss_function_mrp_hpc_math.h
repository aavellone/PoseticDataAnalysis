/**
 * @file hpc_math.h
 * @brief Funzioni matematiche scalari altamente ottimizzate per la teoria dell'informazione e le matrici MRP.
 *
 * @details Questo modulo evita allocazioni dinamiche (es. `std::vector` temporanei)
 * ed effettua calcoli diretti su distribuzioni binarie (p e 1-p), essenziali per
 * massimizzare il throughput delle pipeline HPC.
 *
 * @par Requisiti
 * C++20 o superiore.
 */

#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <vector>

#include <cmath>

namespace hpc_math {
    
    /**
     * @brief Calcola l'entropia di Shannon per una distribuzione a due stati.
     * @param p Probabilità del primo stato.
     * @param inv_log_base Inverso del logaritmo naturale della base desiderata (1 / ln(b)).
     * @return Il valore dell'entropia normalizzato rispetto alla base.
     *
     * @cite Shannon, C. E. (1948). "A Mathematical Theory of Communication". The Bell System Technical Journal.
     */
    [[nodiscard]] inline double Entropy2(double p, double inv_log_base) noexcept {
        if (p <= 1e-15 || p >= 1.0 - 1e-15) return 0.0;
        return -(p * std::log(p) + (1.0 - p) * std::log(1.0 - p)) * inv_log_base;
    }
    
    /**
     * @brief Calcola la divergenza di Kullback-Leibler per due distribuzioni binomiali.
     * @param p Probabilità target (reale).
     * @param q Probabilità approssimata (modello).
     * @param inv_log_base Inverso del logaritmo naturale della base.
     * @return La divergenza KL (Entropia Relativa).
     *
     * @cite Kullback, S., & Leibler, R. A. (1951). "On information and sufficiency". Annals of Mathematical Statistics.
     */
    [[nodiscard]] inline double KLDivergence2(double p, double q, double inv_log_base) noexcept {
        double kl = 0.0;
        if (p > 1e-15 && q > 1e-15) kl += p * std::log(p / q);
        if (p < 1.0 - 1e-15 && q < 1.0 - 1e-15) kl += (1.0 - p) * std::log((1.0 - p) / (1.0 - q));
        return kl * inv_log_base;
    }
    
    /**
     * @brief Calcola la metrica di divergenza di Jensen-Shannon per due distribuzioni binomiali.
     * @details La divergenza JS è una versione simmetrica, limitata e sempre finita della divergenza KL.
     * @param p Prima probabilità.
     * @param q Seconda probabilità.
     * @param inv_log_base Inverso del logaritmo naturale della base.
     * @return Il valore della divergenza JS.
     *
     * @cite Lin, J. (1991). "Divergence measures based on the Shannon entropy". IEEE Transactions on Information Theory.
     */
    [[nodiscard]] inline double JSDivergence2(double p, double q, double inv_log_base) noexcept {
        double m = 0.5 * (p + q);
        return 0.5 * (KLDivergence2(p, m, inv_log_base) + KLDivergence2(q, m, inv_log_base));
    }
}
