
// MSVC: explicit standard includes (not pulled in transitively).
#include <vector>
#include <array>
#include <utility>
#include <cstdint>
#include "first_order_dominance_analysis.h"

#include <cmath>
#include <algorithm>


    FODAnalysis AnalyzeFOD(const Tensor<double, 2>& fod_matrix, double tolerance) {
        const std::uint64_t n = fod_matrix.Extent(0);
        
        FODAnalysis res{
            .fod_closed = MinMaxTransitiveClosure(fod_matrix),
            .approx_cells = Tensor<double, 2>(std::array<std::uint64_t, 2>{n, n}, kUninitialized),
            .approx_tot = 0.0,
            .approx_tot_corr = 0.0,
            .unique_alphas = {}
        };
        
        double sum_orig = 0.0;
        double sum_diff = 0.0;
        
        std::vector<double> alphas;
        alphas.reserve(n * n);
        
        for (std::uint64_t r = 0; r < n; ++r) {
            for (std::uint64_t c = 0; c < n; ++c) {
                const double orig = fod_matrix(r, c);
                const double closed = res.fod_closed(r, c);
                
                // Calcolo dell'errore assoluto per la singola cella.
                const double diff = std::abs(orig - closed);
                res.approx_cells(r, c) = diff;
                
                sum_orig += orig;
                sum_diff += diff;
                
                alphas.push_back(closed);
            }
        }
        
        // Calcolo dell'indice di approssimazione globale.
        res.approx_tot = (sum_orig != 0.0) ? (sum_diff / sum_orig) : 0.0;
        
        const double denom_corr = sum_orig - static_cast<double>(n);
        res.approx_tot_corr = (std::abs(denom_corr) > tolerance) ? (sum_diff / denom_corr) : 0.0;
        
        // Ordinamento dei valori per preparare la rimozione dei duplicati.
        std::sort(alphas.begin(), alphas.end());
        
        // Estrazione dei valori unici. La lambda cattura 'tolerance' per valore
        // e considera "uguali" i double la cui differenza è inferiore alla tolleranza.
        auto last = std::unique(alphas.begin(), alphas.end(), [tolerance](double a, double b) {
            return std::abs(a - b) < tolerance;
        });
        alphas.erase(last, alphas.end());
        
        // Riduciamo la capacità allocata dal vector alla sua dimensione reale
        alphas.shrink_to_fit();
        
        res.unique_alphas = std::move(alphas);
        
        return res;
    }
