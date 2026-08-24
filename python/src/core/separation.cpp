
// MSVC: explicit standard includes (not pulled in transitively).
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <cstdint>
#include <cstddef>
#include "separation.h"

#include <cmath>
#include <numeric>
#include <map>
#include <array>
#include <tuple>

// ===========================================================================
// Namespace Anonimo (Helper HPC Interni & Matematica)
// ===========================================================================
// Tutto ciò che è definito qui dentro ha "internal linkage". Non è visibile
// fuori da questo file, forzando l'inlining ed evitando collisioni di nomi (ODR).
namespace {
    
    // ---------------------------------------------------------------------------
    // Precalcolo dei Fattoriali a Compile-Time (Sostituisce std::tgamma)
    // In HPC, calcolare la funzione Gamma nei cicli annidati è un collo di bottiglia.
    // Usiamo una Look-Up Table (LUT) constexpr per estrarre i fattoriali in O(1).
    // ---------------------------------------------------------------------------
    constexpr std::array<double, 171> GenerateFactorials() {
        std::array<double, 171> facts{};
        facts[0] = 1.0;
        for (int i = 1; i <= 170; ++i) {
            facts[i] = facts[i - 1] * i;
        }
        return facts;
    }
    
    constexpr auto kFactorials = GenerateFactorials();
    
    // Wrapper inline per estrarre in O(1) il fattoriale, con fallback a std::tgamma
    // nel rarissimo caso in cui si superi 170! (che comunque darebbe overflow in un double).
    inline double FastFactorial(std::uint64_t n) noexcept {
        return n < 171 ? kFactorials[n] : std::tgamma(static_cast<double>(n + 1));
    }
    
    // ---------------------------------------------------------------------------
    // LexSeparationPair: Cuore matematico ottimizzato.
    // Elimina la vecchia logica che allocava std::unordered_set nell'heap.
    // Ora esegue tutto interamente nei registri della CPU.
    // ---------------------------------------------------------------------------
    inline std::tuple<double, double, double> LexSeparationPair(
                                                                std::uint64_t k_vars,
                                                                double m_modes,
                                                                const std::vector<std::uint64_t>& p,
                                                                const std::vector<std::uint64_t>& q) noexcept
    {
        
        double k1 = 0.0; // Numero di elementi in cui p > q
        double k2 = 0.0; // Numero di elementi in cui p == q
        double k3 = 0.0; // Numero di elementi in cui p < q
        
        double t_plus = 0.0;
        double t_minus = 0.0;
        
        // Singolo passaggio sui profili per estrarre tutte le metriche.
        // NOTA HPC: L'uso di vector garantisce un accesso sequenziale (cache L1 hit).
        for (std::size_t i = 0; i < p.size(); ++i) {
            if (p[i] < q[i]) {
                k3 += 1.0;
                // Preserviamo la semantica originale: (p < q) produrrà un accumulo negativo in t_minus
                t_minus += (static_cast<double>(p[i]) - static_cast<double>(q[i]));
            } else if (p[i] == q[i]) {
                k2 += 1.0;
            } else {
                k1 += 1.0;
                t_plus += (static_cast<double>(p[i]) - static_cast<double>(q[i]));
            }
        }
        
        double f1 = 0.0;
        double f2 = 0.0;
        double k = static_cast<double>(k_vars);
        
        // Calcolo della parte combinatoria
        if (k2 <= k - 1.0) {
            std::uint64_t u_k2 = static_cast<std::uint64_t>(k2);
            for (std::uint64_t j = 0; j <= u_k2; ++j) {
                double num1 = FastFactorial(k_vars - j - 1);
                double den = FastFactorial(u_k2 - j);
                double prod1 = std::pow(m_modes, k - static_cast<double>(j) - 1.0);
                f1 += ((num1 / den) * prod1);
                
                if (u_k2 <= k_vars - 2) {
                    double num2 = num1 / (k - static_cast<double>(j) - 1.0);
                    double inner_sum = 0.0;
                    for (std::uint64_t h = j + 2; h <= k_vars; ++h) {
                        inner_sum += std::pow(m_modes, k - static_cast<double>(h));
                    }
                    f2 += ((num2 / den) * inner_sum);
                }
            }
        }
        
        double k2_fatt = FastFactorial(static_cast<std::uint64_t>(k2));
        double k_fatt = FastFactorial(k_vars);
        
        // Formule esatte trasposte dal codice originale
        double sep_all = (k2_fatt * (((t_plus - t_minus) * f1) + ((t_plus * (k1 - k3 - 1.0) + t_minus * (k1 - k3 + 1.0)) * f2))) / k_fatt;
        double sep_upper = (k2_fatt * (t_plus * f1 + (t_plus * (k1 - 1.0) + t_minus * k1) * f2)) / k_fatt;
        double sep_lower = (-k2_fatt * (t_minus * f1 + (t_minus * (k3 - 1.0) + t_plus * k3) * f2)) / k_fatt;
        
        return {sep_all, sep_lower, sep_upper};
    }
    
} // namespace anonimo

// ===========================================================================
// Implementazione API Pubblica
// ===========================================================================

// ***********************************************
// LexSeparationEqDeg
// ***********************************************
LexSeparationResult LexSeparationEqDeg(std::uint64_t numero_variabili, std::uint64_t numero_modalita) {
    std::uint64_t numero_profili = static_cast<std::uint64_t>(std::pow(numero_modalita, numero_variabili));
    
    // Inizializzazione basata sulla Move Semantics:
    // nessun puntatore a shared_ptr, ma istanze valide gestite in stack
    LexSeparationResult result{
        .sep_all        = Tensor<double, 2>(std::array<std::uint64_t, 2>{numero_profili, numero_profili}, 0.0),
        .sep_lower      = Tensor<double, 2>(std::array<std::uint64_t, 2>{numero_profili, numero_profili}, 0.0),
        .sep_upper      = Tensor<double, 2>(std::array<std::uint64_t, 2>{numero_profili, numero_profili}, 0.0),
        .sep_vertical   = Tensor<double, 2>(std::array<std::uint64_t, 2>{numero_profili, numero_profili}, 0.0),
        .sep_horizontal = Tensor<double, 2>(std::array<std::uint64_t, 2>{numero_profili, numero_profili}, 0.0),
        .profili        = std::vector<std::vector<std::uint64_t>>(numero_profili, std::vector<std::uint64_t>(numero_variabili, 0))
    };
    
    // Generazione profili: Algoritmo 'Odometer' per mantenere il corretto
    // ordinamento lessicografico richiesto dal modello.
    std::uint64_t modalita = 1;
    for (std::uint64_t profilo_id = 1; profilo_id < numero_profili; ++profilo_id) {
        auto& profilo = result.profili[profilo_id];
        const auto& profilo_prec = result.profili[profilo_id - 1];
        
        if (modalita == numero_modalita) {
            std::int64_t p = static_cast<std::int64_t>(numero_variabili) - 2;
            for (; p >= 0; --p) {
                if (profilo_prec[p] != numero_modalita - 1) {
                    profilo[p] = profilo_prec[p] + 1;
                    for (std::int64_t k = 0; k < p; ++k) {
                        profilo[k] = profilo_prec[k];
                    }
                    break;
                }
            }
            profilo[numero_variabili - 1] = 0;
            modalita = 1;
        } else {
            for (std::uint64_t p = 0; p < numero_variabili - 1; ++p) {
                profilo[p] = profilo_prec[p];
            }
            profilo[numero_variabili - 1] = modalita;
            ++modalita;
        }
    }
    
    double m_val = static_cast<double>(numero_modalita);
    
    // Doppio ciclo per popolare le matrici.
    // L'uso di Structured Binding (C++17) qui estrae la tupla a costo zero.
    for (std::uint64_t pi = 0; pi < numero_profili; ++pi) {
        for (std::uint64_t qi = 0; qi < numero_profili; ++qi) {
            const auto& p = result.profili[pi];
            const auto& q = result.profili[qi];
            
            auto [s_all, s_low, s_up] = LexSeparationPair(numero_variabili, m_val, p, q);
            
            result.sep_all(pi, qi) = s_all;
            result.sep_lower(pi, qi) = s_low;
            result.sep_upper(pi, qi) = s_up;
            result.sep_vertical(pi, qi) = std::abs(s_low - s_up);
            result.sep_horizontal(pi, qi) = s_all - result.sep_vertical(pi, qi);
        }
    }
    
    return result;
}

// ***********************************************
// LexSeparationDeg
// ***********************************************
LexSeparationResult LexSeparationDeg(const std::vector<std::uint64_t>& numero_modalita) {
    std::uint64_t numero_profili = std::accumulate(numero_modalita.begin(), numero_modalita.end(), 1ULL, std::multiplies<std::uint64_t>());
    std::uint64_t numero_variabili = numero_modalita.size();
    
    LexSeparationResult result{
        .sep_all        = Tensor<double, 2>(std::array<std::uint64_t, 2>{numero_profili, numero_profili}, 0.0),
        .sep_lower      = Tensor<double, 2>(std::array<std::uint64_t, 2>{numero_profili, numero_profili}, 0.0),
        .sep_upper      = Tensor<double, 2>(std::array<std::uint64_t, 2>{numero_profili, numero_profili}, 0.0),
        .sep_vertical   = Tensor<double, 2>(std::array<std::uint64_t, 2>{numero_profili, numero_profili}, 0.0),
        .sep_horizontal = Tensor<double, 2>(std::array<std::uint64_t, 2>{numero_profili, numero_profili}, 0.0),
        .profili        = std::vector<std::vector<std::uint64_t>>(numero_profili, std::vector<std::uint64_t>(numero_variabili, 0))
    };
    
    // Lambda per calcolare l'LLE basandosi sulla priorità delle variabili
    auto calcola_lle = [&](const std::vector<std::uint64_t>& var_priority, std::vector<std::vector<std::uint64_t>>& lle) {
        std::uint64_t modalita = 1;
        std::uint64_t last_v = var_priority.back();
        std::fill(lle[0].begin(), lle[0].end(), 0);
        
        for (std::uint64_t profilo_id = 1; profilo_id < lle.size(); ++profilo_id) {
            auto& profilo = lle[profilo_id];
            const auto& profilo_prec = lle[profilo_id - 1];
            
            if (modalita == numero_modalita[last_v]) {
                std::int64_t p = static_cast<std::int64_t>(numero_variabili) - 2;
                for (; p >= 0; --p) {
                    std::uint64_t var = var_priority[p];
                    if (profilo_prec[var] != numero_modalita[var] - 1) {
                        profilo[var] = profilo_prec[var] + 1;
                        for (std::int64_t k = 0; k < p; ++k) {
                            std::uint64_t var_k = var_priority[k];
                            profilo[var_k] = profilo_prec[var_k];
                        }
                        for (std::uint64_t k = p + 1; k < numero_variabili; ++k) {
                            std::uint64_t var_k = var_priority[k];
                            profilo[var_k] = 0;
                        }
                        break;
                    }
                }
                profilo[last_v] = 0;
                modalita = 1;
            } else {
                for (std::uint64_t p = 0; p < numero_variabili - 1; ++p) {
                    std::uint64_t var_p = var_priority[p];
                    profilo[var_p] = profilo_prec[var_p];
                }
                profilo[last_v] = modalita;
                ++modalita;
            }
        }
    };
    
    std::vector<std::uint64_t> var_priority(numero_variabili);
    std::iota(var_priority.begin(), var_priority.end(), 0);
    
    calcola_lle(var_priority, result.profili);
    
    // Mappa Mantenuta per retrocompatibilità. O(log N) lookup.
    // In futuro potrebbe essere sostituita dal calcolo matematico dell'indice Mixed-Radix
    std::map<std::vector<std::uint64_t>, std::uint64_t> posizione_profili;
    for (std::uint64_t p_id = 0; p_id < numero_profili; ++p_id) {
        posizione_profili[result.profili[p_id]] = p_id;
    }
    
    double numero_lex = FastFactorial(numero_variabili);
    std::vector<std::vector<std::uint64_t>> lle(numero_profili, std::vector<std::uint64_t>(numero_variabili, 0));
    
    // Ciclo sulle permutazioni usando l'algoritmo STL standard
    do {
        calcola_lle(var_priority, lle);
        for (std::uint64_t p_id = 0; p_id < numero_profili; ++p_id) {
            for (std::uint64_t q_id = p_id + 1; q_id < numero_profili; ++q_id) {
                double distanza = static_cast<double>(q_id - p_id);
                std::uint64_t sep_p_id = posizione_profili.at(lle[p_id]);
                std::uint64_t sep_q_id = posizione_profili.at(lle[q_id]);
                
                double ratio = distanza / numero_lex;
                result.sep_all(sep_p_id, sep_q_id) += ratio;
                result.sep_all(sep_q_id, sep_p_id) = result.sep_all(sep_p_id, sep_q_id);
                
                result.sep_lower(sep_p_id, sep_q_id) += ratio;
                result.sep_upper(sep_q_id, sep_p_id) += ratio;
            }
        }
    } while (std::next_permutation(var_priority.begin(), var_priority.end()));
    
    for (std::uint64_t p_id = 0; p_id < numero_profili; ++p_id) {
        for (std::uint64_t q_id = 0; q_id < numero_profili; ++q_id) {
            double sep_vert = std::abs(result.sep_lower(p_id, q_id) - result.sep_upper(p_id, q_id));
            result.sep_vertical(p_id, q_id) = sep_vert;
            result.sep_horizontal(p_id, q_id) = result.sep_all(p_id, q_id) - sep_vert;
        }
    }
    
    return result;
}

// ***********************************************
// LexMRP
// ***********************************************
LexMrpResult LexMrp(const std::vector<std::uint64_t>& numero_modalita) {
    std::uint64_t numero_profili = std::accumulate(numero_modalita.begin(), numero_modalita.end(), 1ULL, std::multiplies<std::uint64_t>());
    std::uint64_t numero_variabili = numero_modalita.size();
    
    LexMrpResult result{
        .mrp     = Tensor<double, 2>(std::array<std::uint64_t, 2>{numero_profili, numero_profili}, 0.0),
        .profili = std::vector<std::vector<std::uint64_t>>(numero_profili, std::vector<std::uint64_t>(numero_variabili, 0))
    };
    
    std::uint64_t modalita = 1;
    for (std::uint64_t profilo_id = 1; profilo_id < numero_profili; ++profilo_id) {
        auto& profilo = result.profili[profilo_id];
        const auto& profilo_prec = result.profili[profilo_id - 1];
        
        if (modalita == numero_modalita.back()) {
            std::int64_t p = static_cast<std::int64_t>(numero_variabili) - 2;
            for (; p >= 0; --p) {
                if (profilo_prec[p] != numero_modalita[p] - 1) {
                    profilo[p] = profilo_prec[p] + 1;
                    for (std::int64_t k = 0; k < p; ++k) {
                        profilo[k] = profilo_prec[k];
                    }
                    break;
                }
            }
            profilo.back() = 0;
            modalita = 1;
        } else {
            for (std::uint64_t p = 0; p < numero_variabili - 1; ++p) {
                profilo[p] = profilo_prec[p];
            }
            profilo.back() = modalita;
            ++modalita;
        }
    }
    
    std::uint64_t k = numero_variabili;
    
    for (std::uint64_t p_id = 0; p_id < numero_profili; ++p_id) {
        for (std::uint64_t q_id = p_id + 1; q_id < numero_profili; ++q_id) {
            const auto& p = result.profili[p_id];
            const auto& q = result.profili[q_id];
            
            // In LexMRP i parametri k1 e k2 hanno una semantica leggermente
            // diversa rispetto a LexSeparation (p < q invece di p > q).
            // Li separiamo esplicitamente per non creare confusione.
            std::uint64_t p_minore_q = 0;
            std::uint64_t p_uguale_q = 0;
            
            for (std::size_t i = 0; i < p.size(); ++i) {
                if (p[i] < q[i]) {
                    ++p_minore_q;
                } else if (p[i] == q[i]) {
                    ++p_uguale_q;
                }
            }
            
            double sum_val = 0.0;
            for (std::uint64_t s = 0; s <= p_uguale_q; ++s) {
                double k_s_1_f = FastFactorial(k - s - 1);
                double k2_s_f = FastFactorial(p_uguale_q - s);
                sum_val += (k_s_1_f / k2_s_f);
            }
            
            double k2_f = FastFactorial(p_uguale_q);
            double k_f = FastFactorial(k);
            
            result.mrp(p_id, q_id) = (static_cast<double>(p_minore_q) * k2_f * sum_val) / k_f;
            result.mrp(q_id, p_id) = 1.0 - result.mrp(p_id, q_id);
        }
    }
    
    return result;
}
