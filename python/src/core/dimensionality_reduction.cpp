
// MSVC: explicit standard includes (not pulled in transitively).
#include <string_view>
#include <memory>
#include <cstddef>
#include "display_message.h"
#include "binary_variable_poset.h"
#include "loss_function_mrp.h"
#include "loss_function_mrp_l1.h"
#include "dimensionality_reduction.h"
#include "tensor.h"


#include <utility>
#include <algorithm>
#include <vector>
#include <span>
#include <limits>
#include <numeric>
#include <chrono>
#include <cmath>
#include <thread>
#include <vector>
#include <deque>
#include <atomic>
#include <cstdint>
#include <vector>
#include <format>
#include <ranges>
#include <string>
#include <iterator>
#include <unordered_map>



// =========================================================================
// Costruttore DimensionalityReductionResult
// =========================================================================

DimensionalityReductionResult::DimensionalityReductionResult(
                                                             std::atomic<std::uint64_t>& all_le_elaborated)
    : all_le_elaborated_(all_le_elaborated) {
    
}

// =========================================================================
// Metodi di DimensionalityReductionResult
// =========================================================================

void DimensionalityReductionResult::AddLe(std::vector<std::uint64_t> le, double loss) {
    if (loss < best_loss_) {
        best_loss_ = loss;
        best_permutation_ = le;
    }
    
    le_elaborated_.push_back(std::move(le));
    loss_values_.push_back(loss);
    
    all_le_elaborated_.fetch_add(1, std::memory_order_relaxed);
}


void DimensionalityReductionResult::AddLe(DimensionalityReductionResult& other) {
    // 1. Trasferimento memoria vettori principali (rubiamo la memoria a 'other')
    le_elaborated_.insert(
                          le_elaborated_.end(),
                          std::make_move_iterator(other.le_elaborated_.begin()),
                          std::make_move_iterator(other.le_elaborated_.end())
                          );
    
    loss_values_.insert(
                        loss_values_.end(),
                        other.loss_values_.begin(),
                        other.loss_values_.end()
                        );
    
    // 2. Se il thread concluso ha trovato un minimo globale migliore, aggiorniamo il record
    if (other.best_loss_ < best_loss_) {
        best_loss_ = other.best_loss_;
        
        // Usiamo std::move per rubare la memoria anche a stringhe/vettori pesanti
        best_permutation_ = std::move(other.best_permutation_);
        best_profile_results_ = std::move(other.best_profile_results_);
    }
}

// =========================================================================
// Ricostruzione Risultati Finali
// =========================================================================

void DimensionalityReductionResult::BuildBestProfileResults(
                                                            const std::unordered_map<std::uint64_t, double>& weights,
                                                            LossFunctionMRPV2& loss_func,
                                                            Tensor<double, 2>& matrice,
                                                            const std::vector<std::vector<std::uint64_t>>& list_a,
                                                            const std::vector<std::vector<std::uint64_t>>& list_b) {
    
    // 1. Nessun bisogno di estrarre da un iteratore:
    // best_permutation_ è già un std::vector mantenuto aggiornato.
    const std::size_t num_vars = best_permutation_.size();
    
    // 2. Inversione cache-friendly e standard (Costo O(N), nessuna allocazione manuale)
    std::vector<std::uint64_t> best_permutation_inv(
                                                    best_permutation_.rbegin(),
                                                    best_permutation_.rend()
                                                    );
    
    // 3. Allocazione locale sullo Stack dei vettori
    std::vector<std::uint64_t> le;
    std::vector<std::uint64_t> le_inv;
    
    // Assumo che questa funzione accetti reference standard (std::vector<uint64_t>&)
    DimensionalityReductionBuildLE(num_vars, best_permutation_, best_permutation_inv, le, le_inv);
    
    // 4. Invece di modificare liste con push_front e pop_front (lento e non thread-safe sulle read-only),
    // creiamo vettori temporanei lineari combinando la nuova LE con le precedenti.
    std::vector<std::vector<std::uint64_t>> current_rows;
    current_rows.reserve(1 + list_a.size());
    current_rows.push_back(le);
    current_rows.insert(current_rows.end(), list_a.begin(), list_a.end());
    
    std::vector<std::vector<std::uint64_t>> current_cols;
    current_cols.reserve(1 + list_b.size());
    current_cols.push_back(le_inv);
    current_cols.insert(current_cols.end(), list_b.begin(), list_b.end());
    
    // 5. Valutazione MRP utilizzando vettori invece di liste
    std::vector<std::pair<std::uint64_t, double>> result_error;
    loss_func(matrice, current_rows, current_cols, result_error);
    
    // 6. Lambda ottimizzata (Catturiamo solo il necessario by reference)
    // Usiamo tipi espliciti e costanti per ottimizzare le istruzioni di bit-shift della CPU
    auto calcola_posizione = [&](std::uint64_t profilo) -> std::pair<std::uint64_t, std::uint64_t> {
        std::uint64_t pos_lex = 0;
        std::uint64_t pos_lex_inv = 0;
        
        for (std::size_t i = 0; i < num_vars; ++i) {
            const std::uint64_t offset_dest = num_vars - i - 1;
            
            // Elaborazione riga
            const std::uint64_t var = best_permutation_[i];
            const std::uint64_t val = (profilo >> (num_vars - var - 1)) & 1ULL;
            pos_lex |= (val << offset_dest);
            
            // Elaborazione colonna
            const std::uint64_t var_inv = best_permutation_inv[i];
            const std::uint64_t val_inv = (profilo >> (num_vars - var_inv - 1)) & 1ULL;
            pos_lex_inv |= (val_inv << offset_dest);
        }
        return {pos_lex, pos_lex_inv};
    };
    
    // 7. Loop di costruzione finale con "Structured Bindings" (C++17/20)
    for (const auto& [profilo, errore] : result_error) {
        auto [pos_x, pos_y] = calcola_posizione(profilo);
        
        double peso = weights.at(profilo);
        
        best_profile_results_[profilo] = {pos_x + 1, pos_y + 1, peso, errore};
    }
}


// =========================================================================
// Costruzione Estensione Lineare (Bit Permutation)
// =========================================================================

void DimensionalityReductionBuildLE(
                                    std::uint64_t num_vars,
                                    std::span<const std::uint64_t> permutation,
                                    std::span<const std::uint64_t> permutation_inv,
                                    std::vector<std::uint64_t>& le,
                                    std::vector<std::uint64_t>& le_inv) {
    
    const std::uint64_t domain_size = 1ULL << num_vars;
    
    // 2. Pre-allocazione e riempimento massivo (Vettorializzato dal compilatore)
    const std::uint64_t k_invalid_val = std::numeric_limits<std::uint64_t>::max();
    le.assign(domain_size, k_invalid_val);
    le_inv.assign(domain_size, k_invalid_val);
    
    // 3. Iterazione sugli elementi utilizzati
    for (std::uint64_t n = 0; n < domain_size; ++n) {
        std::uint64_t pos_lex = 0;
        std::uint64_t pos_lex_inv = 0;
        
        // 4. Calcolo delle posizioni con bit-masking branchless
        for (std::size_t i = 0; i < num_vars; ++i) {
            // Isoliamo l'i-esimo bit partendo da sinistra
            const std::uint64_t bit = (n >> (num_vars - i - 1)) & 1ULL;
            
            // Riposizioniamo il bit in base alle permutazioni fornite
            // L'operazione bitwise OR (|) accumula i bit nelle nuove posizioni.
            // Costo CPU: bassissimo, senza salti logici (branchless).
            pos_lex     |= (bit << (num_vars - permutation[i] - 1));
            pos_lex_inv |= (bit << (num_vars - permutation_inv[i] - 1));
        }
        
        // 5. Scrittura diretta in memoria (Accesso casuale O(1))
        le[pos_lex] = n;
        le_inv[pos_lex_inv] = n;
    }
}


void DimensionalityReductionProcessPermutationGroup(
                                                    std::uint64_t lb_value,
                                                    std::uint64_t up_value,
                                                    std::uint64_t num_vars,
                                                    LossFunctionMRPV2& loss_func,
                                                    Tensor<double, 2>& mrp_rf,
                                                    const std::vector<std::vector<std::uint64_t>>& base_rows,
                                                    const std::vector<std::vector<std::uint64_t>>& base_cols,
                                                    std::uint64_t& num_le_elaborated,
                                                    DimensionalityReductionResult& result) {
    
    // 1. Allocazioni fisse spostate FUORI dal ciclo
    std::vector<std::uint64_t> permutation(num_vars);
    std::vector<std::uint64_t> permutation_inv(num_vars);
    std::vector<std::uint64_t> internal_permutation(num_vars - 2);
    
    // Popolamento delle variabili interne
    for (std::uint64_t i = 0, j = 0; i < num_vars; ++i) {
        if (i != lb_value && i != up_value) {
            internal_permutation[j++] = i;
        }
    }
    
    // 2. Vettori per le estensioni lineari: pre-allocati!
    std::vector<std::uint64_t> le;
    std::vector<std::uint64_t> le_inv;
    
    // 3. Flattening delle liste: creiamo array contigui pronti per la loss_func.
    // L'indice 0 ospiterà la LE corrente (simulando il push_front originale).
    std::vector<std::vector<std::uint64_t>> current_rows(1 + base_rows.size());
    std::vector<std::vector<std::uint64_t>> current_cols(1 + base_cols.size());
    
    // Copiamo le funzioni "base" a partire dall'indice 1
    for (std::size_t i = 0; i < base_rows.size(); ++i) current_rows[i + 1] = base_rows[i];
    for (std::size_t i = 0; i < base_cols.size(); ++i) current_cols[i + 1] = base_cols[i];
    
    const std::size_t internal_size = internal_permutation.size();
    
    do {
        // Costruzione permutazione
        permutation[0] = lb_value;
        permutation_inv[0] = up_value;
        
        for (std::size_t p = 0; p < internal_size; ++p) {
            permutation[p + 1] = internal_permutation[p];
            permutation_inv[p + 1] = internal_permutation[internal_size - p - 1];
        }
        
        permutation[num_vars - 1] = up_value;
        permutation_inv[num_vars - 1] = lb_value;
        
        // Costruzione delle LE.
        // Zero allocazioni: .assign() riuserà la capacità esistente di le e le_inv.
        DimensionalityReductionBuildLE(num_vars, permutation, permutation_inv, le, le_inv);

        current_rows[0] = std::move(le);
        current_cols[0] = std::move(le_inv);
        
        // Valutazione matematica
        double loss_value = loss_func(mrp_rf, current_rows, current_cols);
        
        // Salvataggio nel Result (AddLe copierà `permutation`, ed è giusto
        // perché ci serve intatta per il prossimo ciclo std::next_permutation)
        result.AddLe(permutation, loss_value);
        
        // RIPRISTINO: "Riprendiamoci" i vettori vuoti/allocati spostandoli indietro!
        // Così nel ciclo successivo DimensionalityReductionBuildLe avrà già la memoria pronta.
        le = std::move(current_rows[0]);
        le_inv = std::move(current_cols[0]);
        
        ++num_le_elaborated;
        
    } while (std::next_permutation(internal_permutation.begin(), internal_permutation.end()));
}

void ExactDimensionalityReductionThreads(
                                         const std::unordered_map<std::uint64_t, double>& weights,
                                         std::uint64_t num_vars,
                                         LossFunctionMRPV2& loss_func,
                                         Tensor<double, 2>& mrp_rf,
                                         const std::vector<std::vector<std::uint64_t>>& base_rows,
                                         const std::vector<std::vector<std::uint64_t>>& base_cols,
                                         DisplayMessage* display_message,
                                         double thread_percentage,
                                         DimensionalityReductionResult& result) {

    // 1. Inizializzazione e Timing
    auto start_time = std::clock();
    auto t_start = std::chrono::high_resolution_clock::now();
    
    // Lettura dei core hardware disponibili (con fallback a 1 se l'OS non risponde)
    unsigned int hw_concurrency = std::thread::hardware_concurrency();
    const std::uint64_t hw_threads = (hw_concurrency == 0) ? 1ULL : static_cast<std::uint64_t>(hw_concurrency);
    
    // Logica per i thread: se la percentuale è <= 0, forziamo a 1 thread singolo.
    std::uint64_t target_threads = 1;
    if (thread_percentage > 0.0) {
        target_threads = std::max(static_cast<std::uint64_t>(1), static_cast<std::uint64_t>(std::floor(hw_threads * thread_percentage)));
    }
    
    const std::uint64_t total_pairs = (num_vars * (num_vars - 1)) / 2;
    const std::uint64_t permutations_x_pair = static_cast<std::uint64_t>(std::tgamma(num_vars - 1));
    const std::uint64_t total_le = permutations_x_pair * total_pairs;
    
    // Limitiamo i thread al numero di task effettivi
    const std::uint64_t actual_threads = std::min(target_threads, total_pairs);
    
    
    double cpu_seconds = 0.0;
    double wall_seconds = 0.0;
    
    std::jthread output_thread([&](std::stop_token stoken) {
        display_message->Display("Variables: " + std::to_string(num_vars));
        display_message->Display("Linear Extentions: " + std::to_string(total_le));
        display_message->Display("Hardware threads detected: " + std::to_string(hw_threads));
        display_message->Display("Number of threads (Actual): " + std::to_string(actual_threads));
        display_message->Start();
        
        // Pre-calcoliamo i millisecondi
        auto sleep_duration = std::chrono::seconds(display_message->OutputSeconds());
        
        while (!stoken.stop_requested()) {
            display_message->Display();
            std::this_thread::sleep_for(sleep_duration);
        }
        display_message->Stop();
        display_message->Display("------------------------------------------");
        display_message->Display("Tempo totale (CPU):  " + std::to_string(cpu_seconds) + "s");
        display_message->Display("Tempo totale (Wall): " + std::to_string(wall_seconds) + "s");
        
        // Calcoliamo il fattore di parallelizzazione
        if (wall_seconds > 0) {
            double efficiency = (cpu_seconds / wall_seconds);
            display_message->Display("Multi-core Efficiency: " + std::to_string(efficiency) + "x");
        }
        
    });
    
    // =====================================================================
    // GENERAZIONE TASK E CODA CONDIVISA
    // =====================================================================
    
    // Creiamo una lista piatta di tutte le coppie (v1, v2) che devono essere elaborate
    std::vector<std::pair<std::uint64_t, std::uint64_t>> tasks;
    tasks.reserve(total_pairs);
    for (std::uint64_t v1 = 0; v1 < num_vars - 1; ++v1) {
        for (std::uint64_t v2 = v1 + 1; v2 < num_vars; ++v2) {
            tasks.emplace_back(v1, v2);
        }
    }
    
    // Indice atomico globale. Rappresenta il "prossimo lavoro disponibile" nella coda.
    // L'hardware garantirà che nessun thread legga mai lo stesso numero.
    std::atomic<std::size_t> next_task_idx{0};
    
    // =====================================================================
    // ESECUZIONE (WORKER POOL)
    // =====================================================================
    
    std::vector<std::jthread> threads;
    threads.reserve(actual_threads);
    std::deque<std::unique_ptr<DimensionalityReductionResult>> thread_results;
    
    for (std::size_t i = 0; i < actual_threads; ++i) {
        thread_results.push_back(result.SpawnWorkerResult());
        
        threads.emplace_back([&, worker_id = i]() {
            auto& local_result = *(thread_results[worker_id]);
            std::uint64_t local_le_count = 0;
            
            std::vector<std::vector<std::uint64_t>> local_rows = base_rows;
            std::vector<std::vector<std::uint64_t>> local_cols = base_cols;
            while (true) {
                std::size_t my_task_idx = next_task_idx.fetch_add(1, std::memory_order_relaxed);
                
                if (my_task_idx >= tasks.size()) {
                    break;
                }
                
                auto [v1, v2] = tasks[my_task_idx];
                
                // 4. Eseguo il lavoro intensivo!
                DimensionalityReductionProcessPermutationGroup(
                    v1, v2, num_vars, loss_func, mrp_rf,
                    local_rows, local_cols, local_le_count, local_result
                );
            }
        });
    }
    
    // Aspettiamo che tutti i worker finiscano la coda di task
    threads.clear();
    
    // =====================================================================
    // MERGE E CHIUSURA
    // =====================================================================
    
    for (auto& local_result_ptr : thread_results) {
        result.AddLe(*local_result_ptr);
    }
    
    result.BuildBestProfileResults(weights, loss_func, mrp_rf, base_rows, base_cols);
        
    auto end_time = std::clock();
    auto t_end = std::chrono::high_resolution_clock::now();
    
    // 3. Calcolo Tempo CPU (Somma del tempo speso da tutti i core)
    // std::clock() in C++ multithread restituisce il tempo totale di processore
    cpu_seconds = static_cast<double>(end_time - start_time) / CLOCKS_PER_SEC;
    
    // 4. Calcolo Tempo Wall (Tempo effettivo passato "sull'orologio a muro")
    auto duration_wall = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);
    wall_seconds = duration_wall.count() / 1000.0;
    
    output_thread.request_stop();
}

void ExactDimensionalityReduction(
                                  const std::unordered_map<std::uint64_t, double>& weights,
                                  std::uint64_t numero_variabili,
                                  std::string_view loss_str,
                                  int lpom_type,
                                  DisplayMessage* display_message,
                                  double thread_percentage,
                                  DimensionalityReductionResult& result) {
    
    // 1. nomi_variabili per BinaryVariablePOSet
    std::vector<std::string> nomi_variabili;
    nomi_variabili.reserve(numero_variabili);
    for (std::uint64_t i = 0; i < numero_variabili; ++i) {
        nomi_variabili.push_back(std::format("v{}", i + 1));
    }
    
    auto poset = BinaryVariablePOSet::Build(nomi_variabili);
    // 2. Costruzione Matrice di partenza (MRP Start Poset)
    Tensor<double, 2> mrp_start_poset = (lpom_type == 0 ?
                                       poset->BLSDominanceAbsolute():
                                       poset->BLSDominanceRelative());
    
  
    // 3. Inizializzazione della Loss Function
    std::vector<Tensor<double, 2>> reference_matrices;
    reference_matrices.push_back(std::move(mrp_start_poset));
    
    std::shared_ptr<LossFunctionMRPV2> lfmrp = nullptr;
    
    if (loss_str == "LB") {
        lfmrp = std::make_shared<LBMRP2>(reference_matrices, weights);
    }
    else {
        throw MyException(std::format("Loss function error!: {}", loss_str));
    }
    
    // =====================================================================
    // 4. PREPARAZIONE DATI CONDIVISI
    // =====================================================================
    const std::uint64_t num_profiles = 1ULL << numero_variabili;
    
    std::vector<std::uint64_t> le_rf;
    std::vector<std::uint64_t> le_inv_rf;
    std::vector<std::uint64_t> permutazione_rf(numero_variabili);
    std::vector<std::uint64_t> permutazione_inv_rf(numero_variabili);
    
    for (std::uint64_t i = 0; i < numero_variabili; ++i) {
        permutazione_rf[i] = i;
        permutazione_inv_rf[i] = numero_variabili - i - 1;
    }
    
    DimensionalityReductionBuildLE(numero_variabili, permutazione_rf, permutazione_inv_rf, le_rf, le_inv_rf);
    const std::vector<std::vector<std::uint64_t>> base_rows = { le_rf };
    const std::vector<std::vector<std::uint64_t>> base_cols = { le_inv_rf };
    
    // 1. ESTRAZIONE PROFILI ATTIVI
    /*std::vector<std::uint64_t> elements_used;
    elements_used.reserve(weights.size());
    for (const auto& [profilo, peso] : weights) {
        elements_used.push_back(profilo);
    }*/
    std::vector<std::uint64_t> elements_used(num_profiles);
    std::iota(elements_used.begin(), elements_used.end(), 0);
    
    Tensor<double, 2> mrp_rf({num_profiles, num_profiles}, 0);

    if (lpom_type == 0) {
        DimensionalityReductionBuildMRPIntersection<LPOMStrategy::Absolute>(permutazione_rf, le_rf, le_inv_rf, mrp_rf, elements_used);
    } else if (lpom_type == 1) {
        DimensionalityReductionBuildMRPIntersection<LPOMStrategy::Relative>(permutazione_rf, le_rf, le_inv_rf, mrp_rf, elements_used);
    } else {
        throw MyException("LPOM Invalid strategy");
    }
    
    
    // 5. Chiamata al motore di riduzione multithread
    ExactDimensionalityReductionThreads(
                                        weights,
                                        numero_variabili,
                                        *lfmrp,
                                        mrp_rf,
                                        base_rows,
                                        base_cols,
                                        display_message,
                                        thread_percentage,
                                        result
                                        );
}


void BidimentionalPosetRepresentation(const std::unordered_map<std::uint64_t, double>& weights,
                                      std::uint64_t numero_variabili,
                                      std::string_view loss_str,
                                      int lpom_type,
                                      std::vector<std::uint64_t>& variable_priority,
                                      DimensionalityReductionResult& result) {
    
    // 1. nomi_variabili per BinaryVariablePOSet
    std::vector<std::string> nomi_variabili;
    nomi_variabili.reserve(numero_variabili);
    for (std::uint64_t i = 0; i < numero_variabili; ++i) {
        nomi_variabili.push_back(std::format("v{}", i + 1));
    }
    
    auto poset = BinaryVariablePOSet::Build(nomi_variabili);
    // 2. Costruzione Matrice di partenza (MRP Start Poset)
    Tensor<double, 2> mrp_start_poset = (lpom_type == 0 ?
                                       poset->BLSDominanceAbsolute():
                                       poset->BLSDominanceRelative());
    
    // 3. Inizializzazione della Loss Function
    std::vector<Tensor<double, 2>> reference_matrices;
    reference_matrices.push_back(std::move(mrp_start_poset));
    
    std::shared_ptr<LossFunctionMRPV2> loss_func = nullptr;
    
    if (loss_str == "LB") {
        loss_func = std::make_shared<LBMRP2>(reference_matrices, weights);
    }
    else {
        throw MyException(std::format("Loss function error!: {}", loss_str));
    }
    
    // =====================================================================
    // 4. PREPARAZIONE DATI CONDIVISI
    // =====================================================================
    const std::uint64_t num_profiles = 1ULL << numero_variabili;
    
    std::vector<std::uint64_t> le_rf;
    std::vector<std::uint64_t> le_inv_rf;
    std::vector<std::uint64_t> permutazione_rf(numero_variabili);
    std::vector<std::uint64_t> permutazione_inv_rf(numero_variabili);
    
    for (std::uint64_t i = 0; i < numero_variabili; ++i) {
        permutazione_rf[i] = i;
        permutazione_inv_rf[i] = numero_variabili - i - 1;
    }
    
    DimensionalityReductionBuildLE(numero_variabili, permutazione_rf, permutazione_inv_rf, le_rf, le_inv_rf);

    // 1. ESTRAZIONE PROFILI ATTIVI
    /*std::vector<std::uint64_t> elements_used;
    elements_used.reserve(weights.size());
    for (const auto& [profilo, peso] : weights) {
        elements_used.push_back(profilo);
    }*/
    std::vector<std::uint64_t> elements_used(num_profiles);
    std::iota(elements_used.begin(), elements_used.end(), 0);
    
    
    Tensor<double, 2> mrp_rf({num_profiles, num_profiles}, 0.0);

    if (lpom_type == 0) {
        DimensionalityReductionBuildMRPIntersection<LPOMStrategy::Absolute>(permutazione_rf, le_rf, le_inv_rf, mrp_rf, elements_used);
    } else if (lpom_type == 1) {
        DimensionalityReductionBuildMRPIntersection<LPOMStrategy::Relative>(permutazione_rf, le_rf, le_inv_rf, mrp_rf, elements_used);
    } else {
        throw MyException("LPOM Invalid strategy");
    }

    std::vector<std::uint64_t> le;
    std::vector<std::uint64_t> le_inv;
    std::vector<std::uint64_t> variable_priority_inv(variable_priority.size());
    for (std::uint64_t p = 0; p < variable_priority.size(); p++) {
        variable_priority_inv[p] = variable_priority[variable_priority.size() - p - 1];
        
    }
    DimensionalityReductionBuildLE(numero_variabili, variable_priority, variable_priority_inv, le, le_inv);
    
    std::vector<std::vector<std::uint64_t>> base_rows = {le_rf };
    std::vector<std::vector<std::uint64_t>> base_cols = {le_inv_rf };
    const std::vector<std::vector<std::uint64_t>> current_rows = { le, le_rf };
    const std::vector<std::vector<std::uint64_t>> current_cols = { le, le_inv_rf };

    double loss_value = (*loss_func)(mrp_rf, current_rows, current_cols);
    result.AddLe(variable_priority, loss_value);

    result.BuildBestProfileResults(weights, *loss_func, mrp_rf, base_rows, base_cols);
}

