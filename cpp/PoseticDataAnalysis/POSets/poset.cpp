
// MSVC: explicit standard includes (not pulled in transitively).
#include <string>
#include <string_view>
#include <array>
#include <optional>
#include <utility>
#include <cstdint>
#include <cstddef>
#include <functional>
#include <tuple>
#include <format>
#include <vector>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <stack>
#include <queue>
#include <memory>
#include <condition_variable>
#include <chrono>
#include <stdexcept>

#include "tensor.h"
#include "linear_extension.h"
#include "function_linear_extension.h"
#include "poset.h"
#include "linear_extension_generator.h"
#include "linear_extension_generator_tree_of_ideals.h"
#include "linear_extension_generator_bubley_dyer.h"
#include "evaluation_chain_source_bubley_dyer_from_le_matrix.h"
#include "function_linear_extension_mutual_ranking_probability.h"
#include "score.h"

// ***********************************************
// ***********************************************
// ***********************************************

std::unique_ptr<POSet> POSet::Build(const POSet::BEC& p, bool transitiveClosure) {
    std::unique_ptr<POSet> r(new POSet());
    r->FillBaseAttributes(std::get<1>(p), std::get<2>(p), nullptr, transitiveClosure);
    return r;
}

// ***********************************************
// ***********************************************
// ***********************************************

std::unique_ptr<POSet> POSet::Build(const std::vector<std::string>& elements,
                                    const std::vector<std::pair<std::string, std::string>>& comparabilities,
                                    bool transitiveClosure) {
    std::unique_ptr<POSet> r(new POSet());
    // Avendo unificato il metodo, passiamo semplicemente 'nullptr' per il parametro score
    r->FillBaseAttributes(elements, comparabilities, nullptr, transitiveClosure);
    return r;
}


/**
 * @brief Costruttore di Copia di POSet.
 * @details Esegue una copia profonda (deep copy) delle strutture dati contigue (CSR e BitSet)
 * e dei dizionari. I puntatori alle strutture derivate (Tree e Lattice of Ideals)
 * vengono deliberatamente non copiati e reinizializzati a nullptr (Lazy Evaluation).
 */
POSet::POSet(const POSet& other)
    : up_row_ptr(other.up_row_ptr),
    up_col_idx(other.up_col_idx),
    down_row_ptr(other.down_row_ptr),
    down_col_idx(other.down_col_idx),
    upset_bits_(other.upset_bits_),
    downset_bits_(other.downset_bits_),
    eid_vs_ename(other.eid_vs_ename),
    ename_vs_eid(other.ename_vs_eid),
    tree_of_ideals_(nullptr),
    lattice_of_ideals_(nullptr)
{ }

// ***********************************************
// ***********************************************
// ***********************************************

void POSet::ComputeTransitiveClosure(std::vector<std::set<std::uint64_t>>& adj) {
    const std::uint64_t n = adj.size();
    // Algoritmo di Warshall ottimizzato per liste di adiacenza (sets)
    for (std::uint64_t k = 0; k < n; ++k) {
        for (std::uint64_t i = 0; i < n; ++i) {
            // Se i raggiunge k, allora i deve raggiungere tutto ciò che k raggiunge
            if (adj[i].contains(k)) {
                // Inserimento massivo: molto più efficiente del loop manuale
                adj[i].insert(adj[k].begin(), adj[k].end());
            }
        }
    }
}

// ***********************************************
// ***********************************************
// ***********************************************

void POSet::FillBaseAttributes(const std::vector<std::string>& elements,
                               const std::vector<std::pair<std::string, std::string>>& comparabilities,
                               Score* score,
                               bool transitiveClosure) {
    
    // --- 1. RESET E PREPARAZIONE ---
    tree_of_ideals_.reset();
    lattice_of_ideals_.reset();
    
    ename_vs_eid.clear();
    ename_vs_eid.reserve(elements.size());
    eid_vs_ename.clear();
    
    // --- 2. MAPPATURA DEGLI ID ---
    if (score == nullptr) {
        // Assegnazione ID in base all'ordine posizionale nel vettore
        eid_vs_ename = elements;
        for (std::uint64_t i = 0; i < elements.size(); ++i) {
            auto [it, inserted] = ename_vs_eid.try_emplace(elements[i], i);
            if (!inserted) {
                throw std::runtime_error("POSet error: " + elements[i] + " duplicated!");
            }
        }
    } else {
        // Assegnazione ID delegata alla classe Score
        for (const auto& name : elements) {
            std::uint64_t eid = static_cast<std::uint64_t>(score->getEID(name));
            auto [it, inserted] = ename_vs_eid.try_emplace(name, eid);
            if (!inserted) {
                throw std::runtime_error("POSet error: " + name + " duplicated!");
            }
            
            // Garantiamo che il vettore degli ID inversi sia sufficientemente grande,
            // qualora lo Score fornisca ID non strettamente consecutivi o sfasati.
            if (eid >= eid_vs_ename.size()) {
                eid_vs_ename.resize(eid + 1, "");
            }
            eid_vs_ename[eid] = name;
        }
    }
    
    // Il numero massimo di nodi (la dimensione delle nostre matrici/vettori)
    // è determinato dalla grandezza finale del vettore degli ID inversi.
    const std::uint64_t num_nodes = eid_vs_ename.size();
    
    // --- 3. CREAZIONE LISTA ADIACENZA (UP) TEMPORANEA ---
    std::vector<std::set<std::uint64_t>> temp_up(num_nodes);
    
    for (const auto& comp : comparabilities) {
        std::uint64_t lower = GetElementId(comp.first);
        std::uint64_t upper = GetElementId(comp.second);
        
        if (lower == upper) continue; // Ignora auto-anelli espliciti in fase di parsing
        
        temp_up[lower].insert(upper);
        
        // Controllo Antisimmetria rapido sugli archi diretti inseriti
        if (temp_up[upper].contains(lower)) {
            throw std::runtime_error("The binary relation is not antisymmetric due to the comparability " + comp.first + " <= " + comp.second);
        }
    }
    
    // --- 4. CHIUSURA TRANSITIVA (Opzionale ma default) ---
    if (transitiveClosure) {
        ComputeTransitiveClosure(temp_up);
        
        // Controllo robusto post-chiusura: se c'è un ciclo nascosto, un nodo raggiungerà se stesso.
        for (std::uint64_t i = 0; i < num_nodes; ++i) {
            if (temp_up[i].contains(i)) {
                throw std::runtime_error("The binary relation is not a valid POSet (cycle detected after transitive closure).");
            }
        }
    }
    
    // --- 5. COSTRUZIONE COMPRESSED SPARSE ROW (CSR) ---
    
    // A) Costruzione UP CSR
    up_row_ptr.assign(num_nodes + 1, 0);
    up_col_idx.clear();
    // Ottimizzazione: pre-allochiamo la dimensione esatta stimando il numero di archi
    for (std::uint64_t i = 0; i < num_nodes; ++i) {
        up_row_ptr[i + 1] = up_row_ptr[i] + temp_up[i].size();
        for (auto val : temp_up[i]) {
            up_col_idx.push_back(val);
        }
    }
    
    // B) Inversione dei cammini per creare la struttura DOWN temporanea
    std::vector<std::set<std::uint64_t>> temp_down(num_nodes);
    for (std::uint64_t i = 0; i < num_nodes; ++i) {
        for (auto val : temp_up[i]) {
            temp_down[val].insert(i);
        }
    }
    
    // C) Costruzione DOWN CSR
    down_row_ptr.assign(num_nodes + 1, 0);
    down_col_idx.clear();
    for (std::uint64_t i = 0; i < num_nodes; ++i) {
        down_row_ptr[i + 1] = down_row_ptr[i] + temp_down[i].size();
        for (auto val : temp_down[i]) {
            down_col_idx.push_back(val);
        }
    }
    
    // --- 6. COSTRUZIONE BITSETS (UPSET E DOWNSET) ---
    // Il BitSet richiede la capacità al momento della costruzione.
    // Usiamo reserve ed emplace_back per evitare copie inutili (HPC friendly).
    upset_bits_.clear();
    downset_bits_.clear();
    upset_bits_.reserve(num_nodes);
    downset_bits_.reserve(num_nodes);
    
    for (std::uint64_t i = 0; i < num_nodes; ++i) {
        upset_bits_.emplace_back(num_nodes);
        downset_bits_.emplace_back(num_nodes);
    }
    
    for (std::uint64_t i = 0; i < num_nodes; ++i) {
        // Proprietà riflessiva: un nodo è sempre <= e >= a se stesso
        upset_bits_[i].set(i);
        downset_bits_[i].set(i);
        
        // Popola l'upset usando temp_up (elementi j tali che i <= j)
        for (auto j : temp_up[i]) {
            upset_bits_[i].set(j);
        }
        
        // Popola il downset usando temp_down (elementi j tali che j <= i)
        for (auto j : temp_down[i]) {
            downset_bits_[i].set(j);
        }
    }
}

// ***********************************************
// ***********************************************
// ***********************************************

std::string_view POSet::GetElementName(std::uint64_t idx) const {
    return eid_vs_ename[idx];
}

// ***********************************************
// ***********************************************
// ***********************************************

std::string POSet::to_string(char delimiter) const {
    const std::uint64_t n = size();
    std::string result;
    
    // HPC Ottimizzazione: Pre-allochiamo una stima ragionevole di memoria per
    // evitare continue (e lentissime) riallocazioni interne durante il +=.
    // Una matrice densa NxN con un po' di spazio per i nomi.
    result.reserve(n * n * 2 + n * 15);
    
    // Prima riga: intestazioni delle colonne
    for (std::uint64_t eid = 0; eid < n; ++eid) {
        result += delimiter;
        result += GetElementName(eid); // Usiamo il nuovo nome standardizzato
    }
    result += '\n';
    
    // Righe successive: matrice di adiacenza
    for (std::uint64_t out_eid = 0; out_eid < n; ++out_eid) {
        result += GetElementName(out_eid);
        
        // HPC Ottimizzazione: Traversata Lineare del CSR
        // Siccome le colonne del CSR sono state generate da std::set, sono ordinate!
        // Possiamo scorrerle parallelamente al loop in_eid, abbattendo la complessità.
        std::uint64_t row_start = up_row_ptr[out_eid];
        std::uint64_t row_end = up_row_ptr[out_eid + 1];
        std::uint64_t current_up_idx = row_start;
        
        for (std::uint64_t in_eid = 0; in_eid < n; ++in_eid) {
            // Diagonale principale
            if (in_eid == out_eid) {
                result += ";T";
                continue;
            }
            
            // Se l'indice CSR corrente corrisponde alla colonna in_eid, c'è un arco
            if (current_up_idx < row_end && up_col_idx[current_up_idx] == in_eid) {
                result += ";T";
                ++current_up_idx; // Avanziamo al prossimo arco di questa riga
            } else {
                // Non c'è l'arco
                result += delimiter;
            }
        }
        result += '\n';
    }
    
    return result;
}

// ***********************************************
// ***********************************************
// ***********************************************

BitSet POSet::GetImmediatePredecessors(std::uint64_t el) const {
    // Controllo di sicurezza sui limiti
    if (el >= size()) {
        throw std::runtime_error("POSet error: element out of bounds!");
    }
    
    // Inizializziamo un BitSet con capacità pari al numero totale di nodi
    BitSet immediate_preds(size());
    
    // 1. Estraiamo i limiti dal Compressed Sparse Row (CSR) per i predecessori di 'el'
    // Questo ci dà tutti gli elementi 'v' tali che v <= el
    std::uint64_t start = down_row_ptr[el];
    std::uint64_t end = down_row_ptr[el + 1];
    
    // 2. Iteriamo su tutti i candidati predecessori 'v'
    for (std::uint64_t i = start; i < end; ++i) {
        std::uint64_t v = down_col_idx[i];
        
        // Un elemento non è predecessore stretto di se stesso
        if (v == el) continue;
        
        bool is_immediate = true;
        
        // 3. Verifichiamo se 'v' è un predecessore remoto.
        // Lo è se esiste un altro predecessore 'u' di 'el' tale che v < u < el.
        for (std::uint64_t j = start; j < end; ++j) {
            std::uint64_t u = down_col_idx[j];
            
            if (u == el || u == v) continue;
            
            // Se v < u (ovvero 'u' si trova nell'UpSet di 'v'), 'v' è un predecessore
            // remoto rispetto a 'el', quindi NON è immediato.
            // HPC: Usiamo il BitSet per un check istantaneo in O(1) anziché binary_search
            if (upset_bits_[v].test(u)) {
                is_immediate = false;
                break; // Inutile controllare altri 'u', 'v' è scartato
            }
        }
        
        // Se nessun 'u' si è frapposto tra 'v' ed 'el', allora 'v' è un predecessore immediato (copertura)
        if (is_immediate) {
            immediate_preds.set(v);
        }
    }
    
    return immediate_preds;
}

// ***********************************************
// ***********************************************
// ***********************************************

POSet::DATASTORE POSet::GetImmediatePredecessors(const std::vector<std::uint64_t>& els) const {
    
    // Allochiamo il vettore esterno alla dimensione totale dei nodi del POSet.
    // Specifichiamo al std::vector di usare il costruttore BitSet(size()) per
    // istanziare ogni elemento interno.
    // In questo modo result[el] sarà sempre un accesso sicuro e O(1).
    // Gli elementi ignorati rimarranno dei BitSet con tutti i bit a 0.
    POSet::DATASTORE result(size(), BitSet(size()));
    
    // Iteriamo sugli elementi richiesti e deleghiamo il calcolo
    // all'overload altamente ottimizzato (basato su CSR e BitSet in O(1)).
    for (std::uint64_t el : els) {
        result[el] = GetImmediatePredecessors(el);
    }
    
    // Grazie alla Return Value Optimization (RVO), l'intero datastore
    // viene costruito e passato senza alcuna copia in memoria.
    return result;
}

// ***********************************************
// ***********************************************
// ***********************************************

POSet::DATASTORE POSet::GetImmediatePredecessors() const {
    const std::uint64_t n = size();
    
    // Allochiamo l'intero DATASTORE in un colpo solo.
    // Il vettore conterrà esattamente 'n' elementi, e istruiamo il compilatore
    // a costruire ogni elemento come un BitSet di capacità 'n'.
    // In questo modo pre-allochiamo l'intera matrice senza alcuna frammentazione.
    POSet::DATASTORE result(n, BitSet(n));
    
    // Iteriamo sequenzialmente su tutti gli ID
    for (std::uint64_t i = 0; i < n; ++i) {
        // Deleghiamo al metodo ottimizzato per il singolo elemento.
        // L'operatore di assegnamento (operator=) di BitSet copierà
        // istantaneamente i blocchi di memoria a 64 bit.
        result[i] = GetImmediatePredecessors(i);
    }
    
    // Grazie alla Return Value Optimization (RVO), non c'è alcuna copia aggiuntiva
    return result;
}

// ***********************************************
// ***********************************************
// ***********************************************

TreeOfIdeals* POSet::GetTreeOfIdeals() {
    if (!tree_of_ideals_) {
        
        // 1. Recupero dei predecessori immediati
        auto im_pred = GetImmediatePredecessors();
        
        // 2. Generazione della Linear Extension (niente allocazioni dinamiche, va dritta sullo stack)
        LinearExtension linear_ext(size());
        FirstLE(linear_ext);
        
        // 3. Conversione
        // 3.1 Otteniamo il numero di nodi
        std::uint64_t num_nodes = im_pred.size();
        // 3.2 Creazione e allocazione del vettore di BitSet
        std::vector<BitSet> impred_converted;
        impred_converted.reserve(num_nodes);
        for (std::uint64_t i = 0; i < num_nodes; ++i) {
            impred_converted.emplace_back(num_nodes); // Capacità del BitSet = num_nodes
        }
        
        for (std::uint64_t val = 0; val < num_nodes; ++val) {
            std::uint64_t pos_val = linear_ext.GetPos(val);
            
            for (std::uint64_t se : im_pred[val]) {
                std::uint64_t pos_se = linear_ext.GetPos(se);
                // Invece del push_back, impostiamo semplicemente il bit corrispondente in O(1)
                impred_converted[pos_val].set(pos_se);
            }
        }
        
        // 5. Creazione dell'albero
        // Spostiamo i dati localmente nel costruttore usando std::move.
        tree_of_ideals_ = std::make_unique<TreeOfIdeals>(std::move(impred_converted), std::move(linear_ext));
    }
    
    return tree_of_ideals_.get();
}

// ***********************************************
// ***********************************************
// ***********************************************


void POSet::FirstLE(LinearExtension& le) const {
    const std::uint64_t n = size();
    
    // Contatore dei "predecessori rimanenti" per ogni nodo.
    std::vector<std::uint64_t> in_degree(n);
    
    // Min-Heap: coda di priorità che tiene in cima sempre il nodo con l'ID più piccolo.
    std::priority_queue<std::uint64_t, std::vector<std::uint64_t>, std::greater<>> min_heap;
    
    // Inizializzazione O(N):
    for (std::uint64_t i = 0; i < n; ++i) {
        // Il numero di predecessori è esattamente la grandezza del DownSet nel CSR
        in_degree[i] = down_row_ptr[i + 1] - down_row_ptr[i];
        
        // Se non ha predecessori, è un candidato per l'estensione lineare
        if (in_degree[i] == 0) {
            min_heap.push(i);
        }
    }
    
    // Costruzione dell'Estensione Lineare
    std::uint64_t k = 0;
    while (!min_heap.empty() && k < le.size()) {
        // 1. Estraiamo l'elemento minimo disponibile in O(log N)
        std::uint64_t min_el = min_heap.top();
        min_heap.pop();
        
        // 2. Lo assegniamo all'estensione lineare
        // (usiamo il metodo Set in PascalCase che abbiamo aggiornato in LinearExtension)
        le.Set(k++, min_el);
        
        // 3. "Rimuoviamo" min_el dal POSet in O(E).
        // Chi era influenzato da min_el? Esattamente il suo UpSet!
        // Leggiamo i successori in memoria contigua tramite il CSR.
        std::uint64_t up_start = up_row_ptr[min_el];
        std::uint64_t up_end = up_row_ptr[min_el + 1];
        
        for (std::uint64_t idx = up_start; idx < up_end; ++idx) {
            std::uint64_t succ_el = up_col_idx[idx];
            
            // Riduciamo il grado entrante. Se scende a 0, tutti i suoi
            // predecessori sono stati "processati", quindi entra nell'heap.
            if (--in_degree[succ_el] == 0) {
                min_heap.push(succ_el);
            }
        }
    }
}


// ***********************************************
// ***********************************************
// ***********************************************

LatticeOfIdeals* POSet::GetLatticeOfIdeals() {
    // Lazy initialization: costruiamo il reticolo solo alla prima richiesta
    if (!lattice_of_ideals_) {
        
        // 1. Recuperiamo l'albero (usando il getter lazy del TreeOfIdeals)
        TreeOfIdeals* toi = GetTreeOfIdeals();
        
        // 2. Creiamo il reticolo.
        // Passiamo *toi (dereferenziato) perché il costruttore vuole una reference.
        // Usiamo GetRoot() per allinearci alla nuova nomenclatura PascalCase.
        lattice_of_ideals_ = std::make_unique<LatticeOfIdeals>(*toi, toi->GetRoot());
    }
    
    // Restituiamo il raw pointer.
    // La proprietà (ownership) resta del POSet tramite il unique_ptr.
    return lattice_of_ideals_.get();
}

// ***********************************************
// Genera la Matrice di Incidenza del POSet.
// ***********************************************

Tensor<std::uint8_t, 2> POSet::IncidenceMatrix() const {
    const std::size_t n = size(); // cast a size_t per allinearsi al costruttore di Matrice
    
    Tensor<std::uint8_t, 2> result({n, n}, 0);
    
    for (std::size_t i = 0; i < n; ++i) {
        // Relazione riflessiva (a <= a).
        // Usiamo operator() invece di .at() per bypassare il controllo dei limiti (bounds check)
        result(i, i) = 1;
        
        // Recuperiamo i limiti per la riga 'i' dal layout CSR interno del POSet.
        const std::uint64_t start_idx = up_row_ptr[i];
        const std::uint64_t end_idx   = up_row_ptr[i + 1];
        
        // Iteriamo in modo contiguo sulla memoria per riempire la riga della matrice
        for (std::uint64_t k = start_idx; k < end_idx; ++k) {
            const std::size_t j = static_cast<std::size_t>(up_col_idx[k]);
            result(i, j) = 1;
        }
    }
    
    // Il compilatore applica la NRVO (Named Return Value Optimization).
    // Non avverrà nessuna copia e l'unico_ptr interno verrà spostato a costo zero.
    return result;
}

// ***********************************************
// ***********************************************
// ***********************************************

std::vector<std::pair<std::string, std::string>> POSet::OrderRelation() const {
    std::vector<std::pair<std::string, std::string>> result;
    const std::uint64_t n = size();
    
    // Aggiungiamo 'n' per includere la proprietà riflessiva (x <= x per ogni nodo).
    result.reserve(up_col_idx.size() + n);
    
    for (std::uint64_t row = 0; row < n; ++row) {
        std::string_view row_e = GetElementName(row);
        
        // Proprietà riflessiva
        result.emplace_back(row_e, row_e);
        
        const std::uint64_t start = up_row_ptr[row];
        const std::uint64_t end = up_row_ptr[row + 1];
        
        for (std::uint64_t idx = start; idx < end; ++idx) {
            std::uint64_t col = up_col_idx[idx];
            std::string_view col_e = GetElementName(col);
            
            result.emplace_back(row_e, col_e);
        }
    }
    return result;
}

// ***********************************************
// ***********************************************
// ***********************************************

Tensor<std::uint8_t, 2> POSet::CoverMatrix() const {
    const std::uint64_t n = size();
    
    // Matrice flat pre-allocata a 0 (1D array contiguo in memoria)
    Tensor<std::uint8_t, 2> result({n, n}, 0);

    for (std::uint64_t v1 = 0; v1 < n; ++v1) {
        // Estraiamo i successori di v1 usando il CSR
        std::uint64_t start = up_row_ptr[v1];
        std::uint64_t end = up_row_ptr[v1 + 1];
        
        // Iteriamo su tutti i candidati successori 'v2'
        for (std::uint64_t i = start; i < end; ++i) {
            std::uint64_t v2 = up_col_idx[i];
            
            if (v1 == v2) continue; // Un nodo non copre se stesso
            
            bool is_immediate = true;
            
            // Cerchiamo se esiste un nodo v3 frapposto tra v1 e v2
            for (std::uint64_t j = start; j < end; ++j) {
                std::uint64_t v3 = up_col_idx[j];
                
                if (v3 == v1 || v3 == v2) continue;
                
                // Se v1 < v3 < v2, allora v2 è un successore remoto, non una copertura.
                // Controllo O(1) super rapido tramite BitSet
                if (upset_bits_[v3].test(v2)) {
                    is_immediate = false;
                    break;
                }
            }
            
            // Se non abbiamo trovato nodi intermedi, v2 copre v1.
            // Impostiamo il flag nella matrice appiattita.
            if (is_immediate) {
                result(v1, v2) = 1;
            }
        }
    }
    
    return result; // Nessuna copia grazie alla Return Value Optimization (RVO)
}


/**
 * @brief Calcola il Meet (Greatest Lower Bound - Estremo Inferiore) di un sottoinsieme di elementi.
 * * @details
 * L'algoritmo sfrutta un'equivalenza fondamentale della Lattice Theory per mappare
 * la ricerca del Meet in pure operazioni bit a bit (Intersezione e Popcount),
 * abbattendo la complessità da tempo logaritmico/lineare per nodo a O(N/64) cicli di clock.
 * * TEOREMA (Esistenza del Meet tramite ideali principali):
 * Sia P un POSet e S un suo sottoinsieme. L'insieme dei Common Lower Bounds (CLB)
 * di S è dato dall'intersezione dei downset principali di ogni x in S:
 * CLB(S) = ⋂_{x ∈ S} ↓x
 * Se il Meet m = ⋀S esiste, allora l'intero CLB coincide esattamente con il downset di m:
 * CLB(S) = ↓m
 * Poiché per ogni candidato v ∈ CLB(S) vale sempre (↓v ⊆ CLB(S)), ne consegue
 * matematicamente che l'unico elemento 'v' il cui downset ha la STESSA CARDINALITÀ
 * del CLB(S) è esattamente il Meet 'm'.
 * * @cite B. A. Davey, H. A. Priestley, "Introduction to Lattices and Order" (2nd Ed.),
 * Cambridge University Press, 2002. (Capitolo 2: Lattices and lattice homomorphisms).
 */
std::optional<std::uint64_t> POSet::Meet(const std::vector<std::uint64_t>& insieme) const {
    if (insieme.empty()) {
        return std::nullopt; // Il Meet di un insieme vuoto non è definito nel nostro contesto
    }
    
    // 1. Calcoliamo il Common Lower Bound (CLB) = ⋂_{x ∈ S} ↓x
    // Inizializziamo il CLB con il downset del primo elemento dell'insieme.
    BitSet clb = downset_bits_[insieme[0]];
    
    // Intersezione bit a bit (AND logico vettorializzato) con i downset degli altri elementi
    for (size_t i = 1; i < insieme.size(); ++i) {
        clb &= downset_bits_[insieme[i]];
    }
    
    // Calcoliamo la cardinalità |CLB(S)| tramite popcount (istruzione hardware)
    std::uint64_t clb_size = clb.count();
    
    // Se l'intersezione è vuota, non c'è alcun lower bound comune
    if (clb_size == 0) {
        return std::nullopt;
    }
    
    // 2. Cerchiamo l'elemento v ∈ CLB(S) tale che |↓v| == |CLB(S)|
    const std::uint64_t n = size();
    for (std::uint64_t v = 0; v < n; ++v) {
        // Se 'v' fa parte dei lower bound comuni...
        if (clb.test(v)) {
            // ...e il suo downset copre esattamente l'intero CLB, abbiamo trovato il Meet!
            if (downset_bits_[v].count() == clb_size) {
                return v;
            }
        }
    }
    
    // Se nessun elemento soddisfa la proprietà, S ha lower bounds massimali multipli,
    // di conseguenza il Meet unico non esiste.
    return std::nullopt;
}

// ***********************************************
// ***********************************************
// ***********************************************

/**
 * @brief Calcola il Join (Least Upper Bound - Estremo Superiore) di un sottoinsieme di elementi.
 * * @details
 * Per il principio di dualità dei reticoli (Lattice Theory), l'algoritmo è speculare a quello del Meet.
 * L'insieme dei Common Upper Bounds (CUB) di S è dato dall'intersezione degli upset principali:
 * CUB(S) = ⋂_{x ∈ S} ↑x
 * Se il Join j = ⋁S esiste, allora l'intero CUB coincide con l'upset di j:
 * CUB(S) = ↑j
 * Poiché per ogni v ∈ CUB(S) vale sempre (↑v ⊆ CUB(S)), l'unico elemento 'v' il cui upset
 * ha la STESSA CARDINALITÀ del CUB(S) è matematicamente il Join 'j'.
 * * @cite B. A. Davey, H. A. Priestley, "Introduction to Lattices and Order" (2nd Ed.),
 * Cambridge University Press, 2002.
 */
std::optional<std::uint64_t> POSet::Join(const std::vector<std::uint64_t>& insieme) const {
    if (insieme.empty()) {
        return std::nullopt; // Il Join di un insieme vuoto non è definito o è il Bottom
    }
    
    // 1. Calcoliamo il Common Upper Bound (CUB) = ⋂_{x ∈ S} ↑x
    // Inizializziamo il CUB con l'upset del primo elemento dell'insieme.
    BitSet cub = upset_bits_[insieme[0]];
    
    // Intersezione bit a bit (AND vettorializzato SIMD) con gli upset degli altri elementi
    for (size_t i = 1; i < insieme.size(); ++i) {
        cub &= upset_bits_[insieme[i]];
    }
    
    // Calcoliamo la cardinalità |CUB(S)| tramite popcount hardware
    std::uint64_t cub_size = cub.count();
    
    // Se l'intersezione è vuota, non c'è alcun upper bound comune
    if (cub_size == 0) {
        return std::nullopt;
    }
    
    // 2. Cerchiamo l'elemento v ∈ CUB(S) tale che |↑v| == |CUB(S)|
    const std::uint64_t n = size();
    for (std::uint64_t v = 0; v < n; ++v) {
        // Se 'v' fa parte degli upper bound comuni...
        if (cub.test(v)) {
            // ...e il suo upset copre esattamente l'intero CUB, abbiamo trovato il Join!
            if (upset_bits_[v].count() == cub_size) {
                return v;
            }
        }
    }
    
    // Se nessun elemento soddisfa la proprietà, S ha molteplici upper bound minimali,
    // quindi un unico Join non esiste.
    return std::nullopt;
}

// ***********************************************
// ***********************************************
// ***********************************************

const POSet::DATASTORE& POSet::UpSets() const {
    return upset_bits_;
}

/**
 * @brief Calcola il DownSet generato da un sottoinsieme di elementi.
 * @details Restituisce l'unione dei downset principali degli elementi forniti.
 */
BitSet POSet::DownSet(const std::vector<std::uint64_t>& els) const {
    const std::uint64_t n = size();
    BitSet result(n); // Pre-allocato a zero
    
    // Unione logica (OR bit a bit) dei downset di tutti gli elementi
    for (std::uint64_t el : els) {
        // Ignoriamo in modo sicuro eventuali ID fuori limite
        if (el < n) {
            result |= downset_bits_[el];
        }
    }
    
    return result;
}

/**
 * @brief Verifica se un dato sottoinsieme di elementi è un Down-Set (Ideale d'Ordine).
 * * @details
 * Un sottoinsieme S di un POSet P è un down-set se, per ogni elemento x ∈ S,
 * l'intero down-set principale di x è contenuto in S.
 * Matematicamente equivale a dire che il down-set generato da S coincide con S stesso:
 * S è un down-set ⇔ ↓S = S.
 * * L'algoritmo converte l'input S in una maschera di bit in O(|S|) e calcola ↓S in
 * blocchi SIMD sfruttando la funzione DownSet() pre-ottimizzata. Infine, valuta
 * se la differenza insiemistica (↓S \ S) è vuota.
 * * @cite B. A. Davey, H. A. Priestley, "Introduction to Lattices and Order" (2nd Ed.),
 * Cambridge University Press, 2002. (Capitolo 1: Ordered sets).
 */
bool POSet::IsDownSet(const std::vector<std::uint64_t>& els) const {
    const std::uint64_t n = size();
    
    // 1. Rappresentiamo l'insieme S come maschera di bit (O(1) per le ricerche future)
    BitSet input_bits(n);
    for (std::uint64_t el : els) {
        if (el < n) {
            input_bits.set(el);
        }
    }
    
    // 2. Calcoliamo ↓S (l'unione dei downset degli elementi di S)
    // Usiamo il metodo HPC vettorializzato appena scritto
    BitSet generated_downset = DownSet(els);
    
    // 3. Verifichiamo se ↓S \ S == ∅
    // Il metodo set_difference di BitSet esegue (generated_downset AND NOT input_bits)
    generated_downset.SetDifference(generated_downset, input_bits);
    
    // Se non ci sono bit rimanenti accesi, S non ha generato "intrusi" ed è un Down-Set perfetto!
    return generated_downset.empty();
}

/**
 * @brief Calcola l'UpSet (Filtro generato / Upper Set) da un sottoinsieme di elementi.
 * @details Restituisce l'unione degli upset principali degli elementi forniti: ↑S = ⋃_{x ∈ S} ↑x.
 */
BitSet POSet::UpSet(const std::vector<std::uint64_t>& els) const {
    const std::uint64_t n = size();
    
    // Inizializza un BitSet vuoto (tutti i bit a 0) di dimensione n
    BitSet result(n);
    
    // Unione logica (OR bit a bit) degli upset di tutti gli elementi dell'input
    for (std::uint64_t el : els) {
        // Controllo di sicurezza: processiamo solo ID validi
        if (el < n) {
            result |= upset_bits_[el];
        }
    }
    
    // Il compilatore applica la Return Value Optimization (RVO)
    // Trasferisce la proprietà della memoria senza copiare i bit.
    return result;
}

/**
 * @brief Verifica se un dato sottoinsieme di elementi è un Up-Set (Filtro/Upper Set).
 * * @details
 * Un sottoinsieme S di un POSet P è un up-set se, per ogni elemento x ∈ S,
 * l'intero up-set principale di x è contenuto in S.
 * Matematicamente equivale a dire che l'up-set generato da S coincide con S stesso:
 * S è un up-set ⇔ ↑S = S.
 * * L'algoritmo crea una maschera di bit per l'input e calcola ↑S sfruttando
 * l'unione SIMD in O(|S| * N/64) della funzione UpSet(). Successivamente calcola
 * la differenza insiemistica (↑S \ S) tramite bitwise operations. Se il risultato
 * è l'insieme vuoto, S è un up-set valido.
 * * @cite B. A. Davey, H. A. Priestley, "Introduction to Lattices and Order" (2nd Ed.),
 * Cambridge University Press, 2002. (Capitolo 1: Ordered sets).
 */
bool POSet::IsUpSet(const std::vector<std::uint64_t>& els) const {
    const std::uint64_t n = size();
    
    // 1. Rappresentiamo l'insieme in input S come maschera di bit (costo spaziale/temporale bassissimo)
    BitSet input_bits(n);
    for (std::uint64_t el : els) {
        if (el < n) {
            input_bits.set(el);
        }
    }
    
    // 2. Calcoliamo ↑S (l'unione degli upset degli elementi di S)
    // Sfrutta il metodo HPC vettorializzato che abbiamo appena implementato
    BitSet generated_upset = UpSet(els);
    
    // 3. Verifichiamo se ↑S \ S == ∅
    // Il metodo set_difference di BitSet esegue (generated_upset AND NOT input_bits)
    generated_upset.SetDifference(generated_upset, input_bits);
    
    // Se la differenza è vuota, S non ha "intrusi" verso l'alto ed è un perfetto Up-Set!
    return generated_upset.empty();
}

/**
 * @brief Calcola l'insieme di comparabilità di un singolo elemento.
 * * @details
 * In un insieme parzialmente ordinato (POSet), due elementi x ed e sono comparabili
 * se x ≤ e oppure e ≤ x. L'insieme di comparabilità di e è quindi l'unione
 * del suo down-set principale (↓e) e del suo up-set principale (↑e).
 * * Sfruttando l'architettura CSR e i BitSet precalcolati, l'operazione si riduce
 * a una singola operazione hardware SIMD (OR bit a bit) tra le due maschere:
 * ComparabilitySet(e) = ↓e ∪ ↑e.
 * Complessità abbattuta da O(N * log N) a O(N/64) cicli di clock.
 * * @cite B. A. Davey, H. A. Priestley, "Introduction to Lattices and Order" (2nd Ed.),
 * Cambridge University Press, 2002.
 */
BitSet POSet::ComparabilitySetOf(std::uint64_t e) const {
    const std::uint64_t n = size();
    
    // Se l'elemento è fuori dai limiti, restituiamo un BitSet vuoto
    if (e >= n) {
        return BitSet(n);
    }
    
    // Inizializziamo il risultato copiando direttamente l'up-set di 'e'
    BitSet result = upset_bits_[e];
    
    // Unione logica (OR bit a bit) con il down-set di 'e'
    result |= downset_bits_[e];
    
    return result; // Nessuna copia extra grazie alla RVO (Return Value Optimization)
}

/**
 * @brief Calcola l'insieme degli elementi incomparabili a un dato elemento.
 * * @details
 * Due elementi x ed e sono incomparabili se non vale né x <= e né e <= x.
 * Matematicamente, l'insieme di incomparabilità di e è il complemento esatto
 * del suo insieme di comparabilità rispetto all'universo degli elementi:
 * Incomp(e) = U \ (↓e ∪ ↑e).
 * * L'algoritmo sfrutta il ComparabilitySetOf precalcolato e ne inverte i bit
 */
BitSet POSet::IncomparabilitySetOf(std::uint64_t e) const {
    const std::uint64_t n = size();
    
    // Se l'elemento è fuori range, restituiamo un set vuoto
    if (e >= n) {
        return BitSet(n);
    }
    
    // 1. Otteniamo l'insieme di comparabilità (↓e ∪ ↑e)
    BitSet result = ComparabilitySetOf(e);
    
    // 2. Invertiamo tutti i bit (calcolo del complemento)
    result.flip();
    
    // Visto che "e" non è incomparabile a se stesso (e <= e è sempre vero),
    // la flip ha già automaticamente gestito l'auto-esclusione di 'e',
    // perché ComparabilitySetOf includeva giustamente 'e'.
    
    return result; // Nessuna copia extra (RVO)
}

/**
 * @brief Calcola l'insieme degli elementi massimali del POSet.
 * * @details
 * In un insieme parzialmente ordinato, un elemento x è massimale se non esiste
 * alcun y tale che x < y. Di conseguenza, l'up-set stretto di x è vuoto.
 * * L'algoritmo sfrutta il conteggio dei bit accesi (popcount).
 * Se la rappresentazione dei bit è stretta (antiriflessiva), l'upset di un massimale ha 0 bit.
 * Se la rappresentazione è larga (riflessiva), l'upset di un massimale ha 1 bit (se stesso).
 * In entrambi i casi, la condizione (count <= 1) individua univocamente i massimali
 * senza diramazioni condizionali complesse.
 * * @cite B. A. Davey, H. A. Priestley, "Introduction to Lattices and Order" (2nd Ed.),
 * Cambridge University Press, 2002.
 */
BitSet POSet::Maximals() const {
    const std::uint64_t n = size();
    BitSet result(n); // Pre-allocato a zero, pronto per contenere i massimali
    
    for (std::uint64_t v = 0; v < n; ++v) {
        // La condizione <= 1 rende l'algoritmo robusto contro entrambe
        // le rappresentazioni (riflessiva o irriflessiva).
        // Il metodo count() esegue la somma tramite popcount hardware SIMD.
        if (upset_bits_[v].count() <= 1) {
            result.set(v);
        }
    }
    
    return result; // Nessuna copia extra (RVO)
}

/**
 * @brief Calcola l'insieme degli elementi minimali del POSet.
 * @details Un elemento è minimale se il suo down-set non contiene elementi
 * oltre (al massimo) a se stesso. Duale esatto della funzione Maximals().
 */
BitSet POSet::Minimals() const {
    const std::uint64_t n = size();
    BitSet result(n);
    
    for (std::uint64_t v = 0; v < n; ++v) {
        // Stessa logica blindata: <= 1 bit acceso
        if (downset_bits_[v].count() <= 1) {
            result.set(v);
        }
    }
    
    return result;
}

/**
 * @brief Verifica se un singolo elemento del POSet è massimale.
 * * @details
 * Un elemento 'e' è massimale se il suo up-set stretto è vuoto.
 * Come per il calcolo globale dei massimali, l'uso della condizione (count <= 1)
 * garantisce un'esecuzione sicura e corretta a prescindere che la matrice di bit
 * interna sia definita come riflessiva o irriflessiva.
 */
bool POSet::IsMaximal(std::uint64_t e) const {
    // Controllo di sicurezza: se l'elemento non esiste, non può essere massimale
    if (e >= size()) {
        return false;
    }
    
    // Un massimale non ha successori oltre (al massimo) a se stesso.
    return upset_bits_[e].count() <= 1;
}

/**
 * @brief Verifica se un singolo elemento del POSet è minimale.
 * @details Un elemento è minimale se il suo down-set non ha elementi
 * oltre a (al massimo) se stesso. Speculare esatto di IsMaximal().
 */
bool POSet::IsMinimal(std::uint64_t e) const {
    if (e >= size()) {
        return false;
    }
    
    // Un minimale non ha predecessori oltre (al massimo) a se stesso.
    return downset_bits_[e].count() <= 1;
}

/**
 * @brief Calcola la Relazione di Copertura (Cover Relation) del POSet.
 * * @details
 * In un insieme parzialmente ordinato, v1 "copre" v2 (v1 ≺ v2) se v1 < v2
 * e non esiste alcun terzo elemento v3 tale che v1 < v3 < v2.
 * L'insieme delle relazioni di copertura costituisce gli archi del Diagramma di Hasse.
 * * L'algoritmo sfrutta l'architettura CSR e le maschere di bit per trovare
 * i successori immediati in O(1) controlli incrociati. L'utilizzo di std::vector
 * al posto di std::list previene la frammentazione della cache e velocizza l'iterazione.
 * * @cite B. A. Davey, H. A. Priestley, "Introduction to Lattices and Order" (2nd Ed.),
 * Cambridge University Press, 2002. (Capitolo 1: Hasse diagrams).
 */
std::vector<std::pair<std::uint64_t, std::uint64_t>> POSet::CoverRelation() const {
    const std::uint64_t n = size();
    
    // std::vector offre un accesso e un'allocazione enormemente più rapida di std::list
    std::vector<std::pair<std::uint64_t, std::uint64_t>> result;
    
    // Riservare memoria in anticipo aiuta a prevenire riallocazioni dinamiche.
    // Una buona stima per i diagrammi di Hasse sparsi è 2 o 3 archi in media per nodo.
    result.reserve(n * 2);
    
    for (std::uint64_t v1 = 0; v1 < n; ++v1) {
        // Estraiamo il range dei successori di v1 usando gli array CSR
        std::uint64_t start = up_row_ptr[v1];
        std::uint64_t end = up_row_ptr[v1 + 1];
        
        // Iteriamo su tutti i candidati successori 'v2'
        for (std::uint64_t i = start; i < end; ++i) {
            std::uint64_t v2 = up_col_idx[i];
            
            // Un nodo non può coprire se stesso (la cover relation è stretta/irriflessiva)
            if (v1 == v2) continue;
            
            bool is_immediate = true;
            
            // Verifichiamo l'esistenza di un elemento intermedio 'v3'
            for (std::uint64_t j = start; j < end; ++j) {
                std::uint64_t v3 = up_col_idx[j];
                
                if (v3 == v1 || v3 == v2) continue;
                
                // Se v3 >= v1 (già garantito dal ciclo) e v2 >= v3, allora v3 è intermedio.
                // Il controllo è O(1) assoluto grazie ai BitSet precalcolati.
                if (upset_bits_[v3].test(v2)) {
                    is_immediate = false;
                    break;
                }
            }
            
            // Se nessun v3 si interpone tra v1 e v2, allora v1 copre v2.
            if (is_immediate) {
                // emplace_back costruisce la std::pair direttamente in loco, zero copie temporanee
                result.emplace_back(v1, v2);
            }
        }
    }
    
    // Grazie alla Return Value Optimization (RVO), questo vettore non viene mai copiato
    // all'uscita dalla funzione, il suo blocco di memoria viene semplicemente trasferito.
    return result;
}

/**
 * @brief Calcola tutte le coppie uniche di elementi incomparabili nel POSet.
 * * @details
 * Due elementi v1 e v2 sono incomparabili se non vale né v1 <= v2 né v2 <= v1.
 * La funzione restituisce solo combinazioni uniche (v1, v2) con v1 < v2 per evitare duplicati.
 * L'uso dei BitSet precalcolati abbatte i tempi logaritmici di ricerca:
 * interrogare le maschere upset_bits_ e downset_bits_ ha un costo O(1) puro,
 * permettendo l'esplorazione dell'intero spazio (N*(N-1)/2) alla massima velocità hardware.
 */
std::vector<std::pair<std::uint64_t, std::uint64_t>> POSet::Incomparabilities() const {
    const std::uint64_t n = size();
    
    std::vector<std::pair<std::uint64_t, std::uint64_t>> result;
    
    // Riserviamo un blocco di memoria in anticipo per mitigare le riallocazioni dinamiche
    // su POSet sparsi dove l'incomparabilità è statisticamente elevata.
    result.reserve(n * 2);
    
    for (std::uint64_t v1 = 0; v1 < n; ++v1) {
        // Iteriamo da v1 + 1 per evitare di testare (v1, v1) e per
        // ignorare combinazioni già testate (es. se ho testato (1, 2) ignoro (2, 1)).
        for (std::uint64_t v2 = v1 + 1; v2 < n; ++v2) {
            
            // v1 e v2 sono incomparabili se:
            // 1. v2 NON è un successore di v1 (v1 non <= v2)
            // 2. v2 NON è un predecessore di v1 (v2 non <= v1)
            if (!upset_bits_[v1].test(v2) && !downset_bits_[v1].test(v2)) {
                // emplace_back costruisce la coppia direttamente nella memoria del vector
                result.emplace_back(v1, v2);
            }
        }
    }
    
    return result; // Nessuna copia extra (RVO)
}

/**
 * @brief Verifica se il POSet corrente è un'estensione d'ordine del POSet parametrizzato.
 * * @details
 * Un POSet P2 (this) è un'estensione di P1 (p) se e solo se hanno lo stesso insieme
 * di elementi e ogni relazione in P1 esiste anche in P2.
 * Matematicamente, questo significa che per ogni nodo x: ↑x (in P1) ⊆ ↑x (in P2).
 * * L'algoritmo confronta i nomi/ID e poi delega la verifica insiemistica a
 * istruzioni bit a bit hardware (is_subset_of), eliminando migliaia di ricerche su alberi red-black.
 * * @cite B. A. Davey, H. A. Priestley, "Introduction to Lattices and Order" (2nd Ed.).
 */
bool POSet::IsExtensionOf(const POSet& p) const {
    // 1. Controllo base sulle dimensioni
    if (ename_vs_eid.size() != p.ename_vs_eid.size()) {
        return false;
    }
    
    // 2. Controllo identità dei nodi (mappatura nomi).
    // Usiamo std::equal con una lambda per verificare l'uguaglianza delle sole chiavi (first).
    bool same_names = std::equal(
                                 ename_vs_eid.begin(), ename_vs_eid.end(),
                                 p.ename_vs_eid.begin(),
                                 [](const auto& a, const auto& b) {
                                     return a.first == b.first;
                                 }
                                 );
    
    if (!same_names) {
        return false;
    }
    
    // 3. Controllo dell'estensione d'ordine: \leq_p ⊆ \leq_this
    const std::uint64_t n = size();
    for (std::uint64_t v = 0; v < n; ++v) {
        // Se l'upset nel poset parametrizzato NON è sottoinsieme
        // dell'upset del nostro poset, allora non siamo un'estensione.
        if (!p.upset_bits_[v].IsSubsetOf(upset_bits_[v])) {
            return false;
        }
    }
    
    return true;
}

/**
 * @brief Estrae tutte le coppie di elementi comparabili (Relazione d'Ordine).
 * * @details
 * Genera l'insieme di tutte le coppie (x, y) tali che x <= y.
 * Sfruttando la memorizzazione CSR (Compressed Sparse Row) degli upset,
 * il numero totale di relazioni è noto a priori. Questo permette di pre-allocare
 * esattamente la memoria necessaria, riducendo la complessità a O(N + E)
 * con zero riallocazioni (capacity == size) e zero "cache miss".
 */
std::vector<std::pair<std::uint64_t, std::uint64_t>> POSet::Comparabilities() const {
    const std::uint64_t n = size();
    std::vector<std::pair<std::uint64_t, std::uint64_t>> result;
    
    // Sicurezza: se il POSet è vuoto, ritorniamo subito
    if (n == 0) return result;
    
    // In un formato CSR, l'ultimo elemento di row_ptr contiene esattamente
    // il numero totale di elementi (archi) presenti in col_idx.
    // Pre-allochiamo la memoria esatta per evitare riallocazioni durante l'inserimento!
    std::uint64_t total_edges = up_row_ptr[n];
    result.reserve(total_edges);
    
    // Iterazione lineare rapidissima sugli array contigui
    for (std::uint64_t out = 0; out < n; ++out) {
        std::uint64_t start = up_row_ptr[out];
        std::uint64_t end = up_row_ptr[out + 1];
        
        for (std::uint64_t i = start; i < end; ++i) {
            // emplace_back costruisce la std::pair direttamente nel blocco di memoria
            result.emplace_back(out, up_col_idx[i]);
        }
    }
    
    return result; // Return Value Optimization (RVO)
}

/**
 * @brief Crea una copia esatta e indipendente (clone) del POSet corrente.
 * @details Sfrutta il costruttore di copia di default generato dal compilatore C++.
 */
std::shared_ptr<POSet> POSet::Clone() const {
    // Utilizziamo "new" direttamente qui dentro perché, essendo un metodo
    // della classe stessa, ha pieno accesso al costruttore di copia privato.
    // L'oggetto viene clonato istantaneamente in memoria senza ri-eseguire Build().
    return std::shared_ptr<POSet>(new POSet(*this));
}

/**
 * @brief Verifica se il POSet è un Ordine Totale (o Lineare).
 * @details
 * In un ordine totale (o catena), non esistono elementi incomparabili.
 * Di conseguenza, ogni nodo deve essere comparabile con tutti gli altri N - 1 nodi.
 * * Poiché il POSet non ha cicli, l'up-set stretto e il down-set stretto di un nodo
 * sono disgiunti. Sfruttando le maschere di bit, calcoliamo la somma delle loro
 * cardinalità.
 * - Se la matrice è riflessiva, la somma (up + down) è N + 1 (il nodo è contato 2 volte).
 * - Se è irriflessiva, la somma è N - 1 (il nodo non è contato).
 * L'algoritmo short-circuita e ritorna 'false' al primo nodo che non rispetta
 * questa condizione.
 */
bool POSet::IsTotalOrder() const {
    const std::uint64_t n = size();
    
    // Un POSet vuoto o con un solo elemento è, per definizione, un ordine totale
    if (n <= 1) {
        return true;
    }
    
    for (std::uint64_t v = 0; v < n; ++v) {
        // popcount() hardware: calcola quanti bit a 1 ci sono
        std::uint64_t up_count = upset_bits_[v].count();
        std::uint64_t down_count = downset_bits_[v].count();
        
        std::uint64_t related = up_count + down_count;
        
        // Se la somma dei comparabili non coincide con quella attesa
        // per un ordine totale, significa che ci sono "buchi" (incomparabilità).
        if (related != (n + 1) && related != (n - 1)) {
            return false;
        }
    }
    
    return true;
}


/**
 * @brief Aggiorna la matrice dei risultati sommando i valori calcolati da una FunctionLinearExtension.
 * @param fle Funzione di cui estrarre e sommare i risultati (sola lettura).
 */
void POSet::SumUpdate(Tensor<double, 2>& ris, const FunctionLinearExtension& fle) noexcept {
    const std::uint64_t n = fle.ResSize();

    for (std::uint64_t k = 0; k < n; ++k) {
        const std::uint64_t first  = fle.At0(k);
        const std::uint64_t second = fle.At1(k);
        const double add_val       = fle.At2(k);
        ris(first, second) += add_val;
    }
}

// Instanziazione esplicita: garantisce che il motore parallelo venga compilato
// (e devirtualizzato sul tipo concreto) anche in assenza di call site nel core.
template void POSet::evaluation_parallel<ECSBubleyDyerFromLEMatrix>(
                      ECSBubleyDyerFromLEMatrix&,
                      std::vector<std::unique_ptr<Tensor<double, 2>>>&,
                      std::uint64_t&, bool&, DisplayMessage*,
                      POSet::EvaluationUpdateStrategy,
                      std::uint64_t);


/**
 * @brief Computes the probabilistic dominance matrix using the LPOM (Relative Ideals) model.
 * * @details This implementation uses the original Local Partial Order Model formula
 * (Equation 13), developed by Brüggemann, Lerche, and Sørensen.
 * For pairs of incomparable elements, the mutual probability is estimated by evaluating
 * the cardinalities of the relative ideals, meaning the set difference between the
 * Up-sets and Down-sets of the two elements: |U(x) \ U(y)| and |O(x) \ O(y)|.
 * * @see NERI Technical Report No. 479 (2004) - Equation 13.
 * @see Brüggemann, R., Sørensen, P. B., Lerche, D. B., Carlsen, L. (2004).
 * "Estimation of Averaged Ranks by a Local Partial Order Model".
 * * @return Tensor<double, 2> A dense square matrix where element (i, j) represents
 * the probabilistic estimate P(i > j) that element 'i' dominates element 'j'.
 */
Tensor<double, 2> POSet::BLSDominanceRelative() const {
    const std::uint64_t n = size();
    
    // Inizializza una matrice quadrata N x N
    // Usa l'inizializzazione standard che formatta tutto a 0.0
    Tensor<double, 2> matrice({n, n}, 0.0);
    
    // Inizializza la diagonale principale a 1.0 (auto-dominanza)
    for (std::uint64_t k = 0; k < n; ++k) {
        matrice(k, k) = 1.0;
    }
    
    for (std::uint64_t v1 = 0; v1 < n; ++v1) {
        // Cache della riga per evitare letture ripetute
        const BitSet& up_v1 = upset_bits_[v1];
        const BitSet& down_v1 = downset_bits_[v1];
        // Flag per sapere se l'identità è inclusa nel BitSet
        const std::uint64_t self_up_v1 = up_v1.test(v1) ? 1 : 0;
        const std::uint64_t self_down_v1 = down_v1.test(v1) ? 1 : 0;
        
        for (std::uint64_t v2 = v1 + 1; v2 < n; ++v2) {
            
            // Relazioni deterministiche (già comparabili)
            if (up_v1.test(v2)) {
                matrice(v1, v2) = 1.0;
                // matrice(v2, v1) è già 0.0 per l'inizializzazione
            }
            else if (upset_bits_[v2].test(v1)) {
                matrice(v2, v1) = 1.0;
                // matrice(v1, v2) è già 0.0
            }
            // Incomparabili: calcolo probabilistico BLS
            else {
                const BitSet& up_v2 = upset_bits_[v2];
                const BitSet& down_v2 = downset_bits_[v2];
                
                const std::uint64_t self_up_v2 = up_v2.test(v2) ? 1 : 0;
                const std::uint64_t self_down_v2 = down_v2.test(v2) ? 1 : 0;
                
                // upij = | up(v1) \ up(v2) \ {v1, v2} |
                // Sappiamo che v2 non è in up(v1) (incomparabili), quindi togliamo solo v1 se presente.
                const std::uint64_t upij = up_v1.CountDifference(up_v2) - self_up_v1;
                const std::uint64_t upji = up_v2.CountDifference(up_v1) - self_up_v2;
                
                // downij = | down(v1) \ down(v2) \ {v1, v2} |
                const std::uint64_t downij = down_v1.CountDifference(down_v2) - self_down_v1;
                const std::uint64_t downji = down_v2.CountDifference(down_v1) - self_down_v2;
                
                // Formule originali BLS
                const double aij = (static_cast<double>(upij) + 1.0) / (static_cast<double>(downij) + 1.0);
                const double aji = (static_cast<double>(upji) + 1.0) / (static_cast<double>(downji) + 1.0);
                
                const double sum_a = aij + aji;
                
                // Assegnazione simmetrica incrociata
                matrice(v2, v1) = aji / sum_a;
                matrice(v1, v2) = aij / sum_a;
            }
        }
    }
    
    return matrice;
}

/**
 * @brief Computes the probabilistic dominance matrix using the LPOM (Absolute Ideals) model.
 * * @details Implements Equation 13' (the Q'-based variant) from the theory by Brüggemann et al.
 * Unlike the relative variant, this function estimates the degrees of freedom of an element
 * based exclusively on the cardinality of its absolute ideals: |U(x)| and |O(x)|.
  * * @see NERI Technical Report No. 479 (2004) - Equation 13'.
 * * @return Tensor<double, 2> A dense square matrix where element (i, j) represents
 * the probabilistic estimate P(i > j) calculated using absolute ideals.
 */
[[nodiscard]] Tensor<double, 2> POSet::BLSDominanceAbsolute() const {
    const std::uint64_t n = size();
    
    // Inizializza matrice azzerata
    Tensor<double, 2> matrice({n, n}, 0.0);
    
    // Diagonale principale (auto-dominanza)
    for (std::uint64_t k = 0; k < n; ++k) {
        matrice(k, k) = 1.0;
    }
    
    for (std::uint64_t v1 = 0; v1 < n; ++v1) {
        const BitSet& up_v1 = upset_bits_[v1];
        const BitSet& down_v1 = downset_bits_[v1];
        
        const std::uint64_t self_up_v1 = up_v1.test(v1) ? 1 : 0;
        const std::uint64_t self_down_v1 = down_v1.test(v1) ? 1 : 0;
        
        // Cardinalità assolute per v1 (escluso se stesso)
        const std::uint64_t up_size_v1 = up_v1.count() - self_up_v1;
        const std::uint64_t down_size_v1 = down_v1.count() - self_down_v1;
        
        for (std::uint64_t v2 = v1 + 1; v2 < n; ++v2) {
            
            // Relazioni deterministiche
            if (up_v1.test(v2)) {
                matrice(v1, v2) = 1.0;
            }
            else if (upset_bits_[v2].test(v1)) {
                matrice(v2, v1) = 1.0;
            }
            // Incomparabili: Equazione 13'
            else {
                const BitSet& up_v2 = upset_bits_[v2];
                const BitSet& down_v2 = downset_bits_[v2];
                
                const std::uint64_t self_up_v2 = up_v2.test(v2) ? 1 : 0;
                const std::uint64_t self_down_v2 = down_v2.test(v2) ? 1 : 0;
                
                // Cardinalità assolute per v2
                const std::uint64_t up_size_v2 = up_v2.count() - self_up_v2;
                const std::uint64_t down_size_v2 = down_v2.count() - self_down_v2;
                
                // Variabile Q' (Eq 13') basata sugli ideali assoluti
                const double q_v1 = (static_cast<double>(up_size_v1) + 1.0) / (static_cast<double>(down_size_v1) + 1.0);
                const double q_v2 = (static_cast<double>(up_size_v2) + 1.0) / (static_cast<double>(down_size_v2) + 1.0);
                
                const double sum_q = q_v1 + q_v2;
                
                // Assegnazione simmetrica (P(v1 > v2) = q_v2 / (q_v1 + q_v2))
                // riga dominata dalla colonna
                matrice(v2, v1) = q_v2 / sum_q;
                matrice(v1, v2) = q_v1 / sum_q;
            }
        }
    }
    
    return matrice;
}



/**
 * @brief Calcola la Somma Lineare (Ordinal Sum) di due POSet disgiunti (P1 ⊕ P2).
 * @details
 * In una somma lineare, ogni elemento di P1 precede ogni elemento di P2.
 * L'algoritmo verifica la disgiunzione dei due insiemi in O(N) tramite hash map,
 * fonde gli elementi e costruisce gli archi sommando:
 * 1. La relazione d'ordine interna di P1
 * 2. La relazione d'ordine interna di P2
 * 3. Il prodotto cartesiano P1 x P2 (tutti i nodi di P1 verso tutti i nodi di P2)
 * Tutti i vettori sono rigorosamente pre-allocati per azzerare la frammentazione.
 */
std::unique_ptr<POSet> POSet::LinearSum(const POSet& p1, const POSet& p2) {
    const std::uint64_t n1 = p1.size();
    const std::uint64_t n2 = p2.size();
    const std::uint64_t total_nodes = n1 + n2;
    
    std::vector<std::string> elements;
    elements.reserve(total_nodes);
    
    std::unordered_set<std::string> elements_set;
    elements_set.reserve(total_nodes);
    
    
    
    for (std::uint64_t i = 0; i < n1; ++i) {
        std::string_view name = p1.GetElementName(i);
        elements_set.insert(std::string(name));
        elements.emplace_back(name);
    }
    
    for (std::uint64_t i = 0; i < n2; ++i) {
        std::string_view name = p2.GetElementName(i);
        
        if (elements_set.find(std::string(name)) != elements_set.end()) {
            throw MyException(std::format("POSets are not disjoint! Overlapping element: {}", name));
        }
        
        elements_set.insert(std::string(name));
        elements.emplace_back(name);
    }
    
    auto comp1 = p1.Comparabilities();
    auto comp2 = p2.Comparabilities();
    
    // Vettore puro pre-allocato. Velocità di inserimento: massima!
    std::vector<std::pair<std::string, std::string>> edges;
    std::uint64_t total_edges = comp1.size() + comp2.size() + (n1 * n2);
    edges.reserve(total_edges);
    
    for (const auto& pair : comp1) {
        edges.emplace_back(p1.GetElementName(pair.first), p1.GetElementName(pair.second));
    }
    
    for (std::uint64_t i = 0; i < n1; ++i) {
        std::string_view name1 = p1.GetElementName(i);
        
        for (std::uint64_t j = 0; j < n2; ++j) {
            edges.emplace_back(name1, p2.GetElementName(j));
        }
    }
    
    for (const auto& pair : comp2) {
        edges.emplace_back(p2.GetElementName(pair.first), p2.GetElementName(pair.second));
    }
    return POSet::Build(elements, edges, true);
}

/**
 * @brief Calcola la Somma Disgiunta (Parallel Composition) di due POSet (P1 + P2).
 * @details
 * I due POSet vengono uniti mantenendo inalterate le loro relazioni d'ordine interne,
 * senza creare alcuna relazione di comparabilità tra gli elementi di P1 e quelli di P2.
 * L'algoritmo sfrutta l'hash set per la verifica della disgiunzione in O(N) e
 * vettori pre-allocati per massimizzare il throughput di memoria.
 */
std::unique_ptr<POSet> POSet::DisjointSum(const POSet& p1, const POSet& p2) {
    const std::uint64_t n1 = p1.size();
    const std::uint64_t n2 = p2.size();
    const std::uint64_t total_nodes = n1 + n2;
    
    // 1. Unione dei nodi e check disgiunzione in O(N)
    std::vector<std::string> elements;
    elements.reserve(total_nodes);
    
    std::unordered_set<std::string> elements_set;
    elements_set.reserve(total_nodes);
    
    // Nodi di P1
    for (std::uint64_t i = 0; i < n1; ++i) {
        std::string_view name = p1.GetElementName(i);
        elements_set.insert(std::string(name));
        elements.emplace_back(name);
    }
    
    // Nodi di P2 con controllo sovrapposizione
    for (std::uint64_t i = 0; i < n2; ++i) {
        std::string_view name = p2.GetElementName(i);
        
        auto [it, inserted] = elements_set.insert(std::string(name));
        
        if (!inserted) {
            throw MyException(std::format("POSets are not disjoint! Overlapping element: {}", name));
        }
        
        elements.emplace_back(name);
    }
    
    
    // 2. Estrazione istantanea degli archi
    auto comp1 = p1.Comparabilities();
    auto comp2 = p2.Comparabilities();
    
    // 3. Vettore puro pre-allocato esattamente al byte
    std::vector<std::pair<std::string, std::string>> edges;
    edges.reserve(comp1.size() + comp2.size());
    
    // A. Aggiungiamo gli archi interni di P1
    for (const auto& pair : comp1) {
        // Usa push_back({a, b}) se il tuo compilatore fa capricci,
        // altrimenti emplace_back è perfetto.
        edges.emplace_back(p1.GetElementName(pair.first), p1.GetElementName(pair.second));
    }
    
    // B. Aggiungiamo gli archi interni di P2
    // (Nessun ciclo incrociato necessario per la somma disgiunta!)
    for (const auto& pair : comp2) {
        edges.emplace_back(p2.GetElementName(pair.first), p2.GetElementName(pair.second));
    }
    
    // 4. Chiamata al nuovo Factory HPC
    return POSet::Build(elements, edges, true);
}

/**
 * @brief Esegue l'operazione di Lifting aggiungendo un nuovo elemento "bottom".
 * @details
 * Crea un nuovo POSet in cui new_element è strettamente minore di tutti gli elementi
 * precedentemente esistenti.
 * L'algoritmo sfrutta l'hash map interna per validare l'univocità in O(1)
 * e alloca vettori di dimensione esatta per bypassare la riallocazione dinamica.
 */
std::unique_ptr<POSet> POSet::Lifting(const std::string& new_element) const {
    const std::uint64_t n = size();
    
    // 1. Controllo duplicati fulmineo in O(1) sfruttando la Hash Map interna
    // (Assumo che la tua mappa stringa->ID si chiami ename_vs_eid.
    // Altrimenti adegua il nome a quello che usi nel tuo private!)
    if (ename_vs_eid.find(new_element) != ename_vs_eid.end()) {
        throw MyException("POSet bottom element is not new!");
    }
    
    // 2. Pre-allocazione esatta dei Nodi
    std::vector<std::string> elements;
    elements.reserve(n + 1);
    
    // Inseriamo prima i nodi esistenti
    for (std::uint64_t i = 0; i < n; ++i) {
        elements.emplace_back(GetElementName(i));
    }
    
    // Aggiungiamo il nuovo elemento bottom
    elements.push_back(new_element);
    
    // 3. Estrazione immediata degli archi esistenti in RAM sequenziale
    auto comp = Comparabilities();
    
    // 4. Pre-allocazione esatta degli Archi
    std::vector<std::pair<std::string, std::string>> edges;
    edges.reserve(comp.size() + n);
    
    // A. Il nuovo elemento (bottom) precede TUTTI gli elementi preesistenti
    for (std::uint64_t i = 0; i < n; ++i) {
        edges.emplace_back(new_element, GetElementName(i));
        // (Se Clang fa capricci, usa: edges.push_back({new_element, GetElementName(i)}); )
    }
    
    // B. Copiamo le relazioni interne preesistenti
    for (const auto& pair : comp) {
        edges.emplace_back(GetElementName(pair.first), GetElementName(pair.second));
    }
    
    // 5. Deleghiamo la costruzione al Factory HPC
    return POSet::Build(elements, edges, true);
}


/**
 * @brief Calcola il POSet Duale (tutte le relazioni d'ordine invertite).
 * @details
 * Gli elementi rimangono invariati, mentre per ogni arco (x, y) nel POSet
 * originale, viene creato un arco (y, x) nel Duale.
 * Utilizza la rappresentazione piatta CSR (Comparabilities) per invertire
 * il grafo alla massima velocità di bus della memoria.
 */
std::unique_ptr<POSet> POSet::Dual() const {
    const std::uint64_t n = size();
    
    // 1. Estrazione lineare degli elementi in vettore pre-allocato
    std::vector<std::string> elements;
    elements.reserve(n);
    
    for (std::uint64_t i = 0; i < n; ++i) {
        elements.emplace_back(GetElementName(i));
    }
    
    // 2. Estrazione degli archi in RAM sequenziale tramite la tua funzione ottimizzata
    auto comp = Comparabilities();
    
    // 3. Pre-allocazione rigorosa degli archi invertiti
    std::vector<std::pair<std::string, std::string>> edges;
    edges.reserve(comp.size());
    
    // 4. Inversione matematica della relazione d'ordine (y, x)
    for (const auto& pair : comp) {
        // Leggiamo second prima di first!
        edges.emplace_back(GetElementName(pair.second), GetElementName(pair.first));
        
        // (Se il compilatore non digerisce emplace_back, usa:)
        // edges.push_back({GetElementName(pair.second), GetElementName(pair.first)});
    }
    
    // 5. Costruzione tramite il Factory
    // NOTA: il duale di un grafo chiuso transitivamente è a sua volta chiuso transitivamente quindi passiamo false.
    return POSet::Build(elements, edges, false);
}



/**
 * @brief Calcola l'intersezione di due POSET basandosi sui nomi degli elementi.
 * @details La relazione a <= b esiste nel risultato solo se esiste in p1 E in p2.
 */
std::unique_ptr<POSet> POSet::Intersection(const POSet& p1, const POSet& p2) {
    // 1. Validazione: devono avere la stessa cardinalità
    if (p1.size() != p2.size()) {
        throw std::runtime_error("Intersection error: POSets must have the same number of elements.");
    }
    
    const auto& p1_names = p1.eid_vs_ename;  // Vettore nomi p1
    const auto& p2_map   = p2.ename_vs_eid; // Mappa nomi->ID p2
    const auto& p1_ups   = p1.UpSets();      // DATASTORE (vector<BitSet>)
    const auto& p2_ups   = p2.UpSets();      // DATASTORE (vector<BitSet>)
    
    std::vector<std::pair<std::string, std::string>> common_relations;
    // Stima della capacità per evitare riallocazioni
    common_relations.reserve(p1.size());
    
    // 2. Iterazione sugli elementi di p1
    for (std::uint64_t i1 = 0; i1 < p1_names.size(); ++i1) {
        const std::string& name_i = p1_names[i1];
        
        // Trova l'ID corrispondente in p2 tramite il nome
        auto it_i2 = p2_map.find(name_i);
        if (it_i2 == p2_map.end()) {
            throw std::runtime_error("Intersection error: element '" + name_i + "' missing in second POSet.");
        }
        std::uint64_t i2 = it_i2->second;
        
        // Recuperiamo i BitSet degli upset per i due elementi corrispondenti
        const BitSet& upset1 = p1_ups[i1];
        const BitSet& upset2 = p2_ups[i2];
        
        // 3. Confronto delle relazioni
        // Verifichiamo quali elementi 'j' sono maggiori di 'i' in entrambi i POSET
        for (std::uint64_t j1 = 0; j1 < p1_names.size(); ++j1) {
            // Se j1 è nell'upset di i1 in p1...
            if (upset1.test(j1)) {
                const std::string& name_j = p1_names[j1];
                
                // ...controlliamo se name_j è nell'upset di i2 in p2
                auto it_j2 = p2_map.find(name_j);
                if (it_j2 != p2_map.end()) {
                    std::uint64_t j2 = it_j2->second;
                    if (upset2.test(j2)) {
                        // La relazione esiste in entrambi!
                        common_relations.emplace_back(name_i, name_j);
                    }
                }
            }
        }
    }
    
    // 4. Costruzione del nuovo POSET
    return POSet::Build(p1_names, common_relations, true);
}

[[nodiscard]] std::unique_ptr<LinearExtensionGenerator> POSet::CreateLinearExtensionGenerator() {
    LatticeOfIdeals* lattice_ptr = GetLatticeOfIdeals();
    
    if (lattice_ptr == nullptr) {
        throw MyException("Error: Cannot create LEGTreeOfIdeals. LatticeOfIdeals was not generated in the POSet.");
    }
    
    return std::make_unique<LEGTreeOfIdeals>(size(), *lattice_ptr);
}

/**
 * @brief Calcolo della Mutual Ranking Probability (MRP) matrix.
 * @details Metodo generico inserito in POSet.cpp. Utilizza il Factory Method
 * per istanziare il generatore specifico. Non usa le stringhe dei nomi
 * per garantire sicurezza HPC sulla memoria.
 */
[[nodiscard]] Tensor<double, 2> POSet::ComputeMRP() {
    const std::uint64_t n = this->size();
    
    // 1. Guardia sulla memoria (Agnostica)
    if (n > 45000) {
        throw MyException(std::format("Dimensioni POSet ({}) troppo grandi per una matrice MRP densa.", n));
    }
    
    // 2. Preparazione delle Funzioni (FLE)
    std::vector<std::unique_ptr<FunctionLinearExtension>> fles;
    fles.push_back(std::make_unique<FLEMutualRankingProbability>(this));
    
    // 3. Preparazione Risultati
    std::vector<std::unique_ptr<Tensor<double, 2>>> eval_results;
    eval_results.push_back(std::make_unique<Tensor<double, 2>>(std::array<std::uint64_t, 2>{n, n}, 0.0));
    
    // 4. Creazione del Generatore tramite Factory Method
    std::unique_ptr<LinearExtensionGenerator> leg = this->CreateLinearExtensionGenerator();
    
    
    leg->Start(0);
    
    // 5. Variabili di stato
    std::uint64_t le_count = 0;
    bool end_process = false;
    DisplayMessageNull display_message;
    
    // 6. Esecuzione Engine di Valutazione
    POSet::evaluation(
                     fles,
                     *leg,
                     eval_results,
                     le_count,
                     end_process,
                     &display_message,
                     EvaluationUpdateStrategy::Average
                     );
    
    return std::move(*(eval_results[0]));
}



/*
 

 // ***********************************************
 // ***********************************************
 // ***********************************************
 
 std::shared_ptr<BucketPOSet> POSet::BucketPOSetFromMinimal() const {
 std::set<std::uint_fast64_t> not_used;
 std::set<std::uint_fast64_t> used;
 std::vector<std::shared_ptr<std::set<std::uint_fast64_t>>> buckets;
 auto elements = std::make_shared<std::vector<std::string>>();
 std::list<std::pair<std::string, std::string>> comparabilities;
 
 for(std::uint_fast64_t it = 0; it < elementi->size(); ++it) {
 not_used.insert(it);
 elements->push_back(eid_vs_ename->at(it));
 }
 
 while (not_used.size() > 0) {
 auto current_bucket = std::make_shared<std::set<std::uint_fast64_t>>();
 for (auto check_el : not_used) {
 bool not_minimal = false;
 for(std::uint_fast64_t poset_el = 0; poset_el < elementi->size(); ++poset_el) {
 if (poset_el != check_el && used.find(poset_el) == used.end() && elementi->at(poset_el)->find(check_el) != elementi->at(poset_el)->end()) {
 not_minimal = true;
 break;
 }
 }
 if (!not_minimal) {
 current_bucket->insert(check_el);
 }
 }
 used.insert(current_bucket->begin(), current_bucket->end());
 buckets.push_back(current_bucket);
 for (auto e : *current_bucket) {
 not_used.erase(e);
 }
 }
 
 for (std::uint_fast64_t k = 1; k < buckets.size(); ++k) {
 for (auto down_el_idx : *(buckets.at(k-1))) {
 for (auto up_el_idx : *(buckets.at(k))) {
 auto down_el = eid_vs_ename->at(down_el_idx);
 auto up_el = eid_vs_ename->at(up_el_idx);
 comparabilities.push_back(std::make_pair(down_el, up_el));
 }
 }
 }
 
 auto result = BucketPOSet::Build(elements, comparabilities);
 return result;
 }
 
 // ***********************************************
 // ***********************************************
 // ***********************************************
 
 std::shared_ptr<BucketPOSet> POSet::BucketPOSetFromMaximal() const {
 auto remove = [](POSet::DATASTORE& store, std::set<std::uint_fast64_t>& elements) {
 bool ancora_dati = true;
 for (auto element_to_be_removed: elements) {
 store.at(element_to_be_removed) = nullptr;
 ancora_dati = false;
 for (auto d : store) {
 if (d != nullptr) {
 d->erase(element_to_be_removed);
 ancora_dati = true;
 }
 }
 if (ancora_dati == false) {
 break;
 }
 }
 return ancora_dati;
 };
 
 POSet::DATASTORE data_store_clone(size(), nullptr);
 for (std::uint_fast64_t d = 0; d < elementi->size(); ++d) {
 auto upset_d = elementi->at(d);
 data_store_clone.at(d) = std::make_shared<std::set<std::uint_fast64_t>>(upset_d->begin(), upset_d->end());
 }
 
 
 
 std::vector<std::shared_ptr<std::set<std::uint_fast64_t>>> buckets;
 std::list<std::pair<std::string, std::string>> comparabilities;
 
 bool ancora_dati = true;
 while (ancora_dati) {
 auto current_bucket = std::make_shared<std::set<std::uint_fast64_t>>();
 for (std::uint_fast64_t check_el = 0; check_el < data_store_clone.size(); ++check_el) {
 if (data_store_clone.at(check_el) != nullptr && data_store_clone.at(check_el)->size() == 0) {
 current_bucket->insert(check_el);
 }
 }
 buckets.insert(buckets.begin(), current_bucket);
 ancora_dati = remove(data_store_clone, *current_bucket);
 }
 
 auto elements = eid_vs_ename;
 
 for (std::uint_fast64_t k = 1; k < buckets.size(); ++k) {
 for (auto down_el_idx : *(buckets.at(k-1))) {
 for (auto up_el_idx : *(buckets.at(k))) {
 auto down_el = eid_vs_ename->at(down_el_idx);
 auto up_el = eid_vs_ename->at(up_el_idx);
 comparabilities.push_back(std::make_pair(down_el, up_el));
 }
 }
 }
 
 auto result = BucketPOSet::Build(elements, comparabilities);
 return result;
 }
 
 // ***********************************************
 // ***********************************************
 // ***********************************************
 
 
 
 */
