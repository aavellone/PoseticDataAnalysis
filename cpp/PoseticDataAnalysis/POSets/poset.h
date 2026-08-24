#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <utility>
#include <cstddef>

// restrict qualifier portabile: __restrict__ (GCC/Clang) vs __restrict (MSVC).
#if defined(_MSC_VER)
#define POSET_RESTRICT __restrict
#else
#define POSET_RESTRICT __restrict__
#endif

#include "bit_set.h"
#include "tensor.h"
#include "tensor.h"
#include "random.h"
#include "linear_extension.h"
#include "tree_of_ideals.h"
#include "lattice_of_ideals.h"
#include "function_linear_extension.h"
#include "linear_extension_generator.h"
#include "display_message.h"
#include "loss_function_mrp.h"
#include "score.h"

#include <vector>
#include <unordered_map>
#include <set>
#include <tuple>
#include <memory>
#include <thread>
#include <atomic>
#include <chrono>
#include <exception>
#include <condition_variable>
#include <concepts>
#include <cstdint>
#include <format>
#include <ranges>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <string>
#include <functional>
#include <algorithm>
#include <map>
#include <mutex>
#include <optional>

/**
 * @brief Concept for types providing linear extension data.
 * * Ensures that the template argument provides all necessary methods
 * to access matrix indices and values during the average update process.
 */
template <typename T>
concept LinearExtensionProvider = requires(const T t, std::uint64_t k) {
    { t.ResSize() } -> std::convertible_to<std::uint64_t>;
    { t.At0(k) }    -> std::convertible_to<std::uint64_t>;
    { t.At1(k) }    -> std::convertible_to<std::uint64_t>;
    { t.At2(k) }    -> std::convertible_to<double>;
};

/**
 * @brief Task auto-contenuto di una catena della valutazione parallela.
 * * Impacchetta tutto ciò che serve a un worker per elaborare una catena:
 * il generatore (con seed e punto di partenza propri, già Start()-ato),
 * i cloni per-catena delle funzioni (stateful) e i tensori risultato locali
 * (azzerati). Nessuna condivisione tra thread durante il calcolo: il task
 * vive sullo stack del worker e viene distrutto a catena conclusa (i tensori
 * risultato vengono prima fusi nel merge incrementale).
 * @tparam GeneratorType Tipo concreto del generatore (devirtualizzazione).
 */
template <typename GeneratorType>
struct EvaluationChainTask {
    std::unique_ptr<GeneratorType> leg;
    std::vector<std::unique_ptr<FunctionLinearExtension>> fles;
    std::vector<std::unique_ptr<Tensor<double, 2>>> results;
};

/**
 * @brief Concept per le sorgenti lazy di catene di valutazione.
 * * Una sorgente produce EvaluationChainTask on demand, una alla volta
 * (stesso pattern dei LinearExtensionGenerator): HasNext() dice se esiste
 * una catena successiva, Next() la costruisce e la consegna. La sorgente
 * decide come generare i punti iniziali (colonne di una matrice, campioni
 * casuali, strategie adattive, ...) senza che il motore parallelo cambi.
 * * @warning HasNext()/Next() NON devono essere thread-safe: evaluation_parallel
 * li invoca sempre sotto lock. Non devono MAI toccare l'API di R (i worker
 * non sono il main thread).
 */
template <typename S>
concept EvaluationChainSource = requires(S s) {
    typename S::GeneratorType;
    { s.HasNext() } -> std::convertible_to<bool>;
    { s.Next() }    -> std::same_as<EvaluationChainTask<typename S::GeneratorType>>;
};


/**
 * ==============================================================================
 * ARCHITETTURA HPC: COMPRESSED SPARSE ROW (CSR)
 * ==============================================================================
 * Questo POSet utilizza una rappresentazione CSR per massimizzare la Data Locality
 * e sfruttare il prefetching della cache L1/L2 della CPU.
 * * L'intero grafo delle relazioni è "schiacciato" in due array contigui:
 * 1. ROW_PTR (es. up_row_ptr): Ha dimensione N + 1.
 * Indica l'indice di partenza (e di fine) dei successori per ogni elemento.
 * 2. COL_IDX (es. up_col_idx): Contiene tutti gli ID dei successori concatenati.
 * Supponiamo di avere N=5 elementi (da 0 a 4).
 * Relazioni (Elementi superiori):
 * - 0 ha come superiori {1, 2, 3}
 * - 1 ha come superiori {2}
 * - 2 ha come superiori {}      (Nessuno)
 * - 3 ha come superiori {}      (Nessuno)
 * - 4 ha come superiori {2, 3}
 * * Come viene memorizzato nella CSR:
 * up_col_idx = [ 1, 2, 3, 2, 2, 3 ]   <-- Tutti i superiori concatenati
 * up_row_ptr = [ 0, 3, 4, 4, 4, 6 ]   <-- Offset di riga
 * * Spiegazione degli offset:
 * - Riga 0: inizia all'indice 0 e finisce al 3 (elementi 1, 2, 3).
 * - Riga 1: inizia all'indice 3 e finisce al 4 (elemento 2).
 * - Riga 2: inizia all'indice 4 e finisce al 4 (VUOTO).
 * - Riga 3: inizia all'indice 4 e finisce al 4 (VUOTO).
 * - Riga 4: inizia all'indice 4 e finisce al 6 (elementi 2, 3).
 * - L'ultimo valore (6) è la dimensione totale di up_col_idx.
 * * Iterazione super-veloce su 'i':
 * for(uint64_t idx = up_row_ptr[i]; idx < up_row_ptr[i+1]; ++idx) {
 * uint64_t superiore = up_col_idx[idx];
 * }
 * ==============================================================================
 */


class POSet {
public:
    using DATASTORE = std::vector<BitSet>;
    using BEC = std::tuple<std::vector<std::vector<std::uint64_t>>, std::vector<std::string>, std::vector<std::pair<std::string, std::string>>>;

    enum class EvaluationUpdateStrategy {
        Sum,
        Average
    };

protected:
    // Hashing trasparente C++20 per permettere lookup con std::string_view
    // all'interno di una mappa che ha std::string_view come chiavi.
    struct StringHash {
        using is_transparent = void;
        
        [[nodiscard]] std::size_t operator()(std::string_view sv) const noexcept {
            return std::hash<std::string_view>{}(sv);
        }
    };

    
    
    POSet() = default;

    // --- Struttura dati CSR per HPC ---
    std::vector<std::uint64_t> up_row_ptr;
    std::vector<std::uint64_t> up_col_idx;
    
    std::vector<std::uint64_t> down_row_ptr;
    std::vector<std::uint64_t> down_col_idx;
    
    // Permettono IsComparable(a, b) in O(1) invece di fare ricerca binaria su CSR
    std::vector<BitSet> upset_bits_;   // upset_bits_[a].test(b) == vero se a <= b
    std::vector<BitSet> downset_bits_; // downset_bits_[b].test(a) == vero se b >= a
    
    
    std::vector<std::string> eid_vs_ename;
    
    std::unordered_map<std::string, std::uint64_t, StringHash, std::equal_to<>> ename_vs_eid;
    
    std::unique_ptr<TreeOfIdeals> tree_of_ideals_ = nullptr;
    std::unique_ptr<LatticeOfIdeals> lattice_of_ideals_ = nullptr;
    
public:
    // Factory methods
    static std::unique_ptr<POSet> Build(const BEC& p, bool transitiveClosure = true);
    static std::unique_ptr<POSet> Build(const std::vector<std::string>& elements,
                                        const std::vector<std::pair<std::string, std::string>>& comparabilities,
                                        bool transitiveClosure = true);
    
    // Costruttore di Copia esplicito
    POSet(const POSet& other);
    // Distruttore di default
    virtual ~POSet() = default;
    
    void FillBaseAttributes(const std::vector<std::string>& elements,
                            const std::vector<std::pair<std::string, std::string>>& comparabilities,
                            Score* score,
                            bool transitiveClosure);
    
    [[nodiscard]] virtual std::uint64_t GetElementId(std::string_view name) const {
        if (auto it = ename_vs_eid.find(name); it != ename_vs_eid.end()) [[likely]] {
            return it->second;
        }
        throw MyException(std::format("POSet Error: EID not found for element '{}'.", name));
    }
    
    [[nodiscard]] virtual std::string_view GetElementName(std::uint64_t idx) const;
    
    [[nodiscard]] virtual std::uint64_t size() const noexcept {
        return eid_vs_ename.size();
    }
    [[nodiscard]] virtual std::string to_string(char delimiter = ';') const;

    [[nodiscard]] virtual BitSet GetImmediatePredecessors(std::uint64_t el) const;
    [[nodiscard]] virtual POSet::DATASTORE GetImmediatePredecessors(const std::vector<std::uint64_t>& els) const;
    [[nodiscard]] virtual POSet::DATASTORE GetImmediatePredecessors() const;
    
    virtual void FirstLE(LinearExtension& le) const;

    [[nodiscard]] virtual LatticeOfIdeals* GetLatticeOfIdeals();
    [[nodiscard]] virtual TreeOfIdeals* GetTreeOfIdeals();
    void ComputeTransitiveClosure(std::vector<std::set<std::uint64_t>>& adj);

    [[nodiscard]] virtual Tensor<std::uint8_t, 2> IncidenceMatrix() const;

    [[nodiscard]] virtual bool IsLessOrEqual(std::uint64_t a, std::uint64_t b) const noexcept {
        return upset_bits_[a].test(b);
    }
    
    [[nodiscard]] virtual bool GreaterThan(std::uint64_t e1, std::uint64_t e2) const noexcept {
        return (e1 != e2) && upset_bits_[e2].test(e1);
    }
    
    [[nodiscard]] virtual bool IsComparable(std::uint64_t a, std::uint64_t b) const noexcept {
        return upset_bits_[a].test(b) || upset_bits_[b].test(a);
    }
    
    [[nodiscard]] virtual std::vector<std::pair<std::string, std::string>> OrderRelation() const;
    [[nodiscard]] virtual Tensor<std::uint8_t, 2> CoverMatrix() const;
    
    [[nodiscard]] virtual std::optional<std::uint64_t> Meet(const std::vector<std::uint64_t>& insieme) const;
    [[nodiscard]] virtual std::optional<std::uint64_t> Join(const std::vector<std::uint64_t>& insieme) const;
    
    [[nodiscard]] virtual const POSet::DATASTORE& UpSets() const;
    [[nodiscard]] virtual BitSet DownSet(const std::vector<std::uint64_t>& els) const;
    [[nodiscard]] virtual bool IsDownSet(const std::vector<std::uint64_t>& els) const;
    [[nodiscard]] virtual BitSet UpSet(const std::vector<std::uint64_t>& els) const;
    [[nodiscard]] virtual bool IsUpSet(const std::vector<std::uint64_t>& els) const;
    [[nodiscard]] virtual BitSet ComparabilitySetOf(std::uint64_t e) const;
    [[nodiscard]] virtual BitSet IncomparabilitySetOf(std::uint64_t e) const;
    [[nodiscard]] virtual BitSet Maximals() const;
    [[nodiscard]] virtual BitSet Minimals() const;
    [[nodiscard]] virtual bool IsMaximal(std::uint64_t e) const;
    [[nodiscard]] virtual bool IsMinimal(std::uint64_t e) const;
    
    [[nodiscard]] virtual std::vector<std::pair<std::uint64_t, std::uint64_t>> CoverRelation() const;
    [[nodiscard]] virtual std::vector<std::pair<std::uint64_t, std::uint64_t>> Incomparabilities() const;
    [[nodiscard]] virtual bool IsExtensionOf(const POSet& p) const;
    [[nodiscard]] virtual std::vector<std::pair<std::uint64_t, std::uint64_t>> Comparabilities() const;
    
        
    [[nodiscard]] virtual std::shared_ptr<POSet> Clone() const;
    [[nodiscard]] virtual bool IsTotalOrder() const;
    
    [[nodiscard]] virtual Tensor<double, 2> BLSDominanceAbsolute() const;
    [[nodiscard]] virtual Tensor<double, 2> BLSDominanceRelative() const;
    
    static std::unique_ptr<POSet> LinearSum(const POSet& p1, const POSet& p2);
    static std::unique_ptr<POSet> DisjointSum(const POSet& p1, const POSet& p2);
    [[nodiscard]] std::unique_ptr<POSet> Lifting(const std::string& new_element) const;
    [[nodiscard]] virtual std::unique_ptr<POSet> Dual() const;
    
    /**
     * @brief Core HPC engine for evaluating multiple metrics over linear extensions.
     * @param[in] fles Vector of unique pointers to the evaluation functions (observers).
     * @param[in,out] leg The generator engine (TreeOfIdeals, BubleyDyer, etc.).
     * @param[out] eval_results Vector of unique pointers to the output matrices.
     * @param[out] le_count Counter updated with the number of processed extensions.
     * @param[out] end_process Flag set to true if the generator reaches its logical end.
     * @param[in] display_message Raw pointer to the progress display interface (observer).
     * @param[in] strategy Update strategy defining how/when to commit matrix updates.
     * @tparam GeneratorType The concrete type of the generator (deduced at compile time).
     */
    template <typename GeneratorType>
    static void evaluation(
                           const std::vector<std::unique_ptr<FunctionLinearExtension>>& fles,
                           GeneratorType& leg,
                           std::vector<std::unique_ptr<Tensor<double, 2>>>& eval_results,
                           std::uint64_t& le_count,
                           bool& end_process,
                           DisplayMessage* display_message,
                           POSet::EvaluationUpdateStrategy strategy)
    {
        if (strategy == POSet::EvaluationUpdateStrategy::Average) {
            evaluation_impl<POSet::EvaluationUpdateStrategy::Average>(fles, leg, eval_results, le_count, end_process, display_message);
        } else {
            evaluation_impl<POSet::EvaluationUpdateStrategy::Sum>(fles, leg, eval_results, le_count, end_process, display_message);
        }
    }
    
    /**
     * @brief Valutazione parallela multi-catena guidata da una sorgente lazy.
     * * Le catene non esistono tutte insieme: i worker le richiedono una alla
     * volta alla sorgente (Next() sotto lock, una chiamata per catena) e le
     * elaborano con evaluation_impl (stesso motore del caso sequenziale,
     * devirtualizzato sul tipo concreto del generatore). Il numero di catene
     * può quindi non essere noto a priori: si avviano i worker richiesti e chi
     * trova la sorgente esaurita esce subito. Il thread chiamante non calcola:
     * fa polling del progresso e resta l'unico a invocare display_message
     * (requisito API R: output solo dal thread principale).
     * * MERGE INCREMENTALE CON COMMIT ORDINATO: appena una catena termina, i
     * suoi tensori vengono fusi negli accumulatori (eval_results stessi) e
     * distrutti — memoria di picco O(n_threads) catene, non O(K). Il commit
     * avviene rigorosamente in ordine di produzione (le catene completate
     * fuori ordine attendono in un buffer limitato da n_threads): l'ordine
     * delle somme floating-point è quindi deterministico e i risultati sono
     * BIT-IDENTICI tra run a parità di sorgente e seed, qualunque sia
     * n_threads o lo scheduling.
     * * Merge in eval_results:
     * - Average: media pesata, peso = CurrentNumberOfLe() della catena (esatta).
     * - Sum: somma dei tensori locali.
     * * Errori: la prima eccezione (di una catena o della sorgente) ferma la
     * produzione di nuove catene; le catene già in corso terminano, poi
     * l'eccezione viene rilanciata sul thread chiamante.
     * @param source Sorgente di EvaluationChainTask (vedi EvaluationChainSource).
     * @param eval_results Tensori di output fusi.
     * @param le_count Totale LE elaborate in questa esecuzione (somma catene).
     * @param end_process true al completamento di tutte le catene.
     * @param n_threads Numero di thread worker; 0 = automatico. In ogni caso
     *        non più di hardware_concurrency.
     * @tparam SourceType Tipo concreto della sorgente (dedotto a compile time).
     */
    template <EvaluationChainSource SourceType>
    static void evaluation_parallel(
                           SourceType& source,
                           std::vector<std::unique_ptr<Tensor<double, 2>>>& eval_results,
                           std::uint64_t& le_count,
                           bool& end_process,
                           DisplayMessage* display_message,
                           POSet::EvaluationUpdateStrategy strategy,
                           std::uint64_t n_threads = 0);

    [[nodiscard]] virtual std::unique_ptr<LinearExtensionGenerator> CreateLinearExtensionGenerator();
    
    [[nodiscard]] Tensor<double, 2> ComputeMRP();
    
    static std::unique_ptr<POSet> Intersection(const POSet& p1, const POSet& p2);
    
    template<typename... Args>
    static std::unique_ptr<POSet> Intersection(const POSet& p_first, const Args&... p_next) {
        auto p_rest = Intersection(p_next...);
        return Intersection(p_first, *p_rest);
    }
private:
    static void SumUpdate(Tensor<double, 2>& ris, const FunctionLinearExtension& fle) noexcept;
    
    template <LinearExtensionProvider FLE>
    static void AverageUpdate(Tensor<double, 2>& ris, const FLE& fle, std::uint64_t le_number);
    
    template <EvaluationUpdateStrategy Strategy, typename GeneratorType>
    static void evaluation_impl(
                                const std::vector<std::unique_ptr<FunctionLinearExtension>>& fles,
                                GeneratorType& leg,
                                std::vector<std::unique_ptr<Tensor<double, 2>>>& eval_results,
                                std::uint64_t& le_count,
                                bool& end_process,
                                DisplayMessage* displayMessage);
};

/**
 * @brief Performs an incremental moving average update on a result matrix.
 * * This method updates the matrix 'ris' using the data provided by a Linear Extension
 * Provider. It implements an optimized dual-path strategy:
 * - **Dense Path**: Used when the provider size matches the matrix size. Performs
 * a single-pass update.
 * - **Sparse Path**: Used for compact data. First applies a global decay to the
 * matrix (N-1)/N and then adds the new contributions (val/N).
 * * @tparam FLE A type satisfying the LinearExtensionProvider concept.
 * @param ris The matrix (Tensor<double, 2>) to be updated in-place.
 * @param fle The linear extension data source.
 * @param le_number The current iteration count (1-based), used as the N factor.
 * * @note This is a template function to allow the compiler to perform aggressive
 * inlining and SIMD vectorization of the internal loops.
 */
template <LinearExtensionProvider FLE>
void POSet::AverageUpdate(Tensor<double, 2>& ris, const FLE& fle, std::uint64_t le_number) {
    if (le_number == 0) return;
    
    const std::uint64_t entries = fle.ResSize();
    const std::uint64_t total_cells = ris.size();
    const double inv_n = 1.0 / static_cast<double>(le_number);
    
    if (entries == total_cells) {
        for (std::uint64_t k = 0; k < entries; ++k) {
            ris(fle.At0(k), fle.At1(k)) += (fle.At2(k) - ris(fle.At0(k), fle.At1(k))) * inv_n;
        }
    } else {
        const double decay = (static_cast<double>(le_number) - 1.0) * inv_n;
        double* POSET_RESTRICT data = ris.data();
        for (std::uint64_t i = 0; i < total_cells; ++i) data[i] *= decay;
        
        for (std::uint64_t k = 0; k < entries; ++k) {
            ris(fle.At0(k), fle.At1(k)) += (fle.At2(k) * inv_n);
        }
    }
}

/**
 * @brief Motore interno di valutazione (Implementazione Template).
 * il compilatore genera una versione specifica per ogni strategia (Sum, Average, ecc.),
 * eliminando il costo dei rami condizionali (if) all'interno del loop critico.
 * Il template su GeneratorType devirtualizza HasNext()/Next() quando il call site
 * passa il tipo concreto (es. LEGBubleyDyer, LEGTreeOfIdeals — tutti final),
 * permettendo l'inlining dell'intero passo di generazione nel loop.
 */
template <POSet::EvaluationUpdateStrategy Strategy, typename GeneratorType>
void POSet::evaluation_impl(
                     const std::vector<std::unique_ptr<FunctionLinearExtension>>& fles,
                     GeneratorType& leg,
                     std::vector<std::unique_ptr<Tensor<double, 2>>>& eval_results,
                     std::uint64_t& le_count,
                     bool& end_process,
                     DisplayMessage* displayMessage)
{
    const std::uint64_t n_fles = fles.size();
    
    std::vector<FunctionLinearExtension*> raw_fles(n_fles);
    std::vector<Tensor<double, 2>*>         raw_results(n_fles);
    
    for (std::uint64_t k = 0; k < n_fles; ++k) {
        raw_fles[k]    = fles[k].get();
        raw_results[k] = eval_results[k].get();
    }
    
    displayMessage->Start();
    le_count = 1;
    
    while (true) {
        const LinearExtension& leval = leg.Get();
        const std::uint64_t current_le = leg.CurrentNumberOfLe();
        for (std::uint64_t k = 0; k < n_fles; ++k) {
            FunctionLinearExtension* fle = raw_fles[k];
            Tensor<double, 2>* eval_result = raw_results[k];
            
            (*fle)(leval);
            
            if constexpr (Strategy == POSet::EvaluationUpdateStrategy::Average) {
                POSet::AverageUpdate(*eval_result, *fle, current_le);
            }
            else if constexpr (Strategy == POSet::EvaluationUpdateStrategy::Sum) {
                POSet::SumUpdate(*eval_result, *fle);
            }
        }
        
        displayMessage->Display();
        
        if (!leg.HasNext()) break;

        leg.Next();
        ++le_count;
    }

    displayMessage->Stop();
    end_process = true;
}

template <EvaluationChainSource SourceType>
void POSet::evaluation_parallel(
                       SourceType& source,
                       std::vector<std::unique_ptr<Tensor<double, 2>>>& eval_results,
                       std::uint64_t& le_count,
                       bool& end_process,
                       DisplayMessage* display_message,
                       POSet::EvaluationUpdateStrategy strategy,
                       std::uint64_t n_threads)
{
    using GeneratorType = typename SourceType::GeneratorType;

    const std::size_t n_results = eval_results.size();

    // Numero di worker: n_threads esplicito, 0 = automatico; in ogni caso non
    // piu' di hardware_concurrency (fallback a 1 se l'OS non la riporta). Il
    // numero di catene puo' non essere noto a priori: un worker che trova la
    // sorgente esaurita esce subito.
    const unsigned int hw_concurrency = std::thread::hardware_concurrency();
    const std::uint64_t hw_threads = (hw_concurrency == 0) ? 1ULL : static_cast<std::uint64_t>(hw_concurrency);
    const std::uint64_t actual_threads =
        (n_threads == 0) ? hw_threads : std::min<std::uint64_t>(n_threads, hw_threads);

    // Gli output fanno anche da accumulatori del merge incrementale.
    for (auto& result : eval_results) {
        std::fill_n(result->data(), result->size(), 0.0);
    }

    // Stato condiviso protetto da un unico mutex: la sorgente (stateful, non
    // thread-safe) e il merge ordinato. Il lock viene preso una volta per
    // catena — mai per LE — quindi il costo e' trascurabile.
    std::mutex mtx;
    std::uint64_t next_index  = 0;   ///< Prossima catena da produrre.
    std::uint64_t next_commit = 0;   ///< Prossima catena da fondere nel risultato.
    double total_weight = 0.0;
    std::exception_ptr first_error = nullptr;

    // Catene completate fuori ordine, in attesa del proprio turno di commit.
    // Il commit rigorosamente in ordine di produzione rende deterministico
    // l'ordine delle somme floating-point: risultati bit-identici tra run,
    // qualunque sia n_threads o lo scheduling. Ogni worker detiene al piu'
    // un task alla volta, quindi |pending| <= actual_threads.
    struct CompletedChain {
        std::vector<std::unique_ptr<Tensor<double, 2>>> results;
        double weight;
    };
    std::map<std::uint64_t, CompletedChain> pending;

    // Progresso: counter striping. Ogni worker incrementa SOLO il proprio slot,
    // isolato su una cache line dedicata (128 = linea Apple Silicon, multiplo
    // della 64 x86; literal e non hardware_destructive_interference_size perche'
    // quest'ultima non e' disponibile/affidabile su tutte le toolchain R).
    // Un unico atomico condiviso farebbe rimbalzare la sua cache line tra i
    // core a ogni LE; qui la contesa e' zero e il thread principale (l'unico
    // autorizzato all'output R) somma gli slot nel polling.
    struct alignas(128) PaddedTicks {
        std::atomic<std::uint64_t> value{0};
    };
    std::vector<PaddedTicks> tick_slots(actual_threads);
    const auto sum_ticks = [&tick_slots]() noexcept {
        std::uint64_t s = 0;
        for (const auto& slot : tick_slots) s += slot.value.load(std::memory_order_relaxed);
        return s;
    };
    std::atomic<std::size_t> workers_done{0};

    {
        std::vector<std::jthread> workers;
        workers.reserve(actual_threads);

        for (std::uint64_t t = 0; t < actual_threads; ++t) {
            workers.emplace_back([&, t]() {
                DisplayMessageAtomicTick ticker(tick_slots[t].value);
                while (true) {
                    // --- Produzione (sotto lock): una catena dalla sorgente ---
                    std::optional<EvaluationChainTask<GeneratorType>> task;
                    std::uint64_t k = 0;
                    {
                        std::lock_guard lock(mtx);
                        if (first_error != nullptr || !source.HasNext()) break;
                        try {
                            task.emplace(source.Next());
                        } catch (...) {
                            first_error = std::current_exception();
                            break;
                        }
                        k = next_index++;
                    }

                    // --- Calcolo (fuori dal lock): workspace interamente
                    // locale al worker, contatori sullo stack (niente false
                    // sharing). Il progresso live e' affidato a 'ticker'
                    // (una Display() per LE).
                    std::uint64_t local_le_count = 0;
                    bool local_end_process = false;
                    try {
                        if (strategy == POSet::EvaluationUpdateStrategy::Average) {
                            evaluation_impl<POSet::EvaluationUpdateStrategy::Average>(
                                task->fles, *task->leg, task->results,
                                local_le_count, local_end_process, &ticker);
                        } else {
                            evaluation_impl<POSet::EvaluationUpdateStrategy::Sum>(
                                task->fles, *task->leg, task->results,
                                local_le_count, local_end_process, &ticker);
                        }
                    } catch (...) {
                        std::lock_guard lock(mtx);
                        if (first_error == nullptr) {
                            first_error = std::current_exception();
                        }
                        break;
                    }

                    const double weight =
                        (strategy == POSet::EvaluationUpdateStrategy::Average)
                            ? static_cast<double>(task->leg->CurrentNumberOfLe())
                            : 1.0;

                    // --- Commit ordinato (sotto lock): si accoda il risultato
                    // e si fonde tutto cio' che e' diventato contiguo. ---
                    {
                        std::lock_guard lock(mtx);
                        pending.emplace(k, CompletedChain{std::move(task->results), weight});
                        for (auto it = pending.find(next_commit);
                             it != pending.end();
                             it = pending.find(next_commit)) {
                            total_weight += it->second.weight;
                            for (std::size_t r = 0; r < n_results; ++r) {
                                double* POSET_RESTRICT out = eval_results[r]->data();
                                const double* POSET_RESTRICT src = it->second.results[r]->data();
                                const double w = it->second.weight;
                                const std::uint64_t cells = eval_results[r]->size();
                                for (std::uint64_t i = 0; i < cells; ++i) {
                                    out[i] += w * src[i];
                                }
                            }
                            pending.erase(it);
                            ++next_commit;
                        }
                    }
                    // Generatore e cloni delle funzioni muoiono qui, fuori dal lock.
                    task.reset();
                }
                workers_done.fetch_add(1, std::memory_order_release);
            });
        }

        display_message->Display("Hardware threads detected: " + std::to_string(hw_threads));
        display_message->Display("Number of threads (Actual): " + std::to_string(actual_threads));
        display_message->Start();
        while (workers_done.load(std::memory_order_acquire) < actual_threads) {
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            le_count = sum_ticks();
            display_message->Display();
        }
    } // join implicito di tutti i jthread

    // Dopo il join tutte le scritture dei worker sono visibili:
    // la somma degli slot e' il totale esatto delle LE elaborate.
    le_count = sum_ticks();
    display_message->Display("Chains: " + std::to_string(next_commit));
    display_message->Stop();

    if (first_error != nullptr) {
        std::rethrow_exception(first_error);
    }
    if (next_commit == 0) {
        throw MyException("POSet::evaluation_parallel: the source produced no chain.");
    }

    // ------------------------- FINALIZZAZIONE -------------------------
    // Average: gli accumulatori contengono sum_k(w_k * x_k); la divisione per
    // il peso totale completa la media pesata (esatta: peso = LE per catena).
    // Sum: gia' accumulato con peso 1, nulla da fare.
    if (strategy == POSet::EvaluationUpdateStrategy::Average) {
        if (total_weight <= 0.0) {
            throw MyException("POSet::evaluation_parallel: no linear extension processed.");
        }
        const double inv_weight = 1.0 / total_weight;
        for (std::size_t r = 0; r < n_results; ++r) {
            double* POSET_RESTRICT out = eval_results[r]->data();
            const std::uint64_t cells = eval_results[r]->size();
            for (std::uint64_t i = 0; i < cells; ++i) {
                out[i] *= inv_weight;
            }
        }
    }

    end_process = true;
}


