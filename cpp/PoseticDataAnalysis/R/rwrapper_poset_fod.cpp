/**
 * @file rwrapper_poset_fod.cpp
 * @brief R/C++ interface for First Order Dominance analysis on Partially Ordered Sets.
 */
#include <array>
#include <cstdint>
#include <cstdio>
#include <format>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>

#include "dimensionality_reduction.h"
#include "display_message.h"
#include "function_linear_extension.h"
#include "function_linear_extension_average_height.h"
#include "function_linear_extension_dominance.h"
#include "function_linear_extension_mann_whitney_dominance.h"
#include "function_linear_extension_mann_whitney_inferential_dominance.h"
#include "function_linear_extension_mutual_ranking_probability.h"
#include "linear_extension_generator.h"
#include "linear_extension_generator_tree_of_ideals.h"
#include "linear_extension_generator_bubley_dyer.h"
#include "linear_extension_generator_binary_variable.h"
#include "linear_extension_generator_lexicographic_from_linear.h"
#include "linear_generator_wrapper.h"
#include "tensor.h"
#include "my_exception.h"
#include "poset.h"
#include "binary_variable_poset.h"
#include "poset_wrapper.h"
#include "r_display.h"
#include "r_function_linear_extension.h"
#include "r_separation.h"
#include "random.h"
#include "rwrapper.h"
#include "rwrapper_conversion.h"
#include "separation.h"
#include "first_order_dominance_analysis.h"
#include "evaluation_chain_source_bubley_dyer_from_le_matrix.h"

namespace {
/**
 * @brief Configurazione delle metriche da calcolare.
 */
struct MetricsConfig {
    bool use_dominance = false;
    bool use_mann_whitney = false;
    bool use_inferential = false;
    bool use_mrp = false;

    // 0 = Dominance, 1 = Mann-Whitney, 2 = Inferential, 3 = MRP
    std::array<std::uint64_t, 4> sequence = {0, 0, 0, 0};
    std::size_t count = 0;

    void SetFromBooleanMask(std::size_t index, bool is_active) {
        if (!is_active) return;

        sequence[count++] = index;

        if (index == 0) use_dominance = true;
        if (index == 1) use_mann_whitney = true;
        if (index == 2) use_inferential = true;
        if (index == 3) use_mrp = true;
    }
};



/**
 * @brief Splits a string into multiple tokens using a specified delimiter.
 * @details Uses a pre-allocated buffer to achieve zero dynamic memory allocations
 * during repeated calls in high-performance loops.
 */
inline void SplitString(std::string_view s, std::string_view delimiter, std::vector<std::string_view>& tokens) {
    tokens.clear(); // Costo O(1), non libera la memoria heap sottostante
    
    if (delimiter.empty()) {
        tokens.push_back(s);
        return;
    }
    
    std::size_t start = 0;
    std::size_t end = s.find(delimiter);
    
    while (end != std::string_view::npos) {
        tokens.push_back(s.substr(start, end - start));
        start = end + delimiter.length();
        end = s.find(delimiter, start);
    }
    tokens.push_back(s.substr(start));
}

/**
 * @brief Instantiates the requested FunctionLinearExtension metrics based on the provided configuration.
 * * @details This function processes the stack-allocated `MetricsConfig` to determine which
 * metrics (e.g., Dominance, Mann-Whitney) need to be evaluated. It initializes the corresponding
 * objects, maps their names, and pre-allocates the memory-contiguous result matrices
 * required for the High-Performance Computing (HPC) evaluation phase.
 * * @param config The configuration struct defining which metrics to compute and their exact sequence.
 * @param poset_size Total number of elements in the partially ordered set.
 * @param freq_matrix Adjacency frequency matrix, optimized in a row-major layout for cache locality.
 * @param labels Categorical labels corresponding to the rows and columns of the matrix.
 * @param obs_vals Observation values translated into their corresponding global POSet element IDs.
 * @param subpopulation_count Vector containing sub-population sizes, required for inferential metrics.
 * @param total_bins The total number of bins to use for distribution analysis (e.g., in Mann-Whitney).
 * @param poset Pointer to the (single) POSet, required by the MRP metric;
 * nullptr when multiple POSets are provided (MRP not supported in that case).
 * * @return A `std::tuple` containing three aligned vectors:
 * - `std::vector<std::unique_ptr<FunctionLinearExtension>>`: The instantiated metric objects.
 * - `std::vector<std::string>`: The string names of the instantiated metrics.
 * - `std::vector<std::unique_ptr<Tensor<double, 2>>>`: Pre-allocated matrices for evaluation results.
 * * @throws MyException If an invalid or unknown metric identifier is found within the configuration sequence.
 */
std::tuple<
std::vector<std::unique_ptr<FunctionLinearExtension>>,
std::vector<std::string>,
std::vector<std::unique_ptr<Tensor<double, 2>>>
>
InstantiateMetrics(
                   const MetricsConfig& config,
                   std::size_t poset_size,
                   const Tensor<double, 2>& freq_matrix,
                   const Tensor<std::string_view, 1>& labels,
                   const std::vector<std::uint32_t>& obs_vals,
                   const Tensor<double, 1>& subpopulation_count,
                   std::uint64_t total_bins,
                   const POSet* poset)
{
    std::vector<std::unique_ptr<FunctionLinearExtension>> fles;
    std::vector<std::string> functions_name;
    std::vector<std::unique_ptr<Tensor<double, 2>>> eval_results;
    
    fles.reserve(config.count);
    functions_name.reserve(config.count);
    eval_results.reserve(config.count);
    
    for (std::size_t i = 0; i < config.count; ++i) {
        std::uint64_t metric = config.sequence[i];
        std::unique_ptr<FunctionLinearExtension> fle;
        std::string name;
        
        switch (metric) {
            case 0:
                fle = std::make_unique<FLEDominance>(
                                                     poset_size, freq_matrix, labels, obs_vals);
                name = "Dominance";
                break;
            case 1:
                fle = std::make_unique<FLEMannWhitneyDominance>(
                                                                poset_size, freq_matrix, labels, obs_vals, total_bins);
                name = "MannWhitneyDominance";
                break;
            case 2:
                fle = std::make_unique<FLEMannWhitneyInferentialDominance>(
                                                                           poset_size, freq_matrix, labels, subpopulation_count, obs_vals, total_bins);
                name = "MannWhitneyInferentialDominance";
                break;
            case 3:
                if (poset == nullptr) [[unlikely]] {
                    throw MyException("The 'MRP' metric requires a single POSet.");
                }
                fle = std::make_unique<FLEMutualRankingProbability>(poset);
                name = "MRP";
                break;
            default:
                throw MyException("Invalid value in 'metrics_r'. Expected: 0, 1, 2, or 3.");
        }
        
        const auto shape = fle->Shape();
        eval_results.push_back(std::make_unique<Tensor<double, 2>>(std::array<std::uint64_t, 2>{shape.at(0), shape.at(1)}, 0.0));
        
        fles.push_back(std::move(fle));
        functions_name.push_back(std::move(name));
    }
    
    return {std::move(fles), std::move(functions_name), std::move(eval_results)};
}

/**
 * @typedef InstantiateGeneratorVariant
 * @brief A type-safe union containing all possible linear extension generator implementations.
 * * @details This variant enables static polymorphism, allowing the compiler to perform
 * devirtualization and inlining. Unlike a base class pointer, this resides on the stack,
 * improving cache locality and eliminating heap allocation overhead in HPC hot paths.
 */
using InstantiateGeneratorVariant = std::variant<
LEGLexicographicFromLinear,
LEGBinaryVariable,
LEGTreeOfIdeals,
LEGBubleyDyer
>;

/**
 * @struct InstantiateGeneratorResult
 * @brief Outcome of the generator instantiation process.
 * * Bundles the concrete generator instance with its metadata (name) to avoid
 * string allocations or runtime type identification (RTTI).
 */
struct InstantiateGeneratorResult {
    InstantiateGeneratorVariant instance;
    std::string_view name;
    /// Seme effettivamente usato; vuoto per i generatori deterministici
    /// (enumerazione esatta, lessicografico, variabili binarie).
    std::optional<std::uint64_t> seed;
};

/**
 * @brief Creates the appropriate Linear Extension Generator instance (not yet started).
 * * @details Selects a generator based on the input POSet properties and the
 * requested generation mode (exact vs. MCMC). Internal helper: use
 * InstantiateGenerator, which also starts the generator.
 * * @param poset_wraps A vector of pointers to POSet wrappers. If size > 1,
 * LEGLexicographicFromLinear is used.
 * @param count       Optional limit for linear extensions. If provided, triggers
 * MCMC sampling (LEGBubleyDyer).
 * @param seed        Optional seed for stochastic generators.
 * * @return A GeneratorResult containing the move-constructed generator and its name.
 * @throws MyException If 'poset_wraps' is empty, if multiple POSets are provided
 * but are not of type kLin, or if internal POSet casting fails.
 */
[[nodiscard]] InstantiateGeneratorResult MakeGeneratorInstance(
                                                               const std::vector<const POSetWrap*>& poset_wraps,
                                                               std::optional<std::uint64_t> count,
                                                               std::optional<std::uint64_t> seed)
{
    if (poset_wraps.empty()) [[unlikely]] {
        throw MyException("InstantiateGenerator: 'poset_wraps' cannot be empty.");
    }
    
    // --- CASO A: Più POSet (Composizione) ---
    if (poset_wraps.size() > 1) {
        std::vector<std::uint64_t> group_sizes;
        group_sizes.reserve(poset_wraps.size());
        
        for (const auto* pw : poset_wraps) {
            if (pw->GetType() != POSetWrap::PosetType::kLin) [[unlikely]] {
                throw MyException("All POSets must be of type kLin for multiple POSets.");
            }
            group_sizes.push_back(pw->GetPOSet()->size());
        }
        
        return {
            .instance = LEGLexicographicFromLinear(std::move(group_sizes)),
            .name = "Lexicographic"
        };
    }
    
    // --- CASO B: POSet singolo ---
    const POSetWrap* poset_wrap = poset_wraps.front();
    POSet* poset = poset_wrap->GetPOSet();
    
    // Caso senza conteggio limitato (Deterministico/Esatto)
    if (!count.has_value()) {
        if (poset_wrap->GetType() == POSetWrap::PosetType::kBinaryVariable) {
            // HPC Optimization: static_cast (zero overhead rispetto a dynamic_cast)
            auto* binary_poset = static_cast<BinaryVariablePOSet*>(poset);
            return {
                .instance = LEGBinaryVariable(binary_poset->NumberOfVariables()),
                .name = "Binary"
            };
        } else {
            auto lattice = poset->GetLatticeOfIdeals();
            return {
                .instance = LEGTreeOfIdeals(poset->size(), *lattice),
                .name = "TreeOfIdeals"
            };
        }
    }
    // Caso con conteggio (Approssimato/MCMC)
    else {
        if (!seed.has_value()) {
            seed = Random::GLOBAL.RndNextInt(0, std::numeric_limits<std::uint64_t>::max());
        }
        
        // Creazione esplicita e move nel variant
        auto rng = std::make_unique<Random>(*seed);
        return {
            .instance = LEGBubleyDyer(poset, std::move(rng)),
            .name = "BubleyDyer",
            .seed = seed
        };
    }
}

/**
 * @brief Instantiates AND starts the appropriate Linear Extension Generator.
 * * @details After the instance is assigned into the variant, the generator's
 * Start() method is invoked with the value of 'count' (0 when 'count' is
 * not provided, i.e. unlimited/exact generation). std::visit dispatches the
 * call to the concrete generator type at compile time.
 * * @param poset_wraps A vector of pointers to POSet wrappers.
 * @param count       Optional limit for linear extensions; forwarded to Start()
 * (0 if absent). If provided, triggers MCMC sampling (LEGBubleyDyer).
 * @param seed        Optional seed for stochastic generators.
 * * @return A GeneratorResult containing the started generator and its name.
 * @throws MyException See MakeGeneratorInstance.
 */
[[nodiscard]] InstantiateGeneratorResult InstantiateGenerator(
                                                              const std::vector<const POSetWrap*>& poset_wraps,
                                                              std::optional<std::uint64_t> count,
                                                              std::optional<std::uint64_t> seed)
{
    InstantiateGeneratorResult result = MakeGeneratorInstance(poset_wraps, count, seed);
    
    // Avvio del generatore subito dopo l'assegnazione di .instance:
    // se 'count' non ha valore viene passato 0.
    std::visit([&](auto& leg) { leg.Start(count.value_or(0)); }, result.instance);
    
    return result;
}


/**
 * @brief Performs the linear extension evaluation.
 * @details The generator must have already been started (see InstantiateGenerator,
 * which invokes Start() with the requested count). Sets up the display message
 * interval and delegates the core computation to the `POSet::evaluation` method.
 * @param output_interval Optional interval (in seconds) for R console updates.
 * @param le_generator Reference to the active Linear Extension Generator.
 * @param eval_results Output vector where evaluated matrices will be stored.
 * @param fles Vector of functional extensions (metrics) to evaluate.
 * @param le_count Reference to an integer storing the total linear extensions processed.
 * @tparam GeneratorType Concrete type of the generator (deduced at compile time).
 */

template <typename GeneratorType>
void BuildEvaluation(std::optional<std::uint64_t> output_interval,
                     GeneratorType& le_generator,
                     std::vector<std::unique_ptr<Tensor<double, 2>>> &eval_results,
                     std::vector<std::unique_ptr<FunctionLinearExtension>> &fles,
                     std::uint64_t &le_count)
{
    // NB: il generatore è già stato avviato con Start(count) in InstantiateGenerator.
    le_count = 0;
    bool end_process = false;
    
    std::unique_ptr<DisplayMessage> display;
    if (!output_interval.has_value()) {
        display = std::make_unique<DisplayMessageNull>();
    } else {
        display = std::make_unique<DisplayMessageEvaluationExactR>(le_count, output_interval.value());
    }
    
    POSet::evaluation(fles, le_generator, eval_results, le_count, end_process,
                      display.get(), POSet::EvaluationUpdateStrategy::Average);
}

} // namespace

// ===========================================================================
// R-callable functions
// ===========================================================================

extern "C" {
/**
 * @brief Computes First Order Dominance evaluation natively from R parameters.
 * @param poset_r The target POSet as an external pointer (or an R list of external pointers).
 * @param freq_matrix_r The frequency matrix (double or integer) specifying node associations.
 * @param subpopulation_count_r Optional numeric vector of sub-population sizes
 * (mandatory when the inferential metric is requested).
 * @param metrics_r R logical/integer/double vector of length 4 denoting the
 * requested metrics: (Dominance, MannWhitneyDominance,
 * MannWhitneyInferentialDominance, MRP). The MRP metric requires a single
 * POSet and returns only the element-by-element binMatrix (no FOD analysis).
 * @param total_bins_r Positive integer scalar defining the number of bins
 * (mandatory when a Mann-Whitney metric is requested, otherwise it may be NULL).
 * @param count_r Optional integer vector for counting constraints.
 * @param seed_r Optional random seed for generator consistency. When
 * 'linear_extensions_r' is provided it may also be a vector with one seed per chain.
 * @param output_interval_in_sec_r Refresh rate of logging to the R console.
 * @param sep_r Separation parameters.
 * @param linear_extensions_r Optional character matrix (one COLUMN per linear
 * extension, one row per poset element). NULL: comportamento invariato.
 * Altrimenti ogni colonna diventa il punto di partenza di una catena
 * Bubley-Dyer indipendente e la valutazione usa POSet::evaluation_parallel
 * ('count_r' estensioni per catena, obbligatorio). Richiede un POSet singolo.
 * @param n_threads_r Optional positive integer: numero di thread worker per la
 * valutazione multi-catena. NULL o 0 = automatico (hardware_concurrency);
 * in ogni caso non piu' di hardware_concurrency (i worker in eccesso rispetto
 * alle catene escono subito). Ignorato (con warning) se 'linear_extensions_r'
 * e' NULL. Non influenza i risultati numerici (bit-identici a parita' di
 * catene e seed), solo il tempo di esecuzione.
 * @return SEXP An R named list containing the result lists of each evaluated metric.
 */
SEXP FirstOrderDominanceAnalysis(
                                 SEXP poset_r, SEXP freq_matrix_r,
                                 SEXP subpopulation_count_r, SEXP metrics_r,
                                 SEXP total_bins_r, SEXP count_r,
                                 SEXP seed_r, SEXP output_interval_in_sec_r,
                                 SEXP sep_r, SEXP linear_extensions_r,
                                 SEXP n_threads_r) {
    
    // Rf_error esegue un longjmp che salta i distruttori C++: per questo il
    // messaggio d'errore viene copiato in un buffer sullo stack e Rf_error
    // viene invocato solo quando tutti gli oggetti C++ (incluso il
    // RProtectGuard) sono già stati distrutti regolarmente.
    constexpr std::size_t kErrBufSize = 1024;
    char error_buffer[kErrBufSize];
    bool has_error = false;
    SEXP result_r = R_NilValue;
    
    try {
        RProtectGuard guard;
        
        // 1. Estrazione del POSet o della lista di POSet
        std::vector<const POSetWrap*> poset_wraps;
        if (TYPEOF(poset_r) == VECSXP) {
            const int len = Rf_length(poset_r);
            poset_wraps.reserve(static_cast<std::size_t>(len));
            for (int i = 0; i < len; ++i) {
                poset_wraps.push_back(RConvert::ToPOSetWrap(VECTOR_ELT(poset_r, i)));
            }
        } else {
            poset_wraps.push_back(RConvert::ToPOSetWrap(poset_r));
        }
        
        if (poset_wraps.empty()) [[unlikely]] {
            throw MyException("The 'poset_r' parameter cannot be empty.");
        }
        
        // Validazione del parametro sep_r
        std::string_view sep;
        if (poset_wraps.size() > 1) {
            if (Rf_isNull(sep_r) || TYPEOF(sep_r) != STRSXP || Rf_length(sep_r) == 0) [[unlikely]] {
                throw MyException("The separator parameter 'sep_r' is mandatory when multiple POSets are provided.");
            }
            sep = CHAR(STRING_ELT(sep_r, 0));
        }
        
        // 2. Validazione ed estrazione della matrice delle frequenze (NxM).
        // NB: accettiamo sia REALSXP sia INTSXP (es. l'output di table() è integer):
        // chiamare REAL() su una matrice integer provocherebbe un errore R (longjmp).
        if (!Rf_isMatrix(freq_matrix_r) ||
            (TYPEOF(freq_matrix_r) != REALSXP && TYPEOF(freq_matrix_r) != INTSXP)) [[unlikely]] {
            throw MyException("'freq_matrix' must be a numeric (double or integer) matrix.");
        }
        const bool freq_is_real = (TYPEOF(freq_matrix_r) == REALSXP);
        
        const std::size_t freq_nrow = static_cast<std::size_t>(Rf_nrows(freq_matrix_r));
        const std::size_t freq_ncol = static_cast<std::size_t>(Rf_ncols(freq_matrix_r));
        if (freq_nrow == 0 || freq_ncol == 0) [[unlikely]] {
            throw MyException("'freq_matrix' cannot have zero rows or columns.");
        }
        
        SEXP dimnames = Rf_getAttrib(freq_matrix_r, R_DimNamesSymbol);
        SEXP row_names = R_NilValue;
        SEXP col_names = R_NilValue;
        if (dimnames != R_NilValue && Rf_length(dimnames) == 2) {
            row_names = VECTOR_ELT(dimnames, 0);
            col_names = VECTOR_ELT(dimnames, 1);
        }
        
        // 3. Estrazione delle metriche
        const int n_metrics = Rf_length(metrics_r);
        if (n_metrics != 4) [[unlikely]] {
            throw MyException("The 'metrics_r' parameter must be a boolean/integer vector of exactly length 4.");
        }
        
        MetricsConfig config;
        
        // Lettura type-safe: da R può arrivare logical, integer o double (es. c(1, 0, 1)).
        // NA viene rifiutato esplicitamente: NA_LOGICAL/NA_INTEGER valgono INT_MIN,
        // quindi un semplice "!= 0" li interpreterebbe come TRUE.
        for (std::size_t i = 0; i < 4; ++i) {
            bool is_active = false;
            switch (TYPEOF(metrics_r)) {
                case LGLSXP: {
                    const int v = LOGICAL(metrics_r)[i];
                    if (v == NA_LOGICAL) [[unlikely]] {
                        throw MyException("NA is not allowed in 'metrics_r'.");
                    }
                    is_active = (v != 0);
                    break;
                }
                case INTSXP: {
                    const int v = INTEGER(metrics_r)[i];
                    if (v == NA_INTEGER) [[unlikely]] {
                        throw MyException("NA is not allowed in 'metrics_r'.");
                    }
                    is_active = (v != 0);
                    break;
                }
                case REALSXP: {
                    const double v = REAL(metrics_r)[i];
                    if (ISNAN(v)) [[unlikely]] {
                        throw MyException("NA is not allowed in 'metrics_r'.");
                    }
                    is_active = (v != 0.0);
                    break;
                }
                default:
                    throw MyException("'metrics_r' must be a logical or numeric vector of length 4.");
            }
            config.SetFromBooleanMask(i, is_active);
        }
        
        if (config.count == 0) [[unlikely]] {
            throw MyException("At least one metric must be selected (TRUE) in 'metrics_r'.");
        }

        // La MRP e' una matrice elemento x elemento del poset: richiede un
        // oggetto POSet materializzato, che nel caso multi-poset non esiste.
        if (config.use_mrp && poset_wraps.size() > 1) [[unlikely]] {
            throw MyException("The 'MRP' metric is not supported with multiple POSets.");
        }
        
        // Gestione e validazione di subpopulation_count
        bool requires_subpopulation = config.use_inferential;
        
        if (requires_subpopulation && Rf_isNull(subpopulation_count_r)) [[unlikely]] {
            throw MyException("subpopulation_count is mandatory when metric MannWhitneyInferentialDominance is selected.");
        }
        if (!requires_subpopulation && !Rf_isNull(subpopulation_count_r)) [[unlikely]] {
            forward_warning_to_r("'subpopulation_count' is only permitted when the 'MannWhitneyInferentialDominance' metric is requested.");
        }
        const Tensor<double, 1> subpopulation_count = RConvert::ToDouble1DTensor(subpopulation_count_r);
        if (requires_subpopulation && subpopulation_count.size() != freq_ncol) [[unlikely]] {
            throw MyException("subpopulation_count size must match the number of columns in freq_matrix.");
        }
        
        // 4. Estrazione di total_bins: obbligatorio (e > 0) solo quando è
        // richiesta una metrica Mann-Whitney; altrimenti può essere NULL.
        const bool needs_bins = config.use_mann_whitney || config.use_inferential;
        const std::optional<std::uint64_t> total_bins_opt = RConvert::ToOptionalUInt(total_bins_r);
        
        if (needs_bins && (!total_bins_opt.has_value() || *total_bins_opt == 0)) [[unlikely]] {
            throw MyException("'total_bins' must be a positive integer when a Mann-Whitney metric is selected.");
        }
        const std::uint64_t total_bins = total_bins_opt.value_or(0);
        
        /// 5. Estrazione di count, seed, n_threads e output_ogni (opzionali)
        const bool has_le_matrix = !Rf_isNull(linear_extensions_r);
        auto count = RConvert::ToOptionalUInt(count_r);
        std::optional<std::uint64_t> output_interval = RConvert::ToOptionalUInt(output_interval_in_sec_r);

        // n_threads: ha senso solo per la valutazione multi-catena.
        const std::optional<std::uint64_t> n_threads_opt = RConvert::ToOptionalUInt(n_threads_r);
        if (!has_le_matrix && n_threads_opt.has_value()) {
            forward_warning_to_r("'n_threads' is ignored unless 'linear_extensions' is provided.");
        }

        // Nel ramo multi-catena 'seed_r' può essere un vettore (un seed per
        // catena): viene interpretato più avanti da RConvert::ToChainSeeds.
        std::optional<std::uint64_t> seed;
        if (!has_le_matrix) {
            seed = RConvert::ToOptionalSeed(seed_r);
            if (seed.has_value() && !count.has_value()) [[unlikely]] {
                throw MyException("The 'seed_r' parameter can only be provided if 'count_r' is present.");
            }
        }

        // ==============================================================================
        // COSTRUZIONE DEI PARAMETRI
        // ==============================================================================
        if (!Rf_isString(row_names) || !Rf_isString(col_names) ||
            static_cast<std::size_t>(Rf_length(row_names)) != freq_nrow ||
            static_cast<std::size_t>(Rf_length(col_names)) != freq_ncol) [[unlikely]] {
            throw MyException("The frequency matrix must have character labels for rows (observed values) and columns (groups).");
        }
        
        // A. Estrazione labels delle colonne (gruppi)
        Tensor<std::string_view, 1> labels({freq_ncol}, kUninitialized);
        for (std::size_t i = 0; i < freq_ncol; ++i) {
            SEXP char_sexp = STRING_ELT(col_names, static_cast<int>(i));
            labels(i) = std::string_view(CHAR(char_sexp), static_cast<std::size_t>(Rf_xlength(char_sexp)));
        }
        
        
        // B. Mappatura dei valori osservati
        std::vector<std::uint32_t> obs_vals(freq_nrow);
        std::size_t total_poset_size = 1;
        
        if (poset_wraps.size() > 1) {
            // Gli ID globali vengono memorizzati in uint32_t: il prodotto delle
            // dimensioni dei POSet deve quindi stare in 32 bit. Il controllo evita
            // anche l'overflow del prodotto stesso.
            constexpr std::size_t kMaxPosetSize = std::numeric_limits<std::uint32_t>::max();
            for (const auto* pw : poset_wraps) {
                const std::size_t sz = static_cast<std::size_t>(pw->GetPOSet()->size());
                if (sz == 0) [[unlikely]] {
                    throw MyException("POSets with zero elements are not allowed.");
                }
                if (total_poset_size > kMaxPosetSize / sz) [[unlikely]] {
                    throw MyException("The product of the POSet sizes exceeds the 32-bit element id space.");
                }
                total_poset_size *= sz;
            }
            
            std::vector<std::string_view> token_buffer;
            token_buffer.reserve(poset_wraps.size());
            
            for (std::size_t i = 0; i < freq_nrow; ++i) {
                std::string_view elem_name = CHAR(STRING_ELT(row_names, static_cast<int>(i)));
                
                SplitString(elem_name, sep, token_buffer);
                
                if (token_buffer.size() != poset_wraps.size()) [[unlikely]] {
                    throw MyException(std::format("The element label '{}' does not match the expected number of POSets ({}).",
                                                  elem_name, poset_wraps.size()));
                }
                
                // Aritmetica a 64 bit: il controllo sul prodotto delle dimensioni
                // (sopra) garantisce che global_id stia sempre in 32 bit.
                std::uint64_t global_id = 0;
                std::uint64_t multiplier = 1;
                
                for (int j = static_cast<int>(poset_wraps.size()) - 1; j >= 0; --j) {
                    const std::uint64_t id_j = poset_wraps[j]->GetPOSet()->GetElementId(token_buffer[j]);
                    global_id += id_j * multiplier;
                    multiplier *= poset_wraps[j]->GetPOSet()->size();
                }
                obs_vals[i] = static_cast<std::uint32_t>(global_id);
            }
        } else {
            const auto* pw = poset_wraps.front();
            total_poset_size = pw->GetPOSet()->size();
            
            for (std::size_t i = 0; i < freq_nrow; ++i) {
                // Costruito senza std::string_view intermedio per evitare l'allocazione se possibile
                obs_vals[i] = static_cast<std::uint32_t>(pw->GetPOSet()->GetElementId(CHAR(STRING_ELT(row_names, static_cast<int>(i)))));
            }
        }
        
        // C. Conversione matrice: R (col-major) -> C++ (row-major)
        
        Tensor<double, 2> freq_matrix({freq_nrow, freq_ncol}, kUninitialized);
        if (freq_is_real) {
            const double* freq_matrix_ptr = REAL(freq_matrix_r);
            for (std::size_t r = 0; r < freq_nrow; ++r) {
                for (std::size_t c = 0; c < freq_ncol; ++c) {
                    freq_matrix(r, c) = freq_matrix_ptr[r + c * freq_nrow];
                }
            }
        } else {
            const int* freq_matrix_ptr = INTEGER(freq_matrix_r);
            for (std::size_t r = 0; r < freq_nrow; ++r) {
                for (std::size_t c = 0; c < freq_ncol; ++c) {
                    const int v = freq_matrix_ptr[r + c * freq_nrow];
                    freq_matrix(r, c) = (v == NA_INTEGER) ? NA_REAL : static_cast<double>(v);
                }
            }
        }
        
        // ==============================================================================
        // ISTANZIAZIONE DELLE METRICHE, DEI NOMI E DELLE MATRICI DI RISULTATO
        // ==============================================================================
        
        auto instantiate_metrics = InstantiateMetrics(
                                                      config,
                                                      total_poset_size,
                                                      freq_matrix,
                                                      labels,
                                                      obs_vals,
                                                      subpopulation_count,
                                                      total_bins,
                                                      poset_wraps.size() == 1 ? poset_wraps.front()->GetPOSet() : nullptr
                                                      );
        auto& fles = std::get<0>(instantiate_metrics);
        auto& functions_name = std::get<1>(instantiate_metrics);
        auto& eval_results = std::get<2>(instantiate_metrics);
        // ==============================================================================
        // EVALUATION
        // ==============================================================================
        std::uint64_t le_count = 0;
        std::string_view le_type_name;
        // Seme (o semi, uno per catena) effettivamente usato: restituito a R
        // come stringa, resta R_NilValue per i generatori deterministici.
        SEXP seed_used_r = R_NilValue;

        if (!has_le_matrix) {
            // --- Percorso invariato: generatore singolo (esatto o MCMC) ---
            auto instantiate_generator = InstantiateGenerator(poset_wraps, count, seed);
            le_type_name = instantiate_generator.name;
            if (instantiate_generator.seed.has_value()) {
                seed_used_r = RCreate::FromSeed(guard, *instantiate_generator.seed);
            }

            std::visit([&](auto& concrete_leg) {
                BuildEvaluation(output_interval,
                                concrete_leg,
                                eval_results,
                                fles,
                                le_count);
            }, instantiate_generator.instance);
        } else {
            // --- Percorso multi-catena: una catena Bubley-Dyer per colonna ---
            if (poset_wraps.size() > 1) [[unlikely]] {
                throw MyException("'linear_extensions' is not supported with multiple POSets.");
            }
            if (!count.has_value() || *count == 0) [[unlikely]] {
                throw MyException("'count' (number of linear extensions PER CHAIN) is mandatory "
                                  "and must be > 0 when 'linear_extensions' is provided.");
            }
            POSet* poset = poset_wraps.front()->GetPOSet();

            // L'estrazione dai SEXP avviene QUI, sul main thread (l'API di R
            // non e' thread-safe): la sorgente riceve solo dati C++ e puo'
            // essere consumata dai worker.
            auto initial_les = RConvert::ToLinearExtensions(linear_extensions_r, *poset);
            const std::uint64_t n_chains = initial_les.size();
            auto seeds = RConvert::ToChainSeeds(seed_r, n_chains);
            // I semi per-catena vengono riportati a R prima di essere spostati
            // nella sorgente: reimmetterli riproduce esattamente le catene.
            seed_used_r = RCreate::FromSeeds(guard, seeds);

            // Sorgente lazy: ogni Next() costruisce la sola catena successiva
            // (LEGBubleyDyer con validazione della LE iniziale, cloni delle
            // funzioni, tensori locali). Con il merge incrementale di
            // evaluation_parallel la memoria di picco e' O(n_threads) catene.
            ECSBubleyDyerFromLEMatrix source(poset, std::move(initial_les), std::move(seeds),
                                 *count, fles, eval_results);

            std::unique_ptr<DisplayMessage> display;
            if (!output_interval.has_value()) {
                display = std::make_unique<DisplayMessageNull>();
            } else {
                display = std::make_unique<DisplayMessageEvaluationExactR>(le_count, output_interval.value());
            }

            bool end_process = false;
            POSet::evaluation_parallel(source, eval_results, le_count, end_process,
                                       display.get(), POSet::EvaluationUpdateStrategy::Average,
                                       n_threads_opt.value_or(0));
            le_type_name = "BubleyDyerMultiChain";
        }
        
        // ==============================================================================
        // FIRST ORDER DOMINANCE ANALYSIS
        // ==============================================================================
        
        // std::optional: la MRP non ha analisi FOD (la sua matrice e'
        // elemento x elemento del poset, non gruppo x gruppo) e per lei
        // entrambi gli slot restano vuoti (-> R_NilValue nel risultato).
        std::vector<std::optional<Tensor<double, 2>>> fod_analyses;
        std::vector<std::optional<Tensor<double, 2>>> analyze_matricies;
        // 1. Riserviamo lo spazio per evitare i ridimensionamenti dinamici dei vettori
        fod_analyses.reserve(eval_results.size());
        analyze_matricies.reserve(eval_results.size());

        const std::uint64_t n_elements = static_cast<std::uint64_t>(freq_ncol * freq_ncol);

        for (std::size_t k = 0; k < eval_results.size(); ++k) {
            const auto& eval_matrix_ptr = eval_results[k];
            const std::uint64_t current_metric = config.sequence[k];

            if (current_metric == 3) {
                // MRP: solo binMatrix, nessuna analisi FOD.
                analyze_matricies.emplace_back(std::nullopt);
                fod_analyses.emplace_back(std::nullopt);
            }
            else if (current_metric == 0) {
                analyze_matricies.emplace_back(*eval_matrix_ptr);
                fod_analyses.emplace_back(MinMaxTransitiveClosure(*eval_matrix_ptr));
            }
            else {
                // 2. In-Place Construction: allochiamo la matrice DIRETTAMENTE nel vettore finale.
                // Questo azzera l'overhead di allocazioni temporanee sul thread.
                analyze_matricies.emplace_back(std::in_place, std::array<std::uint64_t, 2>{freq_ncol, freq_ncol}, kUninitialized);
                auto& analyze_matrix = *analyze_matricies.back();
                
                double* const __restrict dest = analyze_matrix.data();
                const auto& src_mat = *eval_matrix_ptr;
                
                // Check difensivo: la colonna di indice 'total_bins' e le prime
                // 'n_elements' righe devono esistere nella matrice di valutazione
                // (operator() del Tensor non fa bound checking in release).
                if (src_mat.Extent(0) < n_elements || src_mat.Extent(1) <= total_bins) [[unlikely]] {
                    throw MyException(std::format(
                                                  "Internal inconsistency: evaluation matrix is {}x{}, but {} rows and column index {} are required.",
                                                  src_mat.Extent(0), src_mat.Extent(1), n_elements, total_bins));
                }
                
#if defined(__clang__)
                // Soluzione nativa per Clang (macOS / Xcode)
#pragma clang loop vectorize(enable)
#elif defined(__GNUC__)
                // Soluzione nativa per GCC (Linux / Rtools su Windows)
#pragma GCC ivdep
#elif defined(_MSC_VER)
                // Soluzione nativa per il compilatore Microsoft
#pragma loop(ivdep)
#endif
                for (std::uint64_t idx = 0; idx < n_elements; ++idx) {
                    dest[idx] = src_mat(idx, total_bins);
                }
                
                fod_analyses.emplace_back(MinMaxTransitiveClosure(analyze_matrix));
            }
        }
        
        
        // ==============================================================================
        // PREPARAZIONE RISULTATO FINALE (LISTA DI ANALISI)
        // ==============================================================================
        
        const std::size_t n_results = fod_analyses.size();
        
        // Creiamo la lista principale e il vettore dei nomi delle metriche
        auto [result_list_r, metric_names_r] = RCreate::NamedList(guard, static_cast<int>(n_results + 1));
        
        SEXP labels_r = guard.Protect(Rf_allocVector(STRSXP, labels.size()));
        for (std::size_t j = 0; j < labels.size(); ++j) {
            SET_STRING_ELT(labels_r, j, Rf_mkCharLenCE(labels(j).data(), static_cast<int>(labels(j).size()), CE_UTF8));
        }
        
        for (std::size_t i = 0; i < n_results; ++i) {

            auto [list_r, names_r] = RCreate::NamedList(guard, 3);

            // Per la MRP gli optional sono vuoti: fodClosed/fodMatrix = NULL.
            SEXP fod_closed_r = R_NilValue;
            SEXP fod_matrix_r = R_NilValue;
            if (fod_analyses[i].has_value()) {
                fod_closed_r = RCreate::FromDoubleMatrix(guard, *fod_analyses[i]);
                RCreate::AttachDimNames(guard, fod_closed_r, labels_r, labels_r);

                fod_matrix_r = RCreate::FromDoubleMatrix(guard, *analyze_matricies[i]);
                RCreate::AttachDimNames(guard, fod_matrix_r, labels_r, labels_r);
            }

            SEXP bin_r = RCreate::FromDoubleMatrix(guard, *eval_results[i]);
            
            SEXP bin_rows_r = guard.Protect(Rf_allocVector(STRSXP, eval_results[i]->Extent(0)));
            for (std::size_t j = 0; j < eval_results[i]->Extent(0); ++j) {
                std::string_view row_name = fles[i]->GetRowNameAt(j);
                SET_STRING_ELT(
                               bin_rows_r,
                               j,
                               Rf_mkCharLenCE(row_name.data(), static_cast<int>(row_name.size()), CE_UTF8)
                               );
            }
            
            SEXP bin_cols_r = guard.Protect(Rf_allocVector(STRSXP, eval_results[i]->Extent(1)));
            for (std::size_t j = 0; j < eval_results[i]->Extent(1); ++j) {
                std::string_view col_name = fles[i]->GetColNameAt(j);
                SET_STRING_ELT(
                               bin_cols_r,
                               j,
                               Rf_mkCharLenCE(col_name.data(), static_cast<int>(col_name.size()), CE_UTF8)
                               );
            }
            
            RCreate::AttachDimNames(guard, bin_r, bin_rows_r, bin_cols_r);
            
            RCreate::SetListElement(list_r, names_r, 0, fod_closed_r, "fodClosed");
            RCreate::SetListElement(list_r, names_r, 1, fod_matrix_r, "fodMatrix");
            RCreate::SetListElement(list_r, names_r, 2, bin_r, "binMatrix");
            
            Rf_setAttrib(list_r, R_NamesSymbol, names_r);
            
            RCreate::SetListElement(result_list_r, metric_names_r, static_cast<int>(i),
                                    list_r, functions_name[i].c_str());
        }
        
        // Un std::string_view non è garantito null-terminated: usiamo
        // Rf_mkCharLenCE (che riceve la lunghezza esplicita) invece di Rf_mkString.
        SEXP le_type_r = guard.Protect(Rf_allocVector(STRSXP, 1));
        SET_STRING_ELT(le_type_r, 0,
                       Rf_mkCharLenCE(le_type_name.data(),
                                      static_cast<int>(le_type_name.size()),
                                      CE_UTF8));
        RCreate::SetListElement(
                                result_list_r,
                                metric_names_r,
                                static_cast<int>(n_results),
                                le_type_r,
                                "LEType"
                                );
        
        // Assegniamo i nomi alla lista principale
        Rf_setAttrib(result_list_r, R_NamesSymbol, metric_names_r);
        
        // Aggiungiamo anche le informazioni sul conteggio delle estensioni lineari come attributo
        Rf_setAttrib(result_list_r, Rf_install("le_count"), RCreate::FromDouble(guard, static_cast<double>(le_count)));

        // Seme (o vettore di semi, uno per catena) usato dal campionamento MCMC,
        // come stringa di cifre: reimmetterlo in 'seed' riproduce il risultato.
        // Assente con i generatori deterministici (enumerazione esatta).
        if (seed_used_r != R_NilValue) {
            Rf_setAttrib(result_list_r, Rf_install("seed"), seed_used_r);
        }
        
        // Assegna il risultato da restituire ad R
        result_r = result_list_r;
        
    } catch (const std::exception &ex) {
        std::snprintf(error_buffer, sizeof(error_buffer), "%s", ex.what());
        has_error = true;
    } catch (...) {
        std::snprintf(error_buffer, sizeof(error_buffer), "c++ exception (unknown reason)");
        has_error = true;
    }
    
    if (has_error) {
        // A questo punto tutti gli oggetti C++ sono stati distrutti (il guard ha
        // già eseguito UNPROTECT): il longjmp di Rf_error è quindi sicuro.
        Rf_error("%s", error_buffer);
    }
    
    return result_r;
}
}
