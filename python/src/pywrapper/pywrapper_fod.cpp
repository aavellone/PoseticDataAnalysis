/**
 * @file pywrapper_fod.cpp
 * @brief Entry-point Python per la First Order Dominance analysis.
 *
 * @details Speculare a rwrapper_poset_fod.cpp. Gli helper interni (MetricsConfig,
 * InstantiateMetrics, MakeGeneratorInstance, SplitString) sono ripresi dal
 * modulo R (sono R-free); l'I/O passa da Python e il reporter di progresso usa
 * py_display. Il path multi-catena (linear_extensions) usa evaluation_parallel
 * con il GIL rilasciato.
 */

#include "pywrapper.h"
#include "py_common.h"
#include "py_convert.h"
#include "py_display.h"

#include "binary_variable_poset.h"
#include "evaluation_chain_source_bubley_dyer_from_le_matrix.h"
#include "first_order_dominance_analysis.h"
#include "function_linear_extension.h"
#include "function_linear_extension_dominance.h"
#include "function_linear_extension_mann_whitney_dominance.h"
#include "function_linear_extension_mann_whitney_inferential_dominance.h"
#include "function_linear_extension_mutual_ranking_probability.h"
#include "linear_extension.h"
#include "linear_extension_generator_binary_variable.h"
#include "linear_extension_generator_bubley_dyer.h"
#include "linear_extension_generator_lexicographic_from_linear.h"
#include "linear_extension_generator_tree_of_ideals.h"
#include "my_exception.h"
#include "poset.h"
#include "poset_wrapper.h"
#include "random.h"
#include "separation.h"
#include "tensor.h"

#include <array>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <variant>
#include <vector>

namespace {

struct MetricsConfig {
    std::array<std::uint64_t, 4> sequence = {0, 0, 0, 0};  // 0=Dom,1=MW,2=Inferential,3=MRP
    std::size_t count = 0;
    bool use_mann_whitney = false;
    bool use_inferential = false;
    bool use_mrp = false;

    void Set(std::size_t index, bool active) {
        if (!active) return;
        sequence[count++] = index;
        if (index == 1) use_mann_whitney = true;
        if (index == 2) use_inferential = true;
        if (index == 3) use_mrp = true;
    }
};

void SplitString(std::string_view s, std::string_view delim, std::vector<std::string_view>& out) {
    out.clear();
    if (delim.empty()) {
        out.push_back(s);
        return;
    }
    std::size_t start = 0, end = s.find(delim);
    while (end != std::string_view::npos) {
        out.push_back(s.substr(start, end - start));
        start = end + delim.length();
        end = s.find(delim, start);
    }
    out.push_back(s.substr(start));
}

std::tuple<std::vector<std::unique_ptr<FunctionLinearExtension>>, std::vector<std::string>,
           std::vector<std::unique_ptr<Tensor<double, 2>>>>
InstantiateMetrics(const MetricsConfig& config, std::size_t poset_size,
                   const Tensor<double, 2>& freq_matrix, const Tensor<std::string_view, 1>& labels,
                   const std::vector<std::uint32_t>& obs_vals,
                   const Tensor<double, 1>& subpopulation_count, std::uint64_t total_bins,
                   const POSet* poset) {
    std::vector<std::unique_ptr<FunctionLinearExtension>> fles;
    std::vector<std::string> names;
    std::vector<std::unique_ptr<Tensor<double, 2>>> results;
    for (std::size_t i = 0; i < config.count; ++i) {
        std::unique_ptr<FunctionLinearExtension> fle;
        std::string name;
        switch (config.sequence[i]) {
            case 0:
                fle = std::make_unique<FLEDominance>(poset_size, freq_matrix, labels, obs_vals);
                name = "Dominance";
                break;
            case 1:
                fle = std::make_unique<FLEMannWhitneyDominance>(poset_size, freq_matrix, labels,
                                                                obs_vals, total_bins);
                name = "MannWhitneyDominance";
                break;
            case 2:
                fle = std::make_unique<FLEMannWhitneyInferentialDominance>(
                    poset_size, freq_matrix, labels, subpopulation_count, obs_vals, total_bins);
                name = "MannWhitneyInferentialDominance";
                break;
            case 3:
                if (poset == nullptr) {
                    throw MyException("The 'MRP' metric requires a single POSet.");
                }
                fle = std::make_unique<FLEMutualRankingProbability>(poset);
                name = "MRP";
                break;
            default:
                throw MyException("Invalid metric id.");
        }
        const auto shape = fle->Shape();
        results.push_back(std::make_unique<Tensor<double, 2>>(
            std::array<std::uint64_t, 2>{shape.at(0), shape.at(1)}, 0.0));
        fles.push_back(std::move(fle));
        names.push_back(std::move(name));
    }
    return {std::move(fles), std::move(names), std::move(results)};
}

using GenVariant = std::variant<LEGLexicographicFromLinear, LEGBinaryVariable, LEGTreeOfIdeals,
                                LEGBubleyDyer>;
struct GenResult {
    GenVariant instance;
    std::string_view name;
    std::optional<std::uint64_t> seed;
};

GenResult MakeGenerator(const std::vector<const POSetWrap*>& wraps,
                        std::optional<std::uint64_t> count, std::optional<std::uint64_t> seed) {
    if (wraps.size() > 1) {
        std::vector<std::uint64_t> group_sizes;
        group_sizes.reserve(wraps.size());
        for (const auto* pw : wraps) {
            if (pw->GetType() != POSetWrap::PosetType::kLin) {
                throw MyException("All POSets must be linear (chains) for the multi-POSet case.");
            }
            group_sizes.push_back(pw->GetPOSet()->size());
        }
        return {LEGLexicographicFromLinear(std::move(group_sizes)), "Lexicographic", std::nullopt};
    }
    const POSetWrap* wrap = wraps.front();
    POSet* poset = wrap->GetPOSet();
    if (!count.has_value()) {
        if (wrap->GetType() == POSetWrap::PosetType::kBinaryVariable) {
            auto* bp = static_cast<BinaryVariablePOSet*>(poset);
            return {LEGBinaryVariable(bp->NumberOfVariables()), "Binary", std::nullopt};
        }
        auto lattice = poset->GetLatticeOfIdeals();
        return {LEGTreeOfIdeals(poset->size(), *lattice), "TreeOfIdeals", std::nullopt};
    }
    if (!seed.has_value()) {
        seed = Random::GLOBAL.RndNextInt(0, std::numeric_limits<std::uint64_t>::max());
    }
    return {LEGBubleyDyer(poset, std::make_unique<Random>(*seed)), "BubleyDyer", seed};
}

/// SplitMix64 mix identico a RConvert::ToChainSeeds (per compatibilita' bit-a-bit con R).
std::uint64_t SplitMix64Mix(std::uint64_t x) noexcept {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

std::vector<std::uint64_t> ToChainSeeds(PyObject* seed_obj, std::uint64_t n_chains) {
    // Vettore esplicito di seed (uno per catena)?
    if (seed_obj != nullptr && seed_obj != Py_None && PySequence_Check(seed_obj) &&
        !PyUnicode_Check(seed_obj)) {
        auto seeds_vec = PyConvert::ToSeedVector(seed_obj);
        if (seeds_vec.size() > 1) {
            if (static_cast<std::uint64_t>(seeds_vec.size()) != n_chains) {
                throw std::invalid_argument(
                    "'seed' must be None, a scalar, or one seed per chain.");
            }
            return seeds_vec;
        }
    }
    std::optional<std::uint64_t> base = PyConvert::ToOptionalSeed(seed_obj);
    const std::uint64_t b = base.value_or(
        Random::GLOBAL.RndNextInt(0, std::numeric_limits<std::uint64_t>::max()));
    std::vector<std::uint64_t> seeds;
    seeds.reserve(n_chains);
    for (std::uint64_t k = 0; k < n_chains; ++k) {
        seeds.push_back(SplitMix64Mix(b + k));
    }
    return seeds;
}

/// Converte una list di estensioni (ognuna list di nomi) in vector<LinearExtension>.
std::vector<LinearExtension> ToLinearExtensions(PyObject* les_obj, const POSet& poset) {
    PyRef fast(PySequence_Fast(les_obj, "linear_extensions must be a sequence"));
    if (!fast) {
        throw std::invalid_argument("linear_extensions must be a sequence of extensions.");
    }
    const Py_ssize_t k = PySequence_Fast_GET_SIZE(fast.get());
    const std::uint64_t n = poset.size();
    std::vector<LinearExtension> out;
    out.reserve(static_cast<std::size_t>(k));
    for (Py_ssize_t c = 0; c < k; ++c) {
        auto names = PyConvert::ToStringVector(PySequence_Fast_GET_ITEM(fast.get(), c),
                                               "linear_extensions");
        if (static_cast<std::uint64_t>(names.size()) != n) {
            throw std::invalid_argument(
                "each linear extension must list exactly one name per poset element.");
        }
        LinearExtension le(n);
        for (std::size_t r = 0; r < names.size(); ++r) {
            le.Set(r, poset.GetElementId(names[r]));
        }
        out.push_back(std::move(le));
    }
    return out;
}

/// Costruisce il dict finale per una metrica: {fodClosed, fodMatrix, binMatrix, rows, cols}.
PyObject* MetricToPy(const std::optional<Tensor<double, 2>>& fod_closed,
                     const std::optional<Tensor<double, 2>>& fod_matrix,
                     const Tensor<double, 2>& bin, FunctionLinearExtension* fle) {
    PyRef dict(PyDict_New());
    if (!dict) {
        throw std::bad_alloc();
    }
    PyCreate::DictSetSteal(dict.get(), "fodClosed",
                           fod_closed.has_value() ? PyCreate::FromDoubleMatrix(*fod_closed)
                                                  : (Py_INCREF(Py_None), Py_None));
    PyCreate::DictSetSteal(dict.get(), "fodMatrix",
                           fod_matrix.has_value() ? PyCreate::FromDoubleMatrix(*fod_matrix)
                                                  : (Py_INCREF(Py_None), Py_None));
    PyCreate::DictSetSteal(dict.get(), "binMatrix", PyCreate::FromDoubleMatrix(bin));

    const std::size_t nrow = static_cast<std::size_t>(bin.Extent(0));
    const std::size_t ncol = static_cast<std::size_t>(bin.Extent(1));
    PyRef rows(PyList_New(static_cast<Py_ssize_t>(nrow)));
    PyRef cols(PyList_New(static_cast<Py_ssize_t>(ncol)));
    if (!rows || !cols) {
        throw std::bad_alloc();
    }
    for (std::size_t j = 0; j < nrow; ++j) {
        PyList_SET_ITEM(rows.get(), static_cast<Py_ssize_t>(j),
                        PyCreate::FromString(fle->GetRowNameAt(j)));
    }
    for (std::size_t j = 0; j < ncol; ++j) {
        PyList_SET_ITEM(cols.get(), static_cast<Py_ssize_t>(j),
                        PyCreate::FromString(fle->GetColNameAt(j)));
    }
    PyCreate::DictSetSteal(dict.get(), "rows", rows.release());
    PyCreate::DictSetSteal(dict.get(), "cols", cols.release());
    return dict.release();
}

}  // namespace

// first_order_dominance_analysis(posets, freq_matrix, row_labels, col_labels, metrics,
//   subpopulation_count, total_bins, count, seed, output_every, sep, linear_extensions, n_threads)
PyObject* pyx_first_order_dominance_analysis(PyObject* /*self*/, PyObject* args) {
    PyObject *posets_obj, *freq_obj, *rowlab_obj, *collab_obj, *metrics_obj, *subpop_obj,
        *bins_obj, *count_obj, *seed_obj, *output_obj, *sep_obj, *les_obj, *nthreads_obj;
    if (!PyArg_ParseTuple(args, "OOOOOOOOOOOOO", &posets_obj, &freq_obj, &rowlab_obj, &collab_obj,
                          &metrics_obj, &subpop_obj, &bins_obj, &count_obj, &seed_obj, &output_obj,
                          &sep_obj, &les_obj, &nthreads_obj)) {
        return nullptr;
    }
    try {
        // 1. POSet(s)
        std::vector<const POSetWrap*> wraps;
        {
            PyRef fast(PySequence_Fast(posets_obj, "posets must be a sequence of POSet handles"));
            if (!fast) {
                throw std::invalid_argument("posets must be a sequence of POSet objects.");
            }
            const Py_ssize_t np = PySequence_Fast_GET_SIZE(fast.get());
            if (np == 0) {
                throw std::invalid_argument("at least one POSet is required.");
            }
            for (Py_ssize_t i = 0; i < np; ++i) {
                wraps.push_back(GetPOSetWrap(PySequence_Fast_GET_ITEM(fast.get(), i)));
            }
        }

        std::string sep;
        if (wraps.size() > 1) {
            if (sep_obj == Py_None) {
                throw std::invalid_argument("'sep' is required when multiple POSets are provided.");
            }
            sep = PyConvert::ToString(sep_obj, "sep");
        }

        // 2. freq matrix + labels (owned so string_views stay valid)
        auto freq_matrix = PyConvert::ToDoubleMatrix(freq_obj, "freq_matrix");
        const std::size_t freq_nrow = static_cast<std::size_t>(freq_matrix.Extent(0));
        const std::size_t freq_ncol = static_cast<std::size_t>(freq_matrix.Extent(1));
        auto row_names = PyConvert::ToStringVector(rowlab_obj, "row_labels");
        auto col_names = PyConvert::ToStringVector(collab_obj, "col_labels");
        if (row_names.size() != freq_nrow || col_names.size() != freq_ncol) {
            throw std::invalid_argument(
                "row_labels/col_labels lengths must match freq_matrix dimensions.");
        }
        Tensor<std::string_view, 1> labels({freq_ncol}, kUninitialized);
        for (std::size_t i = 0; i < freq_ncol; ++i) {
            labels(i) = std::string_view(col_names[i]);
        }

        // 3. metrics (length 4)
        auto metrics = PyConvert::ToUIntVector(metrics_obj, "metrics");
        if (metrics.size() != 4) {
            throw std::invalid_argument("'metrics' must have exactly length 4.");
        }
        MetricsConfig config;
        for (std::size_t i = 0; i < 4; ++i) {
            config.Set(i, metrics[i] != 0);
        }
        if (config.count == 0) {
            throw std::invalid_argument("at least one metric must be selected.");
        }
        if (config.use_mrp && wraps.size() > 1) {
            throw std::invalid_argument("the 'MRP' metric is not supported with multiple POSets.");
        }

        // 4. subpopulation_count / total_bins
        Tensor<double, 1> subpopulation_count = [&]() {
            if (subpop_obj == Py_None) return Tensor<double, 1>({0}, 0.0);
            auto v = PyConvert::ToDoubleVector(subpop_obj, "subpopulation_count");
            Tensor<double, 1> t({v.size()}, kUninitialized);
            for (std::size_t i = 0; i < v.size(); ++i) t(i) = v[i];
            return t;
        }();
        if (config.use_inferential) {
            if (subpop_obj == Py_None) {
                throw std::invalid_argument(
                    "subpopulation_count is required for MannWhitneyInferentialDominance.");
            }
            if (subpopulation_count.size() != freq_ncol) {
                throw std::invalid_argument(
                    "subpopulation_count size must match the number of freq_matrix columns.");
            }
        }
        const bool needs_bins = config.use_mann_whitney || config.use_inferential;
        const auto total_bins_opt = PyConvert::ToOptionalUInt(bins_obj);
        if (needs_bins && (!total_bins_opt.has_value() || *total_bins_opt == 0)) {
            throw std::invalid_argument(
                "'total_bins' must be a positive integer for Mann-Whitney metrics.");
        }
        const std::uint64_t total_bins = total_bins_opt.value_or(0);

        // 5. count / seed / output / threads
        const bool has_le_matrix = (les_obj != Py_None);
        auto count = PyConvert::ToOptionalUInt(count_obj);
        auto output_interval = PyConvert::ToOptionalUInt(output_obj);
        const auto n_threads_opt = PyConvert::ToOptionalUInt(nthreads_obj);
        std::optional<std::uint64_t> seed;
        if (!has_le_matrix) {
            seed = PyConvert::ToOptionalSeed(seed_obj);
            if (seed.has_value() && !count.has_value()) {
                throw std::invalid_argument("'seed' can only be provided together with 'count'.");
            }
        }

        // 6. obs_vals mapping
        std::vector<std::uint32_t> obs_vals(freq_nrow);
        std::size_t total_poset_size = 1;
        if (wraps.size() > 1) {
            for (const auto* pw : wraps) {
                total_poset_size *= static_cast<std::size_t>(pw->GetPOSet()->size());
            }
            std::vector<std::string_view> tokens;
            for (std::size_t i = 0; i < freq_nrow; ++i) {
                SplitString(row_names[i], sep, tokens);
                if (tokens.size() != wraps.size()) {
                    throw MyException("A row label does not match the number of POSets.");
                }
                std::uint64_t global_id = 0, mult = 1;
                for (int j = static_cast<int>(wraps.size()) - 1; j >= 0; --j) {
                    global_id += wraps[j]->GetPOSet()->GetElementId(tokens[j]) * mult;
                    mult *= wraps[j]->GetPOSet()->size();
                }
                obs_vals[i] = static_cast<std::uint32_t>(global_id);
            }
        } else {
            total_poset_size = wraps.front()->GetPOSet()->size();
            for (std::size_t i = 0; i < freq_nrow; ++i) {
                obs_vals[i] =
                    static_cast<std::uint32_t>(wraps.front()->GetPOSet()->GetElementId(row_names[i]));
            }
        }

        // 7. Instantiate metrics
        auto [fles, func_names, eval_results] = InstantiateMetrics(
            config, total_poset_size, freq_matrix, labels, obs_vals, subpopulation_count,
            total_bins, wraps.size() == 1 ? wraps.front()->GetPOSet() : nullptr);

        // 8. Evaluation
        std::uint64_t le_count = 0;
        std::string le_type;
        PyObject* seed_used = nullptr;  // new ref or nullptr

        if (!has_le_matrix) {
            auto gen = MakeGenerator(wraps, count, seed);
            le_type = std::string(gen.name);
            std::visit([&](auto& leg) { leg.Start(count.value_or(0)); }, gen.instance);
            if (gen.seed.has_value()) {
                seed_used = PyCreate::FromSeed(*gen.seed);
            }
            auto display = output_interval.has_value()
                               ? std::unique_ptr<DisplayMessage>(std::make_unique<DisplayMessageStdout>(
                                     le_count, std::nullopt, output_interval.value(), "evaluated"))
                               : std::unique_ptr<DisplayMessage>(std::make_unique<DisplayMessageNull>());
            bool end_process = false;
            std::visit(
                [&](auto& leg) {
                    POSet::evaluation(fles, leg, eval_results, le_count, end_process, display.get(),
                                      POSet::EvaluationUpdateStrategy::Average);
                },
                gen.instance);
        } else {
            if (wraps.size() > 1) {
                throw std::invalid_argument("'linear_extensions' is not supported with multiple POSets.");
            }
            if (!count.has_value() || *count == 0) {
                throw std::invalid_argument(
                    "'count' (extensions per chain) is required and must be > 0 with 'linear_extensions'.");
            }
            POSet* poset = wraps.front()->GetPOSet();
            auto initial_les = ToLinearExtensions(les_obj, *poset);
            const std::uint64_t n_chains = initial_les.size();
            auto seeds = ToChainSeeds(seed_obj, n_chains);
            // Report seeds back as a list of strings.
            seed_used = PyList_New(static_cast<Py_ssize_t>(seeds.size()));
            for (std::size_t i = 0; i < seeds.size(); ++i) {
                PyList_SET_ITEM(seed_used, static_cast<Py_ssize_t>(i), PyCreate::FromSeed(seeds[i]));
            }
            ECSBubleyDyerFromLEMatrix source(poset, std::move(initial_les), std::move(seeds),
                                             *count, fles, eval_results);
            auto display = std::make_unique<DisplayMessageNull>();  // parallel path: quiet
            bool end_process = false;
            const std::uint64_t n_threads = n_threads_opt.value_or(0);
            std::string err;
            Py_BEGIN_ALLOW_THREADS
            try {
                POSet::evaluation_parallel(source, eval_results, le_count, end_process,
                                           display.get(), POSet::EvaluationUpdateStrategy::Average,
                                           n_threads);
            } catch (const std::exception& e) {
                err = e.what();
            } catch (...) {
                err = "unknown C++ exception in evaluation_parallel";
            }
            Py_END_ALLOW_THREADS
            if (!err.empty()) {
                Py_XDECREF(seed_used);
                throw std::runtime_error(err);
            }
            le_type = "BubleyDyerMultiChain";
        }

        // 9. FOD analysis per metric
        const std::uint64_t n_elements = static_cast<std::uint64_t>(freq_ncol * freq_ncol);
        PyRef top(PyDict_New());
        if (!top) {
            Py_XDECREF(seed_used);
            throw std::bad_alloc();
        }
        PyRef metric_list(PyList_New(static_cast<Py_ssize_t>(eval_results.size())));
        for (std::size_t k = 0; k < eval_results.size(); ++k) {
            const std::uint64_t metric = config.sequence[k];
            std::optional<Tensor<double, 2>> analyze;
            std::optional<Tensor<double, 2>> fod;
            if (metric == 3) {
                // MRP: only binMatrix.
            } else if (metric == 0) {
                analyze = *eval_results[k];
                fod = MinMaxTransitiveClosure(*eval_results[k]);
            } else {
                Tensor<double, 2> a(std::array<std::uint64_t, 2>{freq_ncol, freq_ncol}, kUninitialized);
                const auto& src = *eval_results[k];
                for (std::uint64_t idx = 0; idx < n_elements; ++idx) {
                    a.data()[idx] = src(idx, total_bins);
                }
                fod = MinMaxTransitiveClosure(a);
                analyze = std::move(a);
            }
            PyObject* entry = MetricToPy(fod, analyze, *eval_results[k], fles[k].get());
            PyCreate::DictSetSteal(entry, "name", PyCreate::FromString(func_names[k]));
            PyList_SET_ITEM(metric_list.get(), static_cast<Py_ssize_t>(k), entry);
        }

        PyCreate::DictSetSteal(top.get(), "metrics", metric_list.release());
        PyCreate::DictSetSteal(top.get(), "LEType", PyCreate::FromString(le_type));
        PyCreate::DictSetSteal(top.get(), "le_count",
                               PyLong_FromUnsignedLongLong(static_cast<unsigned long long>(le_count)));
        PyCreate::DictSetSteal(top.get(), "labels", PyCreate::FromStringVector(col_names));
        if (seed_used != nullptr) {
            PyCreate::DictSetSteal(top.get(), "seed", seed_used);
        } else {
            PyCreate::DictSetSteal(top.get(), "seed", (Py_INCREF(Py_None), Py_None));
        }
        return top.release();
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}
