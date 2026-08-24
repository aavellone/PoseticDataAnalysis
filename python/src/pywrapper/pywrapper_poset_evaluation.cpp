/**
 * @file pywrapper_poset_evaluation.cpp
 * @brief Entry-point Python per le valutazioni sulle estensioni lineari.
 *
 * @details Speculare a rwrapper_poset_evaluation.cpp:
 *  - ExactMRP / BubleyDyerMRP (Mutual Ranking Probability, esatta e MCMC);
 *  - ExactEvaluation / BubleyDyerEvaluation (metriche multiple);
 *  - BLSDominance (Bruggemann-Lerche-Sorensen).
 * Reporting su stdout (py_display.h). Seed a 64 bit restituiti come stringa.
 */

#include "pywrapper.h"
#include "py_common.h"
#include "py_convert.h"
#include "py_display.h"
#include "py_linear_generator_wrapper.h"

#include "function_linear_extension.h"
#include "function_linear_extension_average_height.h"
#include "function_linear_extension_mutual_ranking_probability.h"
#include "function_linear_extension_separation_asymmetric_lower.h"
#include "function_linear_extension_separation_asymmetric_upper.h"
#include "function_linear_extension_separation_symmetric.h"
#include "linear_extension_generator.h"
#include "linear_extension_generator_binary_variable.h"
#include "linear_extension_generator_tree_of_ideals.h"
#include "my_exception.h"
#include "poset.h"
#include "poset_wrapper.h"
#include "random.h"
#include "tensor.h"

#include <array>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace {

using FLEPtr = std::unique_ptr<FunctionLinearExtension>;
using MatPtr = std::unique_ptr<Tensor<double, 2>>;

/// Costruisce una FLE C++ interna dal nome della metrica (nessuna funzione R).
FLEPtr MakeFle(const std::string& name, POSet* poset, std::string& display_name) {
    auto it = POSetWrap::kFunctionLinearMapType.find(name);
    if (it == POSetWrap::kFunctionLinearMapType.end()) {
        throw MyException("Unknown metric name: '" + name +
                          "'. Valid names: MutualRankingProbability, AverageHeight, "
                          "symmetric, asymmetricLower, asymmetricUpper.");
    }
    switch (it->second) {
        case POSetWrap::FunctionLinearType::kMutualRankingProbability:
            display_name = "MutualRankingProbability";
            return std::make_unique<FLEMutualRankingProbability>(poset);
        case POSetWrap::FunctionLinearType::kAverageHeight:
            display_name = "AverageHeight";
            return std::make_unique<FLEAverageHeight>(poset);
        case POSetWrap::FunctionLinearType::kSeparationAsymmetricLower:
            display_name = "asymmetricLower";
            return std::make_unique<FLESeparationAsymmetricLower>(poset);
        case POSetWrap::FunctionLinearType::kSeparationAsymmetricUpper:
            display_name = "asymmetricUpper";
            return std::make_unique<FLESeparationAsymmetricUpper>(poset);
        case POSetWrap::FunctionLinearType::kSeparationSymmetric:
            display_name = "symmetric";
            return std::make_unique<FLESeparationSymmetric>(poset);
        case POSetWrap::FunctionLinearType::kRFunction:
            throw MyException(
                "Custom (RFunction) metrics are not supported in the Python binding.");
        default:
            throw MyException("Metric '" + name + "' is not allowed here.");
    }
}

/// Reporter: DisplayMessageNull se intervallo assente, altrimenti su stdout.
std::unique_ptr<DisplayMessage> MakeDisplay(std::uint64_t& le_count,
                                            std::optional<std::uint64_t> total,
                                            std::optional<std::uint64_t> output_every) {
    if (!output_every.has_value()) {
        return std::make_unique<DisplayMessageNull>();
    }
    return std::make_unique<DisplayMessageStdout>(le_count, total, output_every.value(),
                                                  "evaluated");
}

/// Costruisce list[dict{name,matrix,rows,cols}] dai risultati di valutazione.
PyObject* EvalResultsToPy(const std::vector<MatPtr>& results,
                          const std::vector<FLEPtr>& fles,
                          const std::vector<std::string>& names,
                          std::uint64_t total_le) {
    PyRef top(PyDict_New());
    if (!top) {
        throw std::bad_alloc();
    }
    PyRef list(PyList_New(static_cast<Py_ssize_t>(results.size())));
    if (!list) {
        throw std::bad_alloc();
    }
    for (std::size_t k = 0; k < results.size(); ++k) {
        FunctionLinearExtension* fle = fles[k].get();
        PyObject* entry = PyCreate::LabeledMatrix(
            *results[k],
            [fle](std::size_t i) { return fle->GetRowNameAt(i); },
            [fle](std::size_t j) { return fle->GetColNameAt(j); });
        // aggiunge il nome della metrica al dict
        PyCreate::DictSetSteal(entry, "name", PyCreate::FromString(names[k]));
        PyList_SET_ITEM(list.get(), static_cast<Py_ssize_t>(k), entry);  // ruba entry
    }
    PyCreate::DictSetSteal(top.get(), "results", list.release());
    PyCreate::DictSetSteal(top.get(), "n",
                           PyLong_FromUnsignedLongLong(static_cast<unsigned long long>(total_le)));
    return top.release();
}

/// Costruisce dict{mrp,elements,n} per i risultati MRP.
PyObject* MrpToPy(const Tensor<double, 2>& mrp, POSet* poset, std::uint64_t n_le) {
    PyRef dict(PyDict_New());
    if (!dict) {
        throw std::bad_alloc();
    }
    PyCreate::DictSetSteal(dict.get(), "mrp", PyCreate::FromDoubleMatrix(mrp));

    const std::size_t n = static_cast<std::size_t>(poset->size());
    PyRef els(PyList_New(static_cast<Py_ssize_t>(n)));
    if (!els) {
        throw std::bad_alloc();
    }
    for (std::size_t i = 0; i < n; ++i) {
        PyObject* s = PyCreate::FromString(poset->GetElementName(i));
        if (s == nullptr) {
            throw std::runtime_error("failed to allocate element name.");
        }
        PyList_SET_ITEM(els.get(), static_cast<Py_ssize_t>(i), s);
    }
    PyCreate::DictSetSteal(dict.get(), "elements", els.release());
    PyCreate::DictSetSteal(dict.get(), "n",
                           PyLong_FromUnsignedLongLong(static_cast<unsigned long long>(n_le)));
    return dict.release();
}

}  // namespace

// exact_mrp(poset, output_every) -> dict{mrp, elements, n}
PyObject* pyx_exact_mrp(PyObject* /*self*/, PyObject* args) {
    PyObject* poset_obj = nullptr;
    PyObject* output_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OO", &poset_obj, &output_obj)) {
        return nullptr;
    }
    try {
        const POSetWrap* poset_wrap = GetPOSetWrap(poset_obj);
        POSet* poset = poset_wrap->GetPOSet();
        const auto output_every = PyConvert::ToOptionalUInt(output_obj);

        std::unique_ptr<LinearExtensionGenerator> legenerator;
        if (poset_wrap->GetType() == POSetWrap::PosetType::kBinaryVariable) {
            legenerator = std::make_unique<LEGBinaryVariable>(poset->size());
        } else {
            auto lattice = poset->GetLatticeOfIdeals();
            legenerator = std::make_unique<LEGTreeOfIdeals>(poset->size(), *lattice);
        }
        legenerator->Start(0);

        const std::uint64_t n = poset->size();
        std::vector<MatPtr> eval_results;
        eval_results.push_back(
            std::make_unique<Tensor<double, 2>>(std::array<std::uint64_t, 2>{n, n}, 0.0));
        std::vector<FLEPtr> fles;
        fles.push_back(std::make_unique<FLEMutualRankingProbability>(poset));

        std::uint64_t le_count = 0;
        bool end_process = false;
        auto display = MakeDisplay(le_count, std::nullopt, output_every);
        POSet::evaluation(fles, *legenerator, eval_results, le_count, end_process,
                          display.get(), POSet::EvaluationUpdateStrategy::Average);

        return MrpToPy(*eval_results[0], poset, le_count);
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

// build_bubley_dyer_mrp_generator(poset, seed) -> (capsule, seed_str)
PyObject* pyx_build_bubley_dyer_mrp_generator(PyObject* /*self*/, PyObject* args) {
    PyObject* poset_obj = nullptr;
    PyObject* seed_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OO", &poset_obj, &seed_obj)) {
        return nullptr;
    }
    try {
        const POSetWrap* poset_wrap = GetPOSetWrap(poset_obj);
        std::optional<std::uint64_t> seed = PyConvert::ToOptionalSeed(seed_obj);
        if (!seed.has_value()) {
            seed = Random::GLOBAL.RndNextInt(0, std::numeric_limits<std::uint64_t>::max());
        }
        PyObject* capsule = MakeCapsule(
            BubleyDyerMRPGenerator::BuildBubleyDyerMRPGenerator(poset_wrap, seed.value()),
            kBDMRPCapsuleName);
        if (capsule == nullptr) {
            return nullptr;
        }
        PyObject* seed_str = PyCreate::FromSeed(seed.value());
        PyObject* tup = (seed_str != nullptr) ? PyTuple_New(2) : nullptr;
        if (tup == nullptr) {
            Py_DECREF(capsule);
            Py_XDECREF(seed_str);
            return nullptr;
        }
        PyTuple_SET_ITEM(tup, 0, capsule);
        PyTuple_SET_ITEM(tup, 1, seed_str);
        return tup;
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

// bubley_dyer_mrp(gen, quante, errore, output_every) -> dict{mrp, elements, n}
PyObject* pyx_bubley_dyer_mrp(PyObject* /*self*/, PyObject* args) {
    PyObject* gen_obj = nullptr;
    PyObject* quante_obj = nullptr;
    PyObject* errore_obj = nullptr;
    PyObject* output_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OOOO", &gen_obj, &quante_obj, &errore_obj, &output_obj)) {
        return nullptr;
    }
    try {
        auto* bd = GetCapsule<BubleyDyerMRPGenerator>(gen_obj, kBDMRPCapsuleName,
                                                      "Bubley-Dyer MRP generator");
        auto quante = PyConvert::ToOptionalUInt(quante_obj);
        const auto errore = PyConvert::ToOptionalDouble(errore_obj);
        const auto output_every = PyConvert::ToOptionalUInt(output_obj);

        if (!errore.has_value() && !quante.has_value()) {
            throw MyException("bubley_dyer_mrp: provide either 'n' (count) or a target 'error'.");
        }
        if (errore.has_value() && quante.has_value()) {
            PyDisplay::Warning("bubley_dyer_mrp: 'error' is ignored when 'n' is provided.");
        }

        POSet* poset = bd->poset_;
        bool updated = true;
        if (!bd->used_) {
            if (!quante.has_value()) {
                quante = bd->le_generator_->EvaluateNumberOfIteration(errore.value());
            }
            bd->le_generator_->Start(quante.value());
        } else {
            updated = bd->le_generator_->UpdateCounters(quante, errore);
            if (!updated) {
                PyDisplay::Line("bubley_dyer_mrp: desired error already reached — no new extensions.");
            } else {
                bd->le_generator_->Next();
            }
        }

        std::uint64_t total_le = bd->le_generator_->CurrentNumberOfLe() - 1;
        if (updated) {
            bd->used_ = true;
            std::vector<MatPtr> eval_results;
            eval_results.push_back(std::move(bd->mrp_));
            std::vector<FLEPtr> fles;
            fles.push_back(std::make_unique<FLEMutualRankingProbability>(poset));

            std::uint64_t le_count = 0;
            bool end_process = false;
            const std::uint64_t total_ext =
                (bd->le_generator_->NumberOfLe() + 1) - bd->le_generator_->CurrentNumberOfLe() + 1;
            auto display = MakeDisplay(le_count, total_ext, output_every);
            POSet::evaluation(fles, *(bd->le_generator_), eval_results, le_count, end_process,
                              display.get(), POSet::EvaluationUpdateStrategy::Average);
            bd->mrp_ = std::move(eval_results[0]);
            total_le += le_count;
        }

        return MrpToPy(*(bd->mrp_), poset, total_le);
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

// build_bubley_dyer_evaluation_generator(poset, functions, seed) -> (capsule, seed_str)
PyObject* pyx_build_bubley_dyer_evaluation_generator(PyObject* /*self*/, PyObject* args) {
    PyObject* poset_obj = nullptr;
    PyObject* functions_obj = nullptr;
    PyObject* seed_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OOO", &poset_obj, &functions_obj, &seed_obj)) {
        return nullptr;
    }
    try {
        const POSetWrap* poset_wrap = GetPOSetWrap(poset_obj);
        std::optional<std::uint64_t> seed = PyConvert::ToOptionalSeed(seed_obj);
        if (!seed.has_value()) {
            seed = Random::GLOBAL.RndNextInt(0, std::numeric_limits<std::uint64_t>::max());
        }
        auto functions = PyConvert::ToStringVector(functions_obj, "functions");
        PyObject* capsule = MakeCapsule(
            BubleyDyerEvaluationGenerator::BuildBubleyDyerEvaluationGenerator(
                poset_wrap, seed.value(), functions),
            kBDEvalCapsuleName);
        if (capsule == nullptr) {
            return nullptr;
        }
        PyObject* seed_str = PyCreate::FromSeed(seed.value());
        PyObject* tup = (seed_str != nullptr) ? PyTuple_New(2) : nullptr;
        if (tup == nullptr) {
            Py_DECREF(capsule);
            Py_XDECREF(seed_str);
            return nullptr;
        }
        PyTuple_SET_ITEM(tup, 0, capsule);
        PyTuple_SET_ITEM(tup, 1, seed_str);
        return tup;
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

// bubley_dyer_evaluation(gen, quante, errore, output_every) -> dict{results, n}
PyObject* pyx_bubley_dyer_evaluation(PyObject* /*self*/, PyObject* args) {
    PyObject* gen_obj = nullptr;
    PyObject* quante_obj = nullptr;
    PyObject* errore_obj = nullptr;
    PyObject* output_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OOOO", &gen_obj, &quante_obj, &errore_obj, &output_obj)) {
        return nullptr;
    }
    try {
        auto* gen = GetCapsule<BubleyDyerEvaluationGenerator>(gen_obj, kBDEvalCapsuleName,
                                                              "Bubley-Dyer evaluation generator");
        auto quante = PyConvert::ToOptionalUInt(quante_obj);
        const auto errore = PyConvert::ToOptionalDouble(errore_obj);
        const auto output_every = PyConvert::ToOptionalUInt(output_obj);

        if (!errore.has_value() && !quante.has_value()) {
            throw MyException("bubley_dyer_evaluation: provide either 'n' (count) or a target 'error'.");
        }
        if (errore.has_value() && quante.has_value()) {
            PyDisplay::Warning("bubley_dyer_evaluation: 'error' is ignored when 'n' is provided.");
        }

        bool updated = true;
        if (!gen->used_) {
            if (!quante.has_value()) {
                quante = gen->le_generator_->EvaluateNumberOfIteration(errore.value());
            }
            gen->le_generator_->Start(quante.value());
        } else {
            updated = gen->le_generator_->UpdateCounters(quante, errore);
            if (!updated) {
                PyDisplay::Line("bubley_dyer_evaluation: desired error already reached — no new extensions.");
            } else {
                gen->le_generator_->Next();
            }
        }

        std::uint64_t total_le = gen->le_generator_->CurrentNumberOfLe() - 1;
        if (updated) {
            gen->used_ = true;
            std::uint64_t le_count = 0;
            bool end_process = false;
            const std::uint64_t total_ext =
                (gen->le_generator_->NumberOfLe() + 1) - gen->le_generator_->CurrentNumberOfLe();
            auto display = MakeDisplay(le_count, total_ext, output_every);
            POSet::evaluation(gen->fles_, *(gen->le_generator_), gen->eval_results_, le_count,
                              end_process, display.get(), POSet::EvaluationUpdateStrategy::Average);
            total_le += le_count;
        }

        return EvalResultsToPy(gen->eval_results_, gen->fles_, gen->functions_name_, total_le);
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

// exact_evaluation(poset, functions, output_every) -> dict{results, n}
PyObject* pyx_exact_evaluation(PyObject* /*self*/, PyObject* args) {
    PyObject* poset_obj = nullptr;
    PyObject* functions_obj = nullptr;
    PyObject* output_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OOO", &poset_obj, &functions_obj, &output_obj)) {
        return nullptr;
    }
    try {
        const POSetWrap* poset_wrap = GetPOSetWrap(poset_obj);
        POSet* poset = poset_wrap->GetPOSet();
        const auto output_every = PyConvert::ToOptionalUInt(output_obj);
        auto names = PyConvert::ToStringVector(functions_obj, "functions");

        std::vector<FLEPtr> fles;
        std::vector<MatPtr> eval_results;
        std::vector<std::string> display_names;
        fles.reserve(names.size());
        eval_results.reserve(names.size());
        display_names.reserve(names.size());

        for (const auto& name : names) {
            std::string display_name;
            auto fle = MakeFle(name, poset, display_name);
            auto shape = fle->Shape();
            eval_results.push_back(std::make_unique<Tensor<double, 2>>(
                std::array<std::uint64_t, 2>{shape.at(0), shape.at(1)}, 0.0));
            fles.push_back(std::move(fle));
            display_names.push_back(std::move(display_name));
        }

        auto lattice = poset->GetLatticeOfIdeals();
        LEGTreeOfIdeals le_generator(poset->size(), *lattice);
        le_generator.Start(0);

        std::uint64_t le_count = 0;
        bool end_process = false;
        auto display = MakeDisplay(le_count, std::nullopt, output_every);
        POSet::evaluation(fles, le_generator, eval_results, le_count, end_process,
                          display.get(), POSet::EvaluationUpdateStrategy::Average);

        return EvalResultsToPy(eval_results, fles, display_names, le_count);
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

// bls_dominance(poset, relative) -> dict{matrix, elements}
PyObject* pyx_bls_dominance(PyObject* /*self*/, PyObject* args) {
    PyObject* poset_obj = nullptr;
    PyObject* relative_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OO", &poset_obj, &relative_obj)) {
        return nullptr;
    }
    try {
        POSet* poset = GetPOSetWrap(poset_obj)->GetPOSet();
        const bool relative = PyConvert::ToBool(relative_obj);
        const auto mat = relative ? poset->BLSDominanceRelative() : poset->BLSDominanceAbsolute();

        PyRef dict(PyDict_New());
        if (!dict) {
            throw std::bad_alloc();
        }
        PyCreate::DictSetSteal(dict.get(), "matrix", PyCreate::FromDoubleMatrix(mat));

        const std::size_t n = static_cast<std::size_t>(poset->size());
        PyRef els(PyList_New(static_cast<Py_ssize_t>(n)));
        if (!els) {
            throw std::bad_alloc();
        }
        for (std::size_t i = 0; i < n; ++i) {
            PyObject* s = PyCreate::FromString(poset->GetElementName(i));
            if (s == nullptr) {
                throw std::runtime_error("failed to allocate element name.");
            }
            PyList_SET_ITEM(els.get(), static_cast<Py_ssize_t>(i), s);
        }
        PyCreate::DictSetSteal(dict.get(), "elements", els.release());
        return dict.release();
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}
