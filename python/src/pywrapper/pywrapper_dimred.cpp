/**
 * @file pywrapper_dimred.cpp
 * @brief Entry-point Python per la riduzione dimensionale dei profili.
 *
 * @details Speculare a RunDimensionalityReduction / RunBidimentionalPosetRepresentation
 * di rwrapper_separation.cpp. Il calcolo (multi-thread) del core viene eseguito
 * con il GIL rilasciato; il reporter di progresso usa std::fprintf (thread-safe,
 * indipendente dal GIL).
 */

#include "pywrapper.h"
#include "py_common.h"
#include "py_convert.h"
#include "py_display.h"

#include "dimensionality_reduction.h"
#include "display_message.h"

#include <atomic>
#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace {

/// Legge una matrice profili (righe = profili, colonne = variabili, valori 0/1)
/// e produce la chiave bitwise per riga (prima colonna = bit piu' significativo).
std::vector<std::uint64_t> ProfileKeys(PyObject* profile_obj, std::size_t& out_ncol) {
    PyRef rows(PySequence_Fast(profile_obj, "profile must be a sequence of rows"));
    if (!rows) {
        throw std::invalid_argument("profile must be a matrix (sequence of rows).");
    }
    const Py_ssize_t nrow = PySequence_Fast_GET_SIZE(rows.get());
    std::vector<std::uint64_t> keys(static_cast<std::size_t>(nrow), 0);
    out_ncol = 0;
    for (Py_ssize_t i = 0; i < nrow; ++i) {
        PyRef row(PySequence_Fast(PySequence_Fast_GET_ITEM(rows.get(), i), "row must be a sequence"));
        if (!row) {
            throw std::invalid_argument("profile rows must be sequences.");
        }
        const std::size_t ncol = static_cast<std::size_t>(PySequence_Fast_GET_SIZE(row.get()));
        if (i == 0) {
            out_ncol = ncol;
        } else if (ncol != out_ncol) {
            throw std::invalid_argument("profile must be rectangular.");
        }
        std::uint64_t key = 0;
        for (std::size_t j = 0; j < ncol; ++j) {
            const long v = PyLong_AsLong(PySequence_Fast_GET_ITEM(row.get(), static_cast<Py_ssize_t>(j)));
            if (v == -1 && PyErr_Occurred()) {
                PyErr_Clear();
                throw std::invalid_argument("profile entries must be integers (0/1).");
            }
            key = (key << 1) | static_cast<std::uint64_t>(v != 0 ? 1 : 0);
        }
        keys[static_cast<std::size_t>(i)] = key;
    }
    return keys;
}

std::unordered_map<std::uint64_t, double> WeightsMap(const std::vector<std::uint64_t>& keys,
                                                     PyObject* weights_obj) {
    const auto w = [&]() {
        PyRef fast(PySequence_Fast(weights_obj, "weights must be a sequence"));
        if (!fast) {
            throw std::invalid_argument("weights must be a sequence of numbers.");
        }
        std::vector<double> out;
        const Py_ssize_t n = PySequence_Fast_GET_SIZE(fast.get());
        out.reserve(static_cast<std::size_t>(n));
        for (Py_ssize_t i = 0; i < n; ++i) {
            out.push_back(PyConvert::ToDoubleScalar(PySequence_Fast_GET_ITEM(fast.get(), i), "weight"));
        }
        return out;
    }();
    if (w.size() != keys.size()) {
        throw std::invalid_argument("weights must have one entry per profile row.");
    }
    std::unordered_map<std::uint64_t, double> weights;
    weights.reserve(keys.size());
    for (std::size_t i = 0; i < keys.size(); ++i) {
        weights[keys[i]] = w[i];
    }
    return weights;
}

/// Converte la map {ProfileID -> (rank, rankInv, weight, error)} in dict con
/// liste parallele profiles/x/y/weights/error.
PyObject* RepresentationToPy(
    const std::map<std::uint64_t, std::tuple<std::uint64_t, std::uint64_t, double, double>>& best) {
    const Py_ssize_t n = static_cast<Py_ssize_t>(best.size());
    PyRef profiles(PyList_New(n)), x(PyList_New(n)), y(PyList_New(n)),
        weights(PyList_New(n)), error(PyList_New(n));
    if (!profiles || !x || !y || !weights || !error) {
        throw std::bad_alloc();
    }
    Py_ssize_t idx = 0;
    for (const auto& [prof_id, tup] : best) {
        PyList_SET_ITEM(profiles.get(), idx, PyLong_FromUnsignedLongLong(static_cast<unsigned long long>(prof_id)));
        PyList_SET_ITEM(x.get(), idx, PyLong_FromUnsignedLongLong(static_cast<unsigned long long>(std::get<0>(tup))));
        PyList_SET_ITEM(y.get(), idx, PyLong_FromUnsignedLongLong(static_cast<unsigned long long>(std::get<1>(tup))));
        PyList_SET_ITEM(weights.get(), idx, PyFloat_FromDouble(std::get<2>(tup)));
        PyList_SET_ITEM(error.get(), idx, PyFloat_FromDouble(std::get<3>(tup)));
        ++idx;
    }
    PyRef dict(PyDict_New());
    if (!dict) {
        throw std::bad_alloc();
    }
    PyCreate::DictSetSteal(dict.get(), "profiles", profiles.release());
    PyCreate::DictSetSteal(dict.get(), "x", x.release());
    PyCreate::DictSetSteal(dict.get(), "y", y.release());
    PyCreate::DictSetSteal(dict.get(), "weights", weights.release());
    PyCreate::DictSetSteal(dict.get(), "error", error.release());
    return dict.release();
}

PyObject* IntListPlus1(const std::vector<std::uint64_t>& v) {
    PyRef list(PyList_New(static_cast<Py_ssize_t>(v.size())));
    if (!list) {
        throw std::bad_alloc();
    }
    for (std::size_t i = 0; i < v.size(); ++i) {
        PyList_SET_ITEM(list.get(), static_cast<Py_ssize_t>(i),
                        PyLong_FromUnsignedLongLong(static_cast<unsigned long long>(v[i] + 1)));
    }
    return list.release();
}

}  // namespace

// dimensionality_reduction(profile, weights, loss, lpom_strategy, output_every, thread_percentage)
//   -> {allLoss, variablesPriority, bestLossValue, bestVariablesPriority, bestRepresentation}
PyObject* pyx_dimensionality_reduction(PyObject* /*self*/, PyObject* args) {
    PyObject* profile_obj = nullptr;
    PyObject* weights_obj = nullptr;
    PyObject* loss_obj = nullptr;
    PyObject* lpom_obj = nullptr;
    PyObject* output_obj = nullptr;
    PyObject* thread_pct_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OOOOOO", &profile_obj, &weights_obj, &loss_obj, &lpom_obj,
                          &output_obj, &thread_pct_obj)) {
        return nullptr;
    }
    try {
        std::size_t ncol = 0;
        auto keys = ProfileKeys(profile_obj, ncol);
        auto weights = WeightsMap(keys, weights_obj);
        const std::string loss = PyConvert::ToString(loss_obj, "loss");
        const int lpom = static_cast<int>(PyLong_AsLong(lpom_obj));
        if (lpom == -1 && PyErr_Occurred()) {
            PyErr_Clear();
            throw std::invalid_argument("lpom_strategy must be an integer (0 or 1).");
        }
        const auto output_every = PyConvert::ToOptionalUInt(output_obj);
        const double thread_pct = PyConvert::ToDoubleScalar(thread_pct_obj, "thread_percentage");

        std::atomic<std::uint64_t> all_le_count{0};
        std::unique_ptr<DisplayMessage> display;
        if (!output_every.has_value()) {
            display = std::make_unique<DisplayMessageNull>();
        } else {
            display = std::make_unique<DisplayMessageAtomicStdout>(
                all_le_count, std::nullopt, output_every.value());
        }

        DimensionalityReductionResult result(all_le_count);
        std::string err;
        Py_BEGIN_ALLOW_THREADS
        try {
            ExactDimensionalityReduction(weights, static_cast<std::uint64_t>(ncol),
                                         std::string_view(loss), lpom, display.get(),
                                         thread_pct, result);
        } catch (const std::exception& e) {
            err = e.what();
        } catch (...) {
            err = "unknown C++ exception in dimensionality reduction";
        }
        Py_END_ALLOW_THREADS
        if (!err.empty()) {
            throw std::runtime_error(err);
        }

        auto& le = result.GetLeElaborated();
        auto& loss_vals = result.GetLossValues();

        // allLoss
        PyRef all_loss(PyList_New(static_cast<Py_ssize_t>(loss_vals.size())));
        if (!all_loss) {
            throw std::bad_alloc();
        }
        for (std::size_t i = 0; i < loss_vals.size(); ++i) {
            PyList_SET_ITEM(all_loss.get(), static_cast<Py_ssize_t>(i), PyFloat_FromDouble(loss_vals[i]));
        }

        // variablesPriority: list of permutations (each list of int, 1-indexed)
        PyRef var_prio(PyList_New(static_cast<Py_ssize_t>(le.size())));
        if (!var_prio) {
            throw std::bad_alloc();
        }
        for (std::size_t i = 0; i < le.size(); ++i) {
            PyList_SET_ITEM(var_prio.get(), static_cast<Py_ssize_t>(i), IntListPlus1(le[i]));
        }

        PyRef top(PyDict_New());
        if (!top) {
            throw std::bad_alloc();
        }
        PyCreate::DictSetSteal(top.get(), "allLoss", all_loss.release());
        PyCreate::DictSetSteal(top.get(), "variablesPriority", var_prio.release());
        PyCreate::DictSetSteal(top.get(), "bestLossValue", PyFloat_FromDouble(result.GetBestLoss()));
        PyCreate::DictSetSteal(top.get(), "bestVariablesPriority", IntListPlus1(result.GetBestLe()));
        PyCreate::DictSetSteal(top.get(), "bestRepresentation",
                               RepresentationToPy(result.GetBestProfileResults()));
        return top.release();
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

// bidimensional_poset_representation(profile, weights, loss, lpom_strategy, variable_priority)
//   -> {lossValue, variablesPriority, representation}
PyObject* pyx_bidimensional_poset_representation(PyObject* /*self*/, PyObject* args) {
    PyObject* profile_obj = nullptr;
    PyObject* weights_obj = nullptr;
    PyObject* loss_obj = nullptr;
    PyObject* lpom_obj = nullptr;
    PyObject* var_prio_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OOOOO", &profile_obj, &weights_obj, &loss_obj, &lpom_obj,
                          &var_prio_obj)) {
        return nullptr;
    }
    try {
        std::size_t ncol = 0;
        auto keys = ProfileKeys(profile_obj, ncol);
        auto weights = WeightsMap(keys, weights_obj);
        const std::string loss = PyConvert::ToString(loss_obj, "loss");
        const int lpom = static_cast<int>(PyLong_AsLong(lpom_obj));
        if (lpom == -1 && PyErr_Occurred()) {
            PyErr_Clear();
            throw std::invalid_argument("lpom_strategy must be an integer (0 or 1).");
        }

        // variable_priority (1-indexed from the user) -> 0-indexed for the core.
        auto prio_in = PyConvert::ToUIntVector(var_prio_obj, "variable_priority");
        std::vector<std::uint64_t> variable_priority;
        variable_priority.reserve(prio_in.size());
        for (const auto v : prio_in) {
            if (v == 0) {
                throw std::invalid_argument("variable_priority must be 1-indexed (>= 1).");
            }
            variable_priority.push_back(v - 1);
        }

        std::atomic<std::uint64_t> all_le_count{0};
        DimensionalityReductionResult result(all_le_count);
        std::string err;
        Py_BEGIN_ALLOW_THREADS
        try {
            BidimentionalPosetRepresentation(weights, static_cast<std::uint64_t>(ncol),
                                             std::string_view(loss), lpom, variable_priority,
                                             result);
        } catch (const std::exception& e) {
            err = e.what();
        } catch (...) {
            err = "unknown C++ exception in bidimensional representation";
        }
        Py_END_ALLOW_THREADS
        if (!err.empty()) {
            throw std::runtime_error(err);
        }

        PyRef top(PyDict_New());
        if (!top) {
            throw std::bad_alloc();
        }
        PyCreate::DictSetSteal(top.get(), "lossValue", PyFloat_FromDouble(result.GetBestLoss()));
        PyCreate::DictSetSteal(top.get(), "variablesPriority", IntListPlus1(result.GetBestLe()));
        PyCreate::DictSetSteal(top.get(), "representation",
                               RepresentationToPy(result.GetBestProfileResults()));
        return top.release();
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}
