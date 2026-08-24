/**
 * @file pywrapper_lex_fuzzy.cpp
 * @brief Entry-point Python per separazione lessicografica e misure fuzzy.
 *
 * @details Speculare a rwrapper_separation.cpp (parte lex/fuzzy). Usa le funzioni
 * del core R-free in separation.h. Sono supportate le norme "minimum" e
 * "product"; le norme come funzioni definite dall'utente (path R custom) non
 * sono esposte nel binding Python.
 */

#include "pywrapper.h"
#include "py_common.h"
#include "py_convert.h"

#include "separation.h"
#include "tensor.h"

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

/// Converte una lista Python di liste di str in vector<vector<string>>.
std::vector<std::vector<std::string>> ToModalita(PyObject* obj) {
    PyRef fast(PySequence_Fast(obj, "modalities must be a sequence of variables"));
    if (!fast) {
        throw std::invalid_argument("modalities must be a sequence (one entry per variable).");
    }
    const Py_ssize_t nv = PySequence_Fast_GET_SIZE(fast.get());
    std::vector<std::vector<std::string>> out;
    out.reserve(static_cast<std::size_t>(nv));
    for (Py_ssize_t i = 0; i < nv; ++i) {
        out.push_back(PyConvert::ToStringVector(PySequence_Fast_GET_ITEM(fast.get(), i),
                                                "modalities"));
    }
    return out;
}

/// Etichette dei profili: unisce i nomi di modalità con '_'.
PyObject* ProfileLabels(const std::vector<std::vector<std::uint64_t>>& profili,
                        const std::vector<std::vector<std::string>>& modalita) {
    PyRef list(PyList_New(static_cast<Py_ssize_t>(profili.size())));
    if (!list) {
        throw std::bad_alloc();
    }
    for (std::size_t pid = 0; pid < profili.size(); ++pid) {
        std::string label;
        const auto& profile = profili[pid];
        for (std::size_t vid = 0; vid < modalita.size(); ++vid) {
            label += modalita[vid][static_cast<std::size_t>(profile[vid])];
            if (vid + 1 < modalita.size()) {
                label += "_";
            }
        }
        PyObject* s = PyCreate::FromString(label);
        if (s == nullptr) {
            throw std::runtime_error("failed to build profile label.");
        }
        PyList_SET_ITEM(list.get(), static_cast<Py_ssize_t>(pid), s);
    }
    return list.release();
}

/// Verifica e restituisce i flag (do_all, do_lower, do_upper) dalla lista quali.
void ParseSepFlags(PyObject* quali_obj, bool& do_all, bool& do_lower, bool& do_upper) {
    do_all = do_lower = do_upper = false;
    const auto names = PyConvert::ToStringVector(quali_obj, "quali");
    for (const auto& name : names) {
        if (name == "symmetric") {
            do_all = true;
        } else if (name == "asymmetricLower") {
            do_lower = true;
        } else if (name == "asymmetricUpper") {
            do_upper = true;
        } else {
            throw std::invalid_argument(
                "unknown separation type '" + name +
                "' (use 'symmetric', 'asymmetricLower', 'asymmetricUpper').");
        }
    }
}

/// Costruisce una list 3D (n×n×n) da un tensore fuzzy in-betweenness.
PyObject* From3D(const std::vector<std::vector<std::vector<double>>>& t) {
    PyRef out(PyList_New(static_cast<Py_ssize_t>(t.size())));
    if (!out) {
        throw std::bad_alloc();
    }
    for (std::size_t i = 0; i < t.size(); ++i) {
        PyObject* plane = PyList_New(static_cast<Py_ssize_t>(t[i].size()));
        if (plane == nullptr) {
            throw std::bad_alloc();
        }
        for (std::size_t j = 0; j < t[i].size(); ++j) {
            PyObject* row = PyList_New(static_cast<Py_ssize_t>(t[i][j].size()));
            if (row == nullptr) {
                Py_DECREF(plane);
                throw std::bad_alloc();
            }
            for (std::size_t k = 0; k < t[i][j].size(); ++k) {
                PyList_SET_ITEM(row, static_cast<Py_ssize_t>(k),
                                PyFloat_FromDouble(t[i][j][k]));
            }
            PyList_SET_ITEM(plane, static_cast<Py_ssize_t>(j), row);
        }
        PyList_SET_ITEM(out.get(), static_cast<Py_ssize_t>(i), plane);
    }
    return out.release();
}

}  // namespace

// lex_separation(modalita) -> dict{symmetric, asymmetricLower, asymmetricUpper,
//                                  vertical, horizontal, labels}
PyObject* pyx_lex_separation(PyObject* /*self*/, PyObject* args) {
    PyObject* modalita_obj = nullptr;
    if (!PyArg_ParseTuple(args, "O", &modalita_obj)) {
        return nullptr;
    }
    try {
        auto modalita = ToModalita(modalita_obj);
        std::vector<std::uint64_t> counts;
        counts.reserve(modalita.size());
        bool all_equal = true;
        for (const auto& m : modalita) {
            counts.push_back(static_cast<std::uint64_t>(m.size()));
            if (!counts.empty() && counts.back() != counts.front()) {
                all_equal = false;
            }
        }

        LexSeparationResult result =
            (!all_equal)
                ? LexSeparationDeg(counts)
                : LexSeparationEqDeg(static_cast<std::uint64_t>(modalita.size()),
                                     counts.empty() ? 0 : counts.front());

        PyRef dict(PyDict_New());
        if (!dict) {
            throw std::bad_alloc();
        }
        PyCreate::DictSetSteal(dict.get(), "symmetric", PyCreate::FromDoubleMatrix(result.sep_all));
        PyCreate::DictSetSteal(dict.get(), "asymmetricLower", PyCreate::FromDoubleMatrix(result.sep_lower));
        PyCreate::DictSetSteal(dict.get(), "asymmetricUpper", PyCreate::FromDoubleMatrix(result.sep_upper));
        PyCreate::DictSetSteal(dict.get(), "vertical", PyCreate::FromDoubleMatrix(result.sep_vertical));
        PyCreate::DictSetSteal(dict.get(), "horizontal", PyCreate::FromDoubleMatrix(result.sep_horizontal));
        PyCreate::DictSetSteal(dict.get(), "labels", ProfileLabels(result.profili, modalita));
        return dict.release();
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

// lex_mrp(modalita) -> dict{mrp, labels}
PyObject* pyx_lex_mrp(PyObject* /*self*/, PyObject* args) {
    PyObject* modalita_obj = nullptr;
    if (!PyArg_ParseTuple(args, "O", &modalita_obj)) {
        return nullptr;
    }
    try {
        auto modalita = ToModalita(modalita_obj);
        std::vector<std::uint64_t> counts;
        counts.reserve(modalita.size());
        for (const auto& m : modalita) {
            counts.push_back(static_cast<std::uint64_t>(m.size()));
        }
        auto result = LexMrp(counts);

        PyRef dict(PyDict_New());
        if (!dict) {
            throw std::bad_alloc();
        }
        PyCreate::DictSetSteal(dict.get(), "mrp", PyCreate::FromDoubleMatrix(result.mrp));
        PyCreate::DictSetSteal(dict.get(), "labels", ProfileLabels(result.profili, modalita));
        return dict.release();
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

// fuzzy_separation(dominance, elements, norm, quali) -> dict{<name>: matrix, elements}
PyObject* pyx_fuzzy_separation(PyObject* /*self*/, PyObject* args) {
    PyObject* dom_obj = nullptr;
    PyObject* elems_obj = nullptr;
    PyObject* norm_obj = nullptr;
    PyObject* quali_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OOOO", &dom_obj, &elems_obj, &norm_obj, &quali_obj)) {
        return nullptr;
    }
    try {
        auto dominance = PyConvert::ToDoubleMatrix(dom_obj, "dominance");
        auto elements = PyConvert::ToStringVector(elems_obj, "elements");
        const std::string norm = PyConvert::ToString(norm_obj, "norm");
        bool do_all = false, do_lower = false, do_upper = false;
        ParseSepFlags(quali_obj, do_all, do_lower, do_upper);

        GeneralSeparationResult res =
            (norm == "minimum")
                ? GeneralSeparation(dominance, MinNormConorm{}, MaxNormConorm{}, do_all, do_lower, do_upper)
                : (norm == "product")
                      ? GeneralSeparation(dominance, ProdNormConorm{}, ProbNormConorm{}, do_all, do_lower, do_upper)
                      : throw std::invalid_argument("unknown norm '" + norm +
                                                    "' (use 'minimum' or 'product').");

        PyRef dict(PyDict_New());
        if (!dict) {
            throw std::bad_alloc();
        }
        if (do_all) {
            PyCreate::DictSetSteal(dict.get(), "symmetric", PyCreate::FromDoubleMatrix(res.sep_all));
        }
        if (do_lower) {
            PyCreate::DictSetSteal(dict.get(), "asymmetricLower", PyCreate::FromDoubleMatrix(res.sep_lower));
        }
        if (do_upper) {
            PyCreate::DictSetSteal(dict.get(), "asymmetricUpper", PyCreate::FromDoubleMatrix(res.sep_upper));
        }
        PyCreate::DictSetSteal(dict.get(), "elements", PyCreate::FromStringVector(elements));
        return dict.release();
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

// fuzzy_inbetweenness(dominance, elements, norm, quali) -> dict{<name>: 3d array, elements}
PyObject* pyx_fuzzy_inbetweenness(PyObject* /*self*/, PyObject* args) {
    PyObject* dom_obj = nullptr;
    PyObject* elems_obj = nullptr;
    PyObject* norm_obj = nullptr;
    PyObject* quali_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OOOO", &dom_obj, &elems_obj, &norm_obj, &quali_obj)) {
        return nullptr;
    }
    try {
        auto dominance = PyConvert::ToDoubleMatrix(dom_obj, "dominance");
        auto elements = PyConvert::ToStringVector(elems_obj, "elements");
        const std::string norm = PyConvert::ToString(norm_obj, "norm");
        bool do_all = false, do_lower = false, do_upper = false;
        ParseSepFlags(quali_obj, do_all, do_lower, do_upper);

        const std::size_t n = static_cast<std::size_t>(dominance.Extent(0));
        using Tensor3D = std::vector<std::vector<std::vector<double>>>;
        auto zeros = [n]() {
            return Tensor3D(n, std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0)));
        };
        Tensor3D inbet_all, inbet_lower, inbet_upper;
        if (do_all) inbet_all = zeros();
        if (do_lower) inbet_lower = zeros();
        if (do_upper) inbet_upper = zeros();

        auto compute = [&](auto norm_f, auto conorm_f) {
            for (std::size_t pi = 0; pi < n; ++pi) {
                for (std::size_t qi = pi + 1; qi < n; ++qi) {
                    for (std::size_t ri = 0; ri < n; ++ri) {
                        double finb_prq = 0.0, finb_qrp = 0.0, finbqrp = 0.0;
                        GeneralFuzzyInBetweenness(pi, qi, ri, dominance, norm_f, conorm_f,
                                                  finb_prq, finb_qrp, finbqrp);
                        if (do_lower) {
                            inbet_lower[pi][qi][ri] = finb_prq;
                            inbet_lower[qi][pi][ri] = finb_qrp;
                        }
                        if (do_upper) {
                            inbet_upper[pi][qi][ri] = finb_qrp;
                            inbet_upper[qi][pi][ri] = finb_prq;
                        }
                        if (do_all) {
                            inbet_all[pi][qi][ri] = finbqrp;
                            inbet_all[qi][pi][ri] = finbqrp;
                        }
                    }
                }
            }
        };

        if (norm == "minimum") {
            compute(MinNormConorm{}, MaxNormConorm{});
        } else if (norm == "product") {
            compute(ProdNormConorm{}, ProbNormConorm{});
        } else {
            throw std::invalid_argument("unknown norm '" + norm + "' (use 'minimum' or 'product').");
        }

        PyRef dict(PyDict_New());
        if (!dict) {
            throw std::bad_alloc();
        }
        if (do_all) PyCreate::DictSetSteal(dict.get(), "symmetric", From3D(inbet_all));
        if (do_lower) PyCreate::DictSetSteal(dict.get(), "asymmetricLower", From3D(inbet_lower));
        if (do_upper) PyCreate::DictSetSteal(dict.get(), "asymmetricUpper", From3D(inbet_upper));
        PyCreate::DictSetSteal(dict.get(), "elements", PyCreate::FromStringVector(elements));
        return dict.release();
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}
