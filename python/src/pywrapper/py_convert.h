/**
 * @file py_convert.h
 * @brief Conversioni Python(SEXP-equivalente: PyObject*) <-> C++20.
 *
 * @details Speculare a rwrapper_conversion.h del package R:
 *  - namespace @ref PyConvert — estrazione dati da Python verso C++ (nessuna
 *    ownership acquisita, solo copie o puntatori osservatori).
 *  - namespace @ref PyCreate  — costruzione di oggetti Python da tipi C++.
 *
 * Scelte di mapping (idiomatiche in Python, senza numpy):
 *  - vettore di stringhe R  -> list/tuple di str.
 *  - matrice Nx2 di archi   -> sequenza di coppie (from, to).
 *  - matrice numerica R     -> list di list (righe).
 *  - seed a 64 bit          -> str di cifre decimali (come nel package R),
 *    ma qui accettiamo anche int Python (arbitrariamente grande) senza perdita.
 */

#pragma once

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <array>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include "py_common.h"
#include "tensor.h"

namespace PyConvert {

/**
 * @brief Estrae una std::string UTF-8 da un oggetto Python str.
 * @throws std::invalid_argument se non e' uno str.
 */
[[nodiscard]] inline std::string ToString(PyObject* obj, const char* what = "value") {
    if (!PyUnicode_Check(obj)) {
        throw std::invalid_argument(std::string(what) + " must be a string.");
    }
    Py_ssize_t len = 0;
    const char* data = PyUnicode_AsUTF8AndSize(obj, &len);
    if (data == nullptr) {
        throw std::invalid_argument(std::string("could not decode ") + what + " as UTF-8.");
    }
    return std::string(data, static_cast<std::size_t>(len));
}

/**
 * @brief Converte una sequenza Python (list/tuple) di str in std::vector<std::string>.
 */
[[nodiscard]] inline std::vector<std::string> ToStringVector(PyObject* seq,
                                                             const char* what = "elements") {
    PyRef fast(PySequence_Fast(seq, "expected a sequence of strings"));
    if (!fast) {
        throw std::invalid_argument(std::string(what) + " must be a sequence of strings.");
    }
    const Py_ssize_t n = PySequence_Fast_GET_SIZE(fast.get());
    std::vector<std::string> out;
    out.reserve(static_cast<std::size_t>(n));
    for (Py_ssize_t i = 0; i < n; ++i) {
        PyObject* item = PySequence_Fast_GET_ITEM(fast.get(), i);  // borrowed
        out.push_back(ToString(item, what));
    }
    return out;
}

/**
 * @brief Converte una sequenza di coppie (from, to) di str in
 * std::vector<std::pair<std::string, std::string>>.
 * @details None o sequenza vuota -> vettore vuoto (poset senza relazioni).
 */
[[nodiscard]] inline std::vector<std::pair<std::string, std::string>>
ToPairVector(PyObject* seq) {
    if (seq == nullptr || seq == Py_None) {
        return {};
    }
    PyRef fast(PySequence_Fast(seq, "relations must be a sequence of (from, to) pairs"));
    if (!fast) {
        throw std::invalid_argument("relations must be a sequence of (from, to) pairs.");
    }
    const Py_ssize_t n = PySequence_Fast_GET_SIZE(fast.get());
    std::vector<std::pair<std::string, std::string>> out;
    out.reserve(static_cast<std::size_t>(n));
    for (Py_ssize_t i = 0; i < n; ++i) {
        PyObject* item = PySequence_Fast_GET_ITEM(fast.get(), i);  // borrowed
        PyRef pair(PySequence_Fast(item, "each relation must be a (from, to) pair"));
        if (!pair || PySequence_Fast_GET_SIZE(pair.get()) != 2) {
            throw std::invalid_argument("each relation must be a (from, to) pair of strings.");
        }
        std::string from = ToString(PySequence_Fast_GET_ITEM(pair.get(), 0), "relation endpoint");
        std::string to = ToString(PySequence_Fast_GET_ITEM(pair.get(), 1), "relation endpoint");
        out.emplace_back(std::move(from), std::move(to));
    }
    return out;
}

/**
 * @brief Costruisce una matrice di adiacenza da elementi + coppie (a, b) => a~b.
 * @details Analogo di RConvert::ToAdjacencyMatrix. @p out_elements riceve i nomi.
 * @throws std::invalid_argument su nome duplicato o arco con estremo sconosciuto.
 */
[[nodiscard]] inline Tensor<std::uint8_t, 2> ToAdjacencyMatrix(
    PyObject* elements_obj, PyObject* pairs_obj, std::vector<std::string>& out_elements) {
    out_elements = ToStringVector(elements_obj);
    const auto pairs = ToPairVector(pairs_obj);

    std::unordered_map<std::string, std::size_t> id;
    id.reserve(out_elements.size());
    for (std::size_t i = 0; i < out_elements.size(); ++i) {
        const auto [it, inserted] = id.emplace(out_elements[i], i);
        if (!inserted) {
            throw std::invalid_argument("duplicate element name '" + out_elements[i] + "'.");
        }
    }

    const std::size_t n = out_elements.size();
    Tensor<std::uint8_t, 2> adj(std::array<std::uint64_t, 2>{n, n}, 0);
    for (const auto& p : pairs) {
        auto i1 = id.find(p.first);
        auto i2 = id.find(p.second);
        if (i1 == id.end() || i2 == id.end()) {
            throw std::invalid_argument("relation refers to an unknown element.");
        }
        adj(i1->second, i2->second) = 1;
    }
    return adj;
}

/**
 * @brief Costruisce una matrice di adiacenza inferendo gli elementi dagli archi.
 * @details Analogo di RConvert::ToAdjacencyMatrixFromEdges.
 */
[[nodiscard]] inline Tensor<std::uint8_t, 2> ToAdjacencyMatrixFromEdges(
    PyObject* pairs_obj, std::vector<std::string>& out_elements) {
    const auto pairs = ToPairVector(pairs_obj);
    out_elements.clear();

    std::unordered_map<std::string, std::size_t> id;
    for (const auto& p : pairs) {
        if (id.try_emplace(p.first, out_elements.size()).second) {
            out_elements.push_back(p.first);
        }
        if (id.try_emplace(p.second, out_elements.size()).second) {
            out_elements.push_back(p.second);
        }
    }

    const std::size_t n = out_elements.size();
    Tensor<std::uint8_t, 2> adj(std::array<std::uint64_t, 2>{n, n}, 0);
    for (const auto& p : pairs) {
        adj(id.find(p.first)->second, id.find(p.second)->second) = 1;
    }
    return adj;
}

/**
 * @brief Estrae un vettore di puntatori osservatori POSetWrap da una sequenza
 * di capsule Python (analogo di RConvert::ToPOSetWrapVector).
 */
[[nodiscard]] inline std::vector<const POSetWrap*> ToPOSetWrapVector(PyObject* seq) {
    PyRef fast(PySequence_Fast(seq, "expected a sequence of POSet handles"));
    if (!fast) {
        throw std::invalid_argument("expected a sequence of POSet objects.");
    }
    const Py_ssize_t n = PySequence_Fast_GET_SIZE(fast.get());
    std::vector<const POSetWrap*> out;
    out.reserve(static_cast<std::size_t>(n));
    for (Py_ssize_t i = 0; i < n; ++i) {
        PyObject* item = PySequence_Fast_GET_ITEM(fast.get(), i);  // borrowed
        out.push_back(GetPOSetWrap(item));
    }
    return out;
}

/// Converte una sequenza Python di interi non negativi in vector<uint64>.
[[nodiscard]] inline std::vector<std::uint64_t> ToUIntVector(PyObject* seq, const char* what = "values") {
    PyRef fast(PySequence_Fast(seq, "expected a sequence of integers"));
    if (!fast) {
        throw std::invalid_argument(std::string(what) + " must be a sequence of integers.");
    }
    const Py_ssize_t n = PySequence_Fast_GET_SIZE(fast.get());
    std::vector<std::uint64_t> out;
    out.reserve(static_cast<std::size_t>(n));
    for (Py_ssize_t i = 0; i < n; ++i) {
        PyObject* item = PySequence_Fast_GET_ITEM(fast.get(), i);
        const unsigned long long v = PyLong_AsUnsignedLongLong(item);
        if (v == static_cast<unsigned long long>(-1) && PyErr_Occurred()) {
            PyErr_Clear();
            throw std::invalid_argument(std::string(what) + " must be non-negative integers.");
        }
        out.push_back(static_cast<std::uint64_t>(v));
    }
    return out;
}

/// Converte una sequenza Python di numeri in vector<double>.
[[nodiscard]] inline std::vector<double> ToDoubleVector(PyObject* seq, const char* what = "values") {
    PyRef fast(PySequence_Fast(seq, "expected a sequence of numbers"));
    if (!fast) {
        throw std::invalid_argument(std::string(what) + " must be a sequence of numbers.");
    }
    const Py_ssize_t n = PySequence_Fast_GET_SIZE(fast.get());
    std::vector<double> out;
    out.reserve(static_cast<std::size_t>(n));
    for (Py_ssize_t i = 0; i < n; ++i) {
        const double v = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(fast.get(), i));
        if (v == -1.0 && PyErr_Occurred()) {
            PyErr_Clear();
            throw std::invalid_argument(std::string(what) + " must be numbers.");
        }
        out.push_back(v);
    }
    return out;
}

/// Estrae un double da un oggetto Python (int o float).
[[nodiscard]] inline double ToDoubleScalar(PyObject* obj, const char* what = "value") {
    const double v = PyFloat_AsDouble(obj);
    if (v == -1.0 && PyErr_Occurred()) {
        PyErr_Clear();
        throw std::invalid_argument(std::string(what) + " must be a real number.");
    }
    return v;
}

/**
 * @brief Converte una matrice Python (sequenza di sequenze numeriche) in
 * Tensor<double, 2>. Verifica che sia rettangolare.
 */
[[nodiscard]] inline Tensor<double, 2> ToDoubleMatrix(PyObject* matrix_obj,
                                                      const char* what = "matrix") {
    PyRef rows(PySequence_Fast(matrix_obj, "expected a matrix (sequence of rows)"));
    if (!rows) {
        throw std::invalid_argument(std::string(what) + " must be a sequence of rows.");
    }
    const Py_ssize_t nrow = PySequence_Fast_GET_SIZE(rows.get());
    std::size_t ncol = 0;
    if (nrow > 0) {
        PyRef first(PySequence_Fast(PySequence_Fast_GET_ITEM(rows.get(), 0), "row must be a sequence"));
        if (!first) {
            throw std::invalid_argument(std::string(what) + " rows must be sequences.");
        }
        ncol = static_cast<std::size_t>(PySequence_Fast_GET_SIZE(first.get()));
    }

    Tensor<double, 2> out(std::array<std::uint64_t, 2>{static_cast<std::uint64_t>(nrow),
                                                       static_cast<std::uint64_t>(ncol)},
                          0.0);
    for (Py_ssize_t i = 0; i < nrow; ++i) {
        PyRef row(PySequence_Fast(PySequence_Fast_GET_ITEM(rows.get(), i), "row must be a sequence"));
        if (!row || static_cast<std::size_t>(PySequence_Fast_GET_SIZE(row.get())) != ncol) {
            throw std::invalid_argument(std::string(what) + " must be rectangular.");
        }
        for (std::size_t j = 0; j < ncol; ++j) {
            out(static_cast<std::uint64_t>(i), static_cast<std::uint64_t>(j)) =
                ToDoubleScalar(PySequence_Fast_GET_ITEM(row.get(), static_cast<Py_ssize_t>(j)),
                               "matrix entry");
        }
    }
    return out;
}

/**
 * @brief Estrae un bool da un oggetto Python (usa la truthiness Python).
 */
[[nodiscard]] inline bool ToBool(PyObject* obj) {
    const int v = PyObject_IsTrue(obj);
    if (v < 0) {
        throw std::invalid_argument("could not evaluate value as boolean.");
    }
    return v != 0;
}

/**
 * @brief Estrae un seed a 64 bit da un oggetto Python.
 * @details Canale ufficiale: str di cifre decimali (come nel package R). Qui
 * accettiamo anche int Python: essendo a precisione arbitraria non c'e' perdita
 * (a differenza di R, i cui interi sono a 32 bit). None -> nullopt (seed casuale).
 * @throws std::invalid_argument su tipo errato, negativo o fuori range uint64.
 */
[[nodiscard]] inline std::optional<std::uint64_t> ToOptionalSeed(PyObject* obj) {
    if (obj == nullptr || obj == Py_None) {
        return std::nullopt;
    }
    if (PyUnicode_Check(obj)) {
        const std::string s = ToString(obj, "seed");
        if (s.empty()) {
            throw std::invalid_argument("seed string must contain decimal digits.");
        }
        try {
            std::size_t pos = 0;
            const unsigned long long value = std::stoull(s, &pos, 10);
            if (pos != s.size()) {
                throw std::invalid_argument("seed must be a string of decimal digits.");
            }
            return static_cast<std::uint64_t>(value);
        } catch (const std::out_of_range&) {
            throw std::invalid_argument("seed is out of range for a 64-bit unsigned integer.");
        } catch (const std::invalid_argument&) {
            throw std::invalid_argument("seed must be a string of decimal digits.");
        }
    }
    if (PyLong_Check(obj)) {
        // PyLong_AsUnsignedLongLong raises OverflowError for negative values and
        // for values that do not fit in 64 bits: a single check covers both.
        const unsigned long long value = PyLong_AsUnsignedLongLong(obj);
        if (value == static_cast<unsigned long long>(-1) && PyErr_Occurred()) {
            PyErr_Clear();
            throw std::invalid_argument(
                "seed must be a non-negative integer that fits in 64 bits (0 .. 2**64-1).");
        }
        return static_cast<std::uint64_t>(value);
    }
    throw std::invalid_argument("seed must be an int, a decimal-digit string, or None.");
}

/// Converte una sequenza Python di seed (int o stringhe decimali) in vector<uint64>.
[[nodiscard]] inline std::vector<std::uint64_t> ToSeedVector(PyObject* seq) {
    PyRef fast(PySequence_Fast(seq, "expected a sequence of seeds"));
    if (!fast) {
        throw std::invalid_argument("seeds must be a sequence.");
    }
    const Py_ssize_t n = PySequence_Fast_GET_SIZE(fast.get());
    std::vector<std::uint64_t> out;
    out.reserve(static_cast<std::size_t>(n));
    for (Py_ssize_t i = 0; i < n; ++i) {
        auto s = ToOptionalSeed(PySequence_Fast_GET_ITEM(fast.get(), i));
        if (!s.has_value()) {
            throw std::invalid_argument("a per-chain seed cannot be None.");
        }
        out.push_back(s.value());
    }
    return out;
}

/**
 * @brief Estrae un uint64 opzionale (per parametri tipo 'count'/'linear_extensions').
 */
[[nodiscard]] inline std::optional<std::uint64_t> ToOptionalUInt(PyObject* obj) {
    if (obj == nullptr || obj == Py_None) {
        return std::nullopt;
    }
    if (!PyLong_Check(obj)) {
        throw std::invalid_argument("expected an integer or None.");
    }
    // PyLong_AsUnsignedLongLong raises OverflowError for negative and too-large
    // values alike.
    const unsigned long long value = PyLong_AsUnsignedLongLong(obj);
    if (value == static_cast<unsigned long long>(-1) && PyErr_Occurred()) {
        PyErr_Clear();
        throw std::invalid_argument(
            "value must be a non-negative integer that fits in 64 bits (0 .. 2**64-1).");
    }
    return static_cast<std::uint64_t>(value);
}

/**
 * @brief Estrae un double opzionale (per 'error'/'output_every_sec').
 */
[[nodiscard]] inline std::optional<double> ToOptionalDouble(PyObject* obj) {
    if (obj == nullptr || obj == Py_None) {
        return std::nullopt;
    }
    const double value = PyFloat_AsDouble(obj);
    if (value == -1.0 && PyErr_Occurred()) {
        PyErr_Clear();
        throw std::invalid_argument("expected a real number or None.");
    }
    return value;
}

}  // namespace PyConvert

// ===========================================================================
// PyCreate — costruzione di oggetti Python da tipi C++
// ===========================================================================

namespace PyCreate {

/// Crea uno str Python (UTF-8) da una string_view.
[[nodiscard]] inline PyObject* FromString(std::string_view s) {
    return PyUnicode_FromStringAndSize(s.data(), static_cast<Py_ssize_t>(s.size()));
}

/**
 * @brief Converte un range iterabile di id di elementi (bitset/set) in una
 * list Python di nomi. @p poset serve per la mappa id -> nome.
 */
template <typename Range, typename POSetT>
[[nodiscard]] inline PyObject* NamesFromIds(const Range& ids, const POSetT& poset,
                                            Py_ssize_t count) {
    PyRef list(PyList_New(count));
    if (!list) {
        throw std::bad_alloc();
    }
    Py_ssize_t i = 0;
    for (const auto id : ids) {
        std::string_view name = poset.GetElementName(static_cast<std::size_t>(id));
        PyObject* s = FromString(name);
        if (s == nullptr) {
            throw std::runtime_error("failed to allocate Python string.");
        }
        PyList_SET_ITEM(list.get(), i++, s);  // ruba il riferimento a s
    }
    return list.release();
}

/**
 * @brief Converte un vettore di stringhe C++ in una list Python di str.
 */
[[nodiscard]] inline PyObject* FromStringVector(const std::vector<std::string>& v) {
    PyRef list(PyList_New(static_cast<Py_ssize_t>(v.size())));
    if (!list) {
        throw std::bad_alloc();
    }
    for (std::size_t i = 0; i < v.size(); ++i) {
        PyObject* s = FromString(v[i]);
        if (s == nullptr) {
            throw std::runtime_error("failed to allocate Python string.");
        }
        PyList_SET_ITEM(list.get(), static_cast<Py_ssize_t>(i), s);
    }
    return list.release();
}

/**
 * @brief Converte un contenitore di coppie (id/id o nome/nome) in una list di
 * tuple Python (from, to). @p to_name mappa un elemento della coppia al nome.
 */
template <typename PairRange, typename ToName>
[[nodiscard]] inline PyObject* FromPairs(const PairRange& pairs, Py_ssize_t count,
                                         ToName&& to_name) {
    PyRef list(PyList_New(count));
    if (!list) {
        throw std::bad_alloc();
    }
    Py_ssize_t i = 0;
    for (const auto& p : pairs) {
        PyObject* a = FromString(to_name(p.first));
        PyObject* b = FromString(to_name(p.second));
        if (a == nullptr || b == nullptr) {
            Py_XDECREF(a);
            Py_XDECREF(b);
            throw std::runtime_error("failed to allocate Python string.");
        }
        PyObject* tup = PyTuple_New(2);
        if (tup == nullptr) {
            Py_DECREF(a);
            Py_DECREF(b);
            throw std::bad_alloc();
        }
        PyTuple_SET_ITEM(tup, 0, a);  // ruba i riferimenti
        PyTuple_SET_ITEM(tup, 1, b);
        PyList_SET_ITEM(list.get(), i++, tup);
    }
    return list.release();
}

/**
 * @brief Converte una matrice intera C++ (Tensor 2D) in una list-di-list Python.
 * @tparam MatrixT tipo con Extent(0/1) e operator()(i, j).
 */
template <typename MatrixT>
[[nodiscard]] inline PyObject* FromIntMatrix(const MatrixT& m) {
    const Py_ssize_t nrow = static_cast<Py_ssize_t>(m.Extent(0));
    const Py_ssize_t ncol = static_cast<Py_ssize_t>(m.Extent(1));
    PyRef rows(PyList_New(nrow));
    if (!rows) {
        throw std::bad_alloc();
    }
    for (Py_ssize_t i = 0; i < nrow; ++i) {
        PyObject* row = PyList_New(ncol);
        if (row == nullptr) {
            throw std::bad_alloc();
        }
        for (Py_ssize_t j = 0; j < ncol; ++j) {
            PyObject* v = PyLong_FromLong(static_cast<long>(
                m(static_cast<std::uint64_t>(i), static_cast<std::uint64_t>(j))));
            if (v == nullptr) {
                Py_DECREF(row);
                throw std::runtime_error("failed to allocate Python int.");
            }
            PyList_SET_ITEM(row, j, v);
        }
        PyList_SET_ITEM(rows.get(), i, row);  // ruba il riferimento a row
    }
    return rows.release();
}

/**
 * @brief Converte una matrice di double C++ (Tensor 2D) in una list-di-list Python.
 */
template <typename MatrixT>
[[nodiscard]] inline PyObject* FromDoubleMatrix(const MatrixT& m) {
    const Py_ssize_t nrow = static_cast<Py_ssize_t>(m.Extent(0));
    const Py_ssize_t ncol = static_cast<Py_ssize_t>(m.Extent(1));
    PyRef rows(PyList_New(nrow));
    if (!rows) {
        throw std::bad_alloc();
    }
    for (Py_ssize_t i = 0; i < nrow; ++i) {
        PyObject* row = PyList_New(ncol);
        if (row == nullptr) {
            throw std::bad_alloc();
        }
        for (Py_ssize_t j = 0; j < ncol; ++j) {
            PyObject* v = PyFloat_FromDouble(
                m(static_cast<std::uint64_t>(i), static_cast<std::uint64_t>(j)));
            if (v == nullptr) {
                Py_DECREF(row);
                throw std::runtime_error("failed to allocate Python float.");
            }
            PyList_SET_ITEM(row, j, v);
        }
        PyList_SET_ITEM(rows.get(), i, row);
    }
    return rows.release();
}

/// Converte un seed a 64 bit nella sua rappresentazione decimale (str Python).
[[nodiscard]] inline PyObject* FromSeed(std::uint64_t seed) {
    return PyUnicode_FromString(std::to_string(seed).c_str());
}

/**
 * @brief Inserisce (chiave -> valore) in un dict, consumando il riferimento a
 * @p value (nuovo riferimento) anche in caso di errore.
 * @throws std::bad_alloc / std::runtime_error su fallimento.
 */
inline void DictSetSteal(PyObject* dict, const char* key, PyObject* value) {
    if (value == nullptr) {
        throw std::runtime_error(std::string("failed to build value for key '") + key + "'.");
    }
    const int rc = PyDict_SetItemString(dict, key, value);
    Py_DECREF(value);  // PyDict_SetItemString incrementa il refcount del valore
    if (rc != 0) {
        throw std::runtime_error(std::string("failed to set dict key '") + key + "'.");
    }
}

/**
 * @brief Costruisce una matrice etichettata: dict {matrix, rows, cols}.
 * @tparam MatrixT tipo con Extent(0/1) e operator()(i, j) (Tensor<double, 2>).
 * @param row_name funtore std::size_t -> std::string_view per i nomi di riga.
 * @param col_name funtore std::size_t -> std::string_view per i nomi di colonna.
 */
template <typename MatrixT, typename RowName, typename ColName>
[[nodiscard]] inline PyObject* LabeledMatrix(const MatrixT& m, RowName&& row_name,
                                             ColName&& col_name) {
    const std::size_t nrow = static_cast<std::size_t>(m.Extent(0));
    const std::size_t ncol = static_cast<std::size_t>(m.Extent(1));

    PyRef dict(PyDict_New());
    if (!dict) {
        throw std::bad_alloc();
    }
    DictSetSteal(dict.get(), "matrix", FromDoubleMatrix(m));

    PyRef rows(PyList_New(static_cast<Py_ssize_t>(nrow)));
    if (!rows) {
        throw std::bad_alloc();
    }
    for (std::size_t i = 0; i < nrow; ++i) {
        PyObject* s = FromString(row_name(i));
        if (s == nullptr) {
            throw std::runtime_error("failed to allocate row name.");
        }
        PyList_SET_ITEM(rows.get(), static_cast<Py_ssize_t>(i), s);
    }
    DictSetSteal(dict.get(), "rows", rows.release());

    PyRef cols(PyList_New(static_cast<Py_ssize_t>(ncol)));
    if (!cols) {
        throw std::bad_alloc();
    }
    for (std::size_t j = 0; j < ncol; ++j) {
        PyObject* s = FromString(col_name(j));
        if (s == nullptr) {
            throw std::runtime_error("failed to allocate col name.");
        }
        PyList_SET_ITEM(cols.get(), static_cast<Py_ssize_t>(j), s);
    }
    DictSetSteal(dict.get(), "cols", cols.release());

    return dict.release();
}

}  // namespace PyCreate
