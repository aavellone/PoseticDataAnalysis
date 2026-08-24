/**
 * @file pywrapper_relations.cpp
 * @brief Entry-point Python per le proprietà delle relazioni binarie.
 *
 * @details Speculare a rwrapper_relations.cpp: opera su una relazione data come
 * (elementi, coppie) oppure solo coppie, costruisce la matrice di adiacenza e
 * interroga le funzioni generic::. Le chiusure restituiscono list[tuple].
 */

#include "pywrapper.h"
#include "py_common.h"
#include "py_convert.h"

#include "generic_functions.h"
#include "tensor.h"

#include <cstdint>
#include <string>
#include <vector>

namespace {

/// Restituisce gli archi presenti (adj(i,j)==1) come list[tuple(str, str)].
PyObject* EdgesToPy(const Tensor<std::uint8_t, 2>& adj,
                    const std::vector<std::string>& elements) {
    const std::uint64_t n = static_cast<std::uint64_t>(elements.size());
    PyRef list(PyList_New(0));
    if (!list) {
        throw std::bad_alloc();
    }
    for (std::uint64_t i = 0; i < n; ++i) {
        for (std::uint64_t j = 0; j < n; ++j) {
            if (adj(i, j) == 1) {
                PyObject* a = PyCreate::FromString(elements[i]);
                PyObject* b = PyCreate::FromString(elements[j]);
                PyObject* tup = (a && b) ? PyTuple_New(2) : nullptr;
                if (tup == nullptr) {
                    Py_XDECREF(a);
                    Py_XDECREF(b);
                    throw std::runtime_error("failed to build edge tuple.");
                }
                PyTuple_SET_ITEM(tup, 0, a);
                PyTuple_SET_ITEM(tup, 1, b);
                if (PyList_Append(list.get(), tup) != 0) {
                    Py_DECREF(tup);
                    throw std::runtime_error("failed to append edge.");
                }
                Py_DECREF(tup);  // PyList_Append incrementa il refcount
            }
        }
    }
    return list.release();
}

}  // namespace

// --- (set, rel) --------------------------------------------------------------

PyObject* pyx_is_reflexive(PyObject* /*self*/, PyObject* args) {
    PyObject* set_obj = nullptr;
    PyObject* rel_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OO", &set_obj, &rel_obj)) {
        return nullptr;
    }
    try {
        std::vector<std::string> elements;
        auto adj = PyConvert::ToAdjacencyMatrix(set_obj, rel_obj, elements);
        return PyBool_FromLong(generic::IsReflexive(elements.size(), adj) ? 1 : 0);
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

PyObject* pyx_is_preorder(PyObject* /*self*/, PyObject* args) {
    PyObject* set_obj = nullptr;
    PyObject* rel_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OO", &set_obj, &rel_obj)) {
        return nullptr;
    }
    try {
        std::vector<std::string> elements;
        auto adj = PyConvert::ToAdjacencyMatrix(set_obj, rel_obj, elements);
        return PyBool_FromLong(generic::IsPreorder(elements.size(), adj) ? 1 : 0);
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

PyObject* pyx_is_partial_order(PyObject* /*self*/, PyObject* args) {
    PyObject* set_obj = nullptr;
    PyObject* rel_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OO", &set_obj, &rel_obj)) {
        return nullptr;
    }
    try {
        std::vector<std::string> elements;
        auto adj = PyConvert::ToAdjacencyMatrix(set_obj, rel_obj, elements);
        return PyBool_FromLong(generic::IsPartialOrder(elements.size(), adj) ? 1 : 0);
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

PyObject* pyx_reflexive_closure(PyObject* /*self*/, PyObject* args) {
    PyObject* set_obj = nullptr;
    PyObject* rel_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OO", &set_obj, &rel_obj)) {
        return nullptr;
    }
    try {
        std::vector<std::string> elements;
        auto adj = PyConvert::ToAdjacencyMatrix(set_obj, rel_obj, elements);
        generic::ReflexiveClosureInPlace(elements.size(), adj);
        return EdgesToPy(adj, elements);
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

// --- (rel) only --------------------------------------------------------------

PyObject* pyx_is_symmetric(PyObject* /*self*/, PyObject* args) {
    PyObject* rel_obj = nullptr;
    if (!PyArg_ParseTuple(args, "O", &rel_obj)) {
        return nullptr;
    }
    try {
        std::vector<std::string> elements;
        auto adj = PyConvert::ToAdjacencyMatrixFromEdges(rel_obj, elements);
        return PyBool_FromLong(generic::IsSymmetric(elements.size(), adj) ? 1 : 0);
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

PyObject* pyx_is_antisymmetric(PyObject* /*self*/, PyObject* args) {
    PyObject* rel_obj = nullptr;
    if (!PyArg_ParseTuple(args, "O", &rel_obj)) {
        return nullptr;
    }
    try {
        std::vector<std::string> elements;
        auto adj = PyConvert::ToAdjacencyMatrixFromEdges(rel_obj, elements);
        return PyBool_FromLong(generic::IsAntisymmetric(elements.size(), adj) ? 1 : 0);
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

PyObject* pyx_is_transitive(PyObject* /*self*/, PyObject* args) {
    PyObject* rel_obj = nullptr;
    if (!PyArg_ParseTuple(args, "O", &rel_obj)) {
        return nullptr;
    }
    try {
        std::vector<std::string> elements;
        auto adj = PyConvert::ToAdjacencyMatrixFromEdges(rel_obj, elements);
        return PyBool_FromLong(generic::IsTransitive(elements.size(), adj) ? 1 : 0);
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

PyObject* pyx_transitive_closure(PyObject* /*self*/, PyObject* args) {
    PyObject* rel_obj = nullptr;
    if (!PyArg_ParseTuple(args, "O", &rel_obj)) {
        return nullptr;
    }
    try {
        std::vector<std::string> elements;
        auto adj = PyConvert::ToAdjacencyMatrixFromEdges(rel_obj, elements);
        generic::TransitiveClosureInPlace(elements.size(), adj);
        return EdgesToPy(adj, elements);
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}
