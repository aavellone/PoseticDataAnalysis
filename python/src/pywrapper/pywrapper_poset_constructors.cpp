/**
 * @file pywrapper_poset_constructors.cpp
 * @brief Entry-point Python per i costruttori di POSet.
 *
 * @details Speculare a rwrapper_poset_constructors.cpp: stesso pattern in tre fasi
 *  1. Estrazione (PyConvert) degli argomenti Python in strutture C++20.
 *  2. Elaborazione tramite i Factory Method di POSetWrap.
 *  3. Restituzione: la ownership dello std::unique_ptr passa al GC di Python
 *     tramite una capsule (MakePOSetCapsule), come l'ExternalPtr in R.
 */

#include "pywrapper.h"
#include "py_common.h"
#include "py_convert.h"

#include "poset_wrapper.h"

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

// build_poset(elements, relations) -> capsule
PyObject* pyx_build_poset(PyObject* /*self*/, PyObject* args) {
    PyObject* elements_obj = nullptr;
    PyObject* relations_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OO", &elements_obj, &relations_obj)) {
        return nullptr;
    }
    try {
        auto elements = PyConvert::ToStringVector(elements_obj);
        auto relations = PyConvert::ToPairVector(relations_obj);
        return MakePOSetCapsule(POSetWrap::BuildPoset(elements, relations));
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

// build_linear_poset(elements) -> capsule
PyObject* pyx_build_linear_poset(PyObject* /*self*/, PyObject* args) {
    PyObject* elements_obj = nullptr;
    if (!PyArg_ParseTuple(args, "O", &elements_obj)) {
        return nullptr;
    }
    try {
        auto elements = PyConvert::ToStringVector(elements_obj);
        return MakePOSetCapsule(POSetWrap::BuildLinear(elements));
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

// Helper comune per i costruttori che ricevono una sequenza di POSet.
namespace {
template <typename Factory>
PyObject* BuildFromPosetList(PyObject* args, Factory&& factory) {
    PyObject* posets_obj = nullptr;
    if (!PyArg_ParseTuple(args, "O", &posets_obj)) {
        return nullptr;
    }
    try {
        auto posets = PyConvert::ToPOSetWrapVector(posets_obj);
        return MakePOSetCapsule(factory(posets));
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}
}  // namespace

PyObject* pyx_build_product_poset(PyObject* /*self*/, PyObject* args) {
    return BuildFromPosetList(args, POSetWrap::BuildProduct);
}

PyObject* pyx_build_lexicographic_product_poset(PyObject* /*self*/, PyObject* args) {
    return BuildFromPosetList(args, POSetWrap::BuildLexicographicProduct);
}

PyObject* pyx_build_intersection_poset(PyObject* /*self*/, PyObject* args) {
    return BuildFromPosetList(args, POSetWrap::BuildIntersection);
}

PyObject* pyx_build_linear_sum_poset(PyObject* /*self*/, PyObject* args) {
    return BuildFromPosetList(args, POSetWrap::BuildLinearSum);
}

PyObject* pyx_build_disjoint_sum_poset(PyObject* /*self*/, PyObject* args) {
    return BuildFromPosetList(args, POSetWrap::BuildDisjointSum);
}

// build_lifting_poset(poset, new_element) -> capsule
PyObject* pyx_build_lifting_poset(PyObject* /*self*/, PyObject* args) {
    PyObject* poset_obj = nullptr;
    PyObject* new_element_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OO", &poset_obj, &new_element_obj)) {
        return nullptr;
    }
    try {
        const POSetWrap* base = GetPOSetWrap(poset_obj);
        const std::string new_element = PyConvert::ToString(new_element_obj, "new_element");
        return MakePOSetCapsule(POSetWrap::BuildLiftingPOSet(base, new_element));
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

// build_binary_variable_poset(variables) -> capsule
PyObject* pyx_build_binary_variable_poset(PyObject* /*self*/, PyObject* args) {
    PyObject* variables_obj = nullptr;
    if (!PyArg_ParseTuple(args, "O", &variables_obj)) {
        return nullptr;
    }
    try {
        auto variables = PyConvert::ToStringVector(variables_obj, "variables");
        return MakePOSetCapsule(POSetWrap::BuildBinaryVariablePOSet(variables));
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

// build_fence_poset(elements, orientation) -> capsule
PyObject* pyx_build_fence_poset(PyObject* /*self*/, PyObject* args) {
    PyObject* elements_obj = nullptr;
    PyObject* orientation_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OO", &elements_obj, &orientation_obj)) {
        return nullptr;
    }
    try {
        auto elements = PyConvert::ToStringVector(elements_obj);
        const bool orientation = PyConvert::ToBool(orientation_obj);
        return MakePOSetCapsule(POSetWrap::BuildFencePOSet(elements, orientation));
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

// build_crown_poset(elements_1, elements_2) -> capsule
PyObject* pyx_build_crown_poset(PyObject* /*self*/, PyObject* args) {
    PyObject* e1_obj = nullptr;
    PyObject* e2_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OO", &e1_obj, &e2_obj)) {
        return nullptr;
    }
    try {
        auto elements_1 = PyConvert::ToStringVector(e1_obj);
        auto elements_2 = PyConvert::ToStringVector(e2_obj);
        const std::size_t n = elements_1.size();
        if (elements_2.size() != n) {
            throw std::invalid_argument(
                "crown POSet: the two element lists must have the same length.");
        }

        std::vector<std::string> all_elements;
        all_elements.reserve(2 * n);
        for (std::size_t k = 0; k < n; ++k) {
            all_elements.push_back(elements_1[k]);
            all_elements.push_back(elements_2[k]);
        }

        std::vector<std::pair<std::string, std::string>> comparabilities;
        if (n > 0) {
            comparabilities.reserve(n * (n - 1));
        }
        for (std::size_t k = 0; k < n; ++k) {
            for (std::size_t h = 0; h < n; ++h) {
                if (k != h) {
                    comparabilities.emplace_back(elements_1[k], elements_2[h]);
                }
            }
        }
        return MakePOSetCapsule(POSetWrap::BuildPoset(all_elements, comparabilities));
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

// build_dual_poset(poset) -> capsule
PyObject* pyx_build_dual_poset(PyObject* /*self*/, PyObject* args) {
    PyObject* poset_obj = nullptr;
    if (!PyArg_ParseTuple(args, "O", &poset_obj)) {
        return nullptr;
    }
    try {
        const POSetWrap* base = GetPOSetWrap(poset_obj);
        return MakePOSetCapsule(POSetWrap::BuildDualPOSet(base));
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}
