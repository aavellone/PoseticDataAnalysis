/**
 * @file pywrapper_poset_properties.cpp
 * @brief Entry-point Python per le query di proprietà e relazioni del POSet.
 *
 * @details Speculare a rwrapper_poset_properties.cpp. Ogni funzione:
 *  1. estrae il POSetWrap dalla capsule (GetPOSetWrap, con validazione),
 *  2. interroga il core C++,
 *  3. converte il risultato in un oggetto Python nativo (PyCreate).
 *
 * Le funzioni pairwise (dominates, is_dominated_by, ...) accettano DUE sequenze
 * di nomi e restituiscono una list di bool elemento per elemento, come le
 * corrispondenti funzioni vettoriali di R.
 */

#include "pywrapper.h"
#include "py_common.h"
#include "py_convert.h"

#include "poset_wrapper.h"
#include "poset.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace {

/// Estrae il POSet* dalla capsule passata come primo (e unico) argomento.
POSet* PosetFromArgs(PyObject* args) {
    PyObject* poset_obj = nullptr;
    if (!PyArg_ParseTuple(args, "O", &poset_obj)) {
        return nullptr;  // errore Python gia' impostato
    }
    return GetPOSetWrap(poset_obj)->GetPOSet();
}

/// Converte una sequenza Python di nomi in vettore di id del poset.
std::vector<std::uint64_t> NamesToIds(POSet* poset, PyObject* names_obj) {
    const auto names = PyConvert::ToStringVector(names_obj, "elements");
    std::vector<std::uint64_t> ids;
    ids.reserve(names.size());
    for (const auto& name : names) {
        ids.push_back(poset->GetElementId(name));
    }
    return ids;
}

/// Applica un predicato pairwise a due sequenze di nomi -> list di bool.
template <typename Pred>
PyObject* Pairwise(PyObject* args, Pred&& pred) {
    PyObject* poset_obj = nullptr;
    PyObject* a_obj = nullptr;
    PyObject* b_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OOO", &poset_obj, &a_obj, &b_obj)) {
        return nullptr;
    }
    try {
        POSet* poset = GetPOSetWrap(poset_obj)->GetPOSet();
        const auto a = PyConvert::ToStringVector(a_obj, "elements");
        const auto b = PyConvert::ToStringVector(b_obj, "elements");
        if (a.size() != b.size()) {
            throw std::invalid_argument("the two element sequences must have the same length.");
        }
        PyRef list(PyList_New(static_cast<Py_ssize_t>(a.size())));
        if (!list) {
            throw std::bad_alloc();
        }
        for (std::size_t k = 0; k < a.size(); ++k) {
            const std::uint64_t n1 = poset->GetElementId(a[k]);
            const std::uint64_t n2 = poset->GetElementId(b[k]);
            PyObject* v = PyBool_FromLong(pred(poset, n1, n2) ? 1 : 0);
            PyList_SET_ITEM(list.get(), static_cast<Py_ssize_t>(k), v);
        }
        return list.release();
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

}  // namespace

// --- Struttura di base ------------------------------------------------------

PyObject* pyx_size(PyObject* /*self*/, PyObject* args) {
    try {
        POSet* poset = PosetFromArgs(args);
        if (poset == nullptr) {
            return nullptr;
        }
        return PyLong_FromUnsignedLongLong(static_cast<unsigned long long>(poset->size()));
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

PyObject* pyx_elements(PyObject* /*self*/, PyObject* args) {
    try {
        POSet* poset = PosetFromArgs(args);
        if (poset == nullptr) {
            return nullptr;
        }
        const std::size_t n = static_cast<std::size_t>(poset->size());
        PyRef list(PyList_New(static_cast<Py_ssize_t>(n)));
        if (!list) {
            throw std::bad_alloc();
        }
        for (std::size_t i = 0; i < n; ++i) {
            PyObject* s = PyCreate::FromString(poset->GetElementName(i));
            PyList_SET_ITEM(list.get(), static_cast<Py_ssize_t>(i), s);
        }
        return list.release();
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

PyObject* pyx_incidence_matrix(PyObject* /*self*/, PyObject* args) {
    try {
        POSet* poset = PosetFromArgs(args);
        if (poset == nullptr) {
            return nullptr;
        }
        return PyCreate::FromIntMatrix(poset->IncidenceMatrix());
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

PyObject* pyx_cover_matrix(PyObject* /*self*/, PyObject* args) {
    try {
        POSet* poset = PosetFromArgs(args);
        if (poset == nullptr) {
            return nullptr;
        }
        return PyCreate::FromIntMatrix(poset->CoverMatrix());
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

PyObject* pyx_order_relation(PyObject* /*self*/, PyObject* args) {
    try {
        POSet* poset = PosetFromArgs(args);
        if (poset == nullptr) {
            return nullptr;
        }
        auto relations = poset->OrderRelation();  // list<pair<string,string>>
        return PyCreate::FromPairs(
            relations, static_cast<Py_ssize_t>(relations.size()),
            [](const std::string& s) -> std::string_view { return s; });
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

PyObject* pyx_cover_relation(PyObject* /*self*/, PyObject* args) {
    try {
        POSet* poset = PosetFromArgs(args);
        if (poset == nullptr) {
            return nullptr;
        }
        auto relations = poset->CoverRelation();  // list<pair<uint64,uint64>>
        return PyCreate::FromPairs(
            relations, static_cast<Py_ssize_t>(relations.size()),
            [poset](std::uint64_t id) -> std::string_view {
                return poset->GetElementName(static_cast<std::size_t>(id));
            });
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

PyObject* pyx_incomparabilities(PyObject* /*self*/, PyObject* args) {
    try {
        POSet* poset = PosetFromArgs(args);
        if (poset == nullptr) {
            return nullptr;
        }
        auto pairs = poset->Incomparabilities();  // list<pair<uint64,uint64>>
        return PyCreate::FromPairs(
            pairs, static_cast<Py_ssize_t>(pairs.size()),
            [poset](std::uint64_t id) -> std::string_view {
                return poset->GetElementName(static_cast<std::size_t>(id));
            });
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

// --- Confronti pairwise -----------------------------------------------------

PyObject* pyx_is_dominated_by(PyObject* /*self*/, PyObject* args) {
    // v1 dominato da v2: v1 == v2 oppure v2 > v1
    return Pairwise(args, [](POSet* p, std::uint64_t n1, std::uint64_t n2) {
        return n1 == n2 || p->GreaterThan(n2, n1);
    });
}

PyObject* pyx_dominates(PyObject* /*self*/, PyObject* args) {
    // v1 domina v2: v1 == v2 oppure v1 > v2
    return Pairwise(args, [](POSet* p, std::uint64_t n1, std::uint64_t n2) {
        return n1 == n2 || p->GreaterThan(n1, n2);
    });
}

PyObject* pyx_is_comparable_with(PyObject* /*self*/, PyObject* args) {
    // comparabili: uguali oppure uno domina l'altro (logica corretta;
    // il wrapper R aveva una doppia assegnazione che sovrascriveva questo caso).
    return Pairwise(args, [](POSet* p, std::uint64_t n1, std::uint64_t n2) {
        return n1 == n2 || p->GreaterThan(n1, n2) || p->GreaterThan(n2, n1);
    });
}

PyObject* pyx_is_incomparable_with(PyObject* /*self*/, PyObject* args) {
    // incomparabili: diversi e nessuna relazione d'ordine nei due versi
    return Pairwise(args, [](POSet* p, std::uint64_t n1, std::uint64_t n2) {
        return n1 != n2 && !p->GreaterThan(n1, n2) && !p->GreaterThan(n2, n1);
    });
}

// --- Upset / Downset --------------------------------------------------------

PyObject* pyx_upset_of(PyObject* /*self*/, PyObject* args) {
    PyObject* poset_obj = nullptr;
    PyObject* names_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OO", &poset_obj, &names_obj)) {
        return nullptr;
    }
    try {
        POSet* poset = GetPOSetWrap(poset_obj)->GetPOSet();
        auto result = poset->UpSet(NamesToIds(poset, names_obj));
        return PyCreate::NamesFromIds(result, *poset,
                                      static_cast<Py_ssize_t>(result.count()));
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

PyObject* pyx_downset_of(PyObject* /*self*/, PyObject* args) {
    PyObject* poset_obj = nullptr;
    PyObject* names_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OO", &poset_obj, &names_obj)) {
        return nullptr;
    }
    try {
        POSet* poset = GetPOSetWrap(poset_obj)->GetPOSet();
        auto result = poset->DownSet(NamesToIds(poset, names_obj));
        return PyCreate::NamesFromIds(result, *poset,
                                      static_cast<Py_ssize_t>(result.count()));
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

PyObject* pyx_is_upset(PyObject* /*self*/, PyObject* args) {
    PyObject* poset_obj = nullptr;
    PyObject* names_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OO", &poset_obj, &names_obj)) {
        return nullptr;
    }
    try {
        POSet* poset = GetPOSetWrap(poset_obj)->GetPOSet();
        return PyBool_FromLong(poset->IsUpSet(NamesToIds(poset, names_obj)) ? 1 : 0);
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

PyObject* pyx_is_downset(PyObject* /*self*/, PyObject* args) {
    PyObject* poset_obj = nullptr;
    PyObject* names_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OO", &poset_obj, &names_obj)) {
        return nullptr;
    }
    try {
        POSet* poset = GetPOSetWrap(poset_obj)->GetPOSet();
        return PyBool_FromLong(poset->IsDownSet(NamesToIds(poset, names_obj)) ? 1 : 0);
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

// --- Set di (in)comparabilità di un elemento --------------------------------

PyObject* pyx_comparability_set_of(PyObject* /*self*/, PyObject* args) {
    PyObject* poset_obj = nullptr;
    PyObject* element_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OO", &poset_obj, &element_obj)) {
        return nullptr;
    }
    try {
        POSet* poset = GetPOSetWrap(poset_obj)->GetPOSet();
        const std::string name = PyConvert::ToString(element_obj, "element");
        auto result = poset->ComparabilitySetOf(poset->GetElementId(name));
        return PyCreate::NamesFromIds(result, *poset,
                                      static_cast<Py_ssize_t>(result.count()));
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

PyObject* pyx_incomparability_set_of(PyObject* /*self*/, PyObject* args) {
    PyObject* poset_obj = nullptr;
    PyObject* element_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OO", &poset_obj, &element_obj)) {
        return nullptr;
    }
    try {
        POSet* poset = GetPOSetWrap(poset_obj)->GetPOSet();
        const std::string name = PyConvert::ToString(element_obj, "element");
        auto result = poset->IncomparabilitySetOf(poset->GetElementId(name));
        return PyCreate::NamesFromIds(result, *poset,
                                      static_cast<Py_ssize_t>(result.count()));
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

// --- Elementi estremali -----------------------------------------------------

PyObject* pyx_maximals(PyObject* /*self*/, PyObject* args) {
    try {
        POSet* poset = PosetFromArgs(args);
        if (poset == nullptr) {
            return nullptr;
        }
        auto result = poset->Maximals();
        return PyCreate::NamesFromIds(result, *poset,
                                      static_cast<Py_ssize_t>(result.count()));
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

PyObject* pyx_minimals(PyObject* /*self*/, PyObject* args) {
    try {
        POSet* poset = PosetFromArgs(args);
        if (poset == nullptr) {
            return nullptr;
        }
        auto result = poset->Minimals();
        return PyCreate::NamesFromIds(result, *poset,
                                      static_cast<Py_ssize_t>(result.count()));
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

PyObject* pyx_is_maximal(PyObject* /*self*/, PyObject* args) {
    PyObject* poset_obj = nullptr;
    PyObject* element_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OO", &poset_obj, &element_obj)) {
        return nullptr;
    }
    try {
        POSet* poset = GetPOSetWrap(poset_obj)->GetPOSet();
        const std::string name = PyConvert::ToString(element_obj, "element");
        return PyBool_FromLong(poset->IsMaximal(poset->GetElementId(name)) ? 1 : 0);
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

PyObject* pyx_is_minimal(PyObject* /*self*/, PyObject* args) {
    PyObject* poset_obj = nullptr;
    PyObject* element_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OO", &poset_obj, &element_obj)) {
        return nullptr;
    }
    try {
        POSet* poset = GetPOSetWrap(poset_obj)->GetPOSet();
        const std::string name = PyConvert::ToString(element_obj, "element");
        return PyBool_FromLong(poset->IsMinimal(poset->GetElementId(name)) ? 1 : 0);
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

// --- Meet / Join ------------------------------------------------------------

PyObject* pyx_meet(PyObject* /*self*/, PyObject* args) {
    PyObject* poset_obj = nullptr;
    PyObject* names_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OO", &poset_obj, &names_obj)) {
        return nullptr;
    }
    try {
        POSet* poset = GetPOSetWrap(poset_obj)->GetPOSet();
        auto result = poset->Meet(NamesToIds(poset, names_obj));
        if (!result.has_value()) {
            Py_RETURN_NONE;
        }
        return PyCreate::FromString(poset->GetElementName(static_cast<std::size_t>(result.value())));
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

PyObject* pyx_join(PyObject* /*self*/, PyObject* args) {
    PyObject* poset_obj = nullptr;
    PyObject* names_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OO", &poset_obj, &names_obj)) {
        return nullptr;
    }
    try {
        POSet* poset = GetPOSetWrap(poset_obj)->GetPOSet();
        auto result = poset->Join(NamesToIds(poset, names_obj));
        if (!result.has_value()) {
            Py_RETURN_NONE;
        }
        return PyCreate::FromString(poset->GetElementName(static_cast<std::size_t>(result.value())));
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

// --- Estensione -------------------------------------------------------------

PyObject* pyx_is_extension_of(PyObject* /*self*/, PyObject* args) {
    PyObject* poset_obj_1 = nullptr;
    PyObject* poset_obj_2 = nullptr;
    if (!PyArg_ParseTuple(args, "OO", &poset_obj_1, &poset_obj_2)) {
        return nullptr;
    }
    try {
        POSet* poset_1 = GetPOSetWrap(poset_obj_1)->GetPOSet();
        POSet* poset_2 = GetPOSetWrap(poset_obj_2)->GetPOSet();
        return PyBool_FromLong(poset_1->IsExtensionOf(*poset_2) ? 1 : 0);
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}
