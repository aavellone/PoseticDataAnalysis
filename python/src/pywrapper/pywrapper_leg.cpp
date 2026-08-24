/**
 * @file pywrapper_leg.cpp
 * @brief Entry-point Python per i generatori di estensioni lineari (LEG).
 *
 * @details Speculare a rwrapper_leg.cpp:
 *  - generatore esatto (Tree-of-Ideals): build_le_generator / leg_get;
 *  - campionatore MCMC Bubley-Dyer: build_bubley_dyer_le_generator / leg_bubley_dyer_get.
 * I seed a 64 bit sono restituiti come stringa (reimmetterli riproduce la stessa
 * sequenza).
 */

#include "pywrapper.h"
#include "py_common.h"
#include "py_convert.h"
#include "py_linear_generator_wrapper.h"

#include "poset_wrapper.h"
#include "random.h"

#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace {

/// Converte un vettore di estensioni (ognuna vettore di nomi) in list[list[str]].
PyObject* ExtensionsToPy(const std::vector<std::vector<std::string>>& les) {
    PyRef out(PyList_New(static_cast<Py_ssize_t>(les.size())));
    if (!out) {
        throw std::bad_alloc();
    }
    for (std::size_t c = 0; c < les.size(); ++c) {
        const auto& le = les[c];
        PyObject* col = PyList_New(static_cast<Py_ssize_t>(le.size()));
        if (col == nullptr) {
            throw std::bad_alloc();
        }
        for (std::size_t r = 0; r < le.size(); ++r) {
            PyObject* s = PyCreate::FromString(le[r]);
            if (s == nullptr) {
                Py_DECREF(col);
                throw std::runtime_error("failed to allocate Python string.");
            }
            PyList_SET_ITEM(col, static_cast<Py_ssize_t>(r), s);
        }
        PyList_SET_ITEM(out.get(), static_cast<Py_ssize_t>(c), col);
    }
    return out.release();
}

}  // namespace

// build_le_generator(poset) -> LEG capsule (exact Tree-of-Ideals)
PyObject* pyx_build_le_generator(PyObject* /*self*/, PyObject* args) {
    PyObject* poset_obj = nullptr;
    if (!PyArg_ParseTuple(args, "O", &poset_obj)) {
        return nullptr;
    }
    try {
        const POSetWrap* poset_wrap = GetPOSetWrap(poset_obj);
        return MakeCapsule(LinearGeneratorWrap::BuildLEGenerator(poset_wrap), kLEGCapsuleName);
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

// build_bubley_dyer_le_generator(poset, seed) -> (LEG capsule, seed_str)
PyObject* pyx_build_bubley_dyer_le_generator(PyObject* /*self*/, PyObject* args) {
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
            LinearGeneratorWrap::BuildBubleyDyerGenerator(poset_wrap, seed.value()),
            kLEGCapsuleName);
        if (capsule == nullptr) {
            return nullptr;
        }
        PyObject* seed_str = PyCreate::FromSeed(seed.value());
        if (seed_str == nullptr) {
            Py_DECREF(capsule);
            return nullptr;
        }
        PyObject* tup = PyTuple_New(2);
        if (tup == nullptr) {
            Py_DECREF(capsule);
            Py_DECREF(seed_str);
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

// leg_get(gen, from_start, quante, output_every) -> list[list[str]]  (exact)
PyObject* pyx_leg_get(PyObject* /*self*/, PyObject* args) {
    PyObject* gen_obj = nullptr;
    PyObject* from_start_obj = nullptr;
    PyObject* quante_obj = nullptr;
    PyObject* output_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OOOO", &gen_obj, &from_start_obj, &quante_obj, &output_obj)) {
        return nullptr;
    }
    try {
        auto* gen = GetCapsule<LinearGeneratorWrap>(gen_obj, kLEGCapsuleName,
                                                    "linear extension generator");
        const bool from_start = PyConvert::ToBool(from_start_obj);
        const auto quante = PyConvert::ToOptionalUInt(quante_obj);
        const auto output_every = PyConvert::ToOptionalUInt(output_obj);
        const auto les = gen->GetFromLE(from_start, quante, output_every);
        return ExtensionsToPy(les);
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}

// leg_bubley_dyer_get(gen, from_start, quante, errore, output_every) -> list[list[str]]
PyObject* pyx_leg_bubley_dyer_get(PyObject* /*self*/, PyObject* args) {
    PyObject* gen_obj = nullptr;
    PyObject* from_start_obj = nullptr;
    PyObject* quante_obj = nullptr;
    PyObject* errore_obj = nullptr;
    PyObject* output_obj = nullptr;
    if (!PyArg_ParseTuple(args, "OOOOO", &gen_obj, &from_start_obj, &quante_obj,
                          &errore_obj, &output_obj)) {
        return nullptr;
    }
    try {
        auto* gen = GetCapsule<LinearGeneratorWrap>(gen_obj, kLEGCapsuleName,
                                                    "linear extension generator");
        const bool from_start = PyConvert::ToBool(from_start_obj);
        const auto quante = PyConvert::ToOptionalUInt(quante_obj);
        const auto errore = PyConvert::ToOptionalDouble(errore_obj);
        const auto output_every = PyConvert::ToOptionalUInt(output_obj);
        const auto les = gen->GetFromBubleyDyer(from_start, quante, errore, output_every);
        return ExtensionsToPy(les);
    } catch (...) {
        TranslateException();
        return nullptr;
    }
}
