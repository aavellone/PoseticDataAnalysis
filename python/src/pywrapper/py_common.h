/**
 * @file py_common.h
 * @brief Fondamenta del layer Python/C++ nativo (CPython C-API puro).
 *
 * @details Analogo di rwrapper.h del package R, ma per Python:
 *  - @ref PyRef            — RAII per i riferimenti PyObject* (come RProtectGuard).
 *  - @ref TranslateException — inoltra un'eccezione C++ come eccezione Python.
 *  - capsule per POSetWrap — equivalente dell'ExternalPtr R (puntatore opaco +
 *    distruttore che libera l'oggetto C++ quando il GC di Python lo raccoglie).
 *
 * Nessuna dipendenza esterna: solo <Python.h> e la C++20 standard library.
 *
 * @author Alessandro Avellone
 */

#pragma once

#define PY_SSIZE_T_CLEAN
#include <Python.h>

// Explicit standard includes (MSVC does not pull these in transitively).
// Included here, after Python.h, so every translation unit that includes
// py_common.h gets them in the correct order.
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include "poset_wrapper.h"

// ===========================================================================
// PyRef — RAII per i riferimenti forti PyObject*
// ===========================================================================

/**
 * @class PyRef
 * @brief Possiede un riferimento forte a un PyObject* e lo rilascia (Py_XDECREF)
 * all'uscita di scope, anche sul percorso di eccezione. Non copiabile.
 */
class PyRef {
  public:
    PyRef() noexcept : obj_(nullptr) {}
    explicit PyRef(PyObject* obj) noexcept : obj_(obj) {}

    ~PyRef() { Py_XDECREF(obj_); }

    PyRef(const PyRef&) = delete;
    PyRef& operator=(const PyRef&) = delete;

    PyRef(PyRef&& other) noexcept : obj_(other.obj_) { other.obj_ = nullptr; }
    PyRef& operator=(PyRef&& other) noexcept {
        if (this != &other) {
            Py_XDECREF(obj_);
            obj_ = other.obj_;
            other.obj_ = nullptr;
        }
        return *this;
    }

    [[nodiscard]] PyObject* get() const noexcept { return obj_; }
    [[nodiscard]] explicit operator bool() const noexcept { return obj_ != nullptr; }

    /// Cede la ownership al chiamante (per il return finale verso Python).
    [[nodiscard]] PyObject* release() noexcept {
        PyObject* tmp = obj_;
        obj_ = nullptr;
        return tmp;
    }

  private:
    PyObject* obj_;
};

// ===========================================================================
// Inoltro delle eccezioni C++ verso Python
// ===========================================================================

/**
 * @brief Converte l'eccezione C++ "in volo" nella corrispondente eccezione
 * Python e imposta lo stato d'errore dell'interprete.
 * @details Da chiamare dentro un blocco catch(...) { TranslateException(); return nullptr; }.
 * std::invalid_argument -> ValueError, std::out_of_range -> IndexError,
 * ogni altra std::exception -> RuntimeError.
 */
inline void TranslateException() {
    // Se un errore Python e' gia' impostato (es. propagato da una conversione),
    // non sovrascriverlo.
    if (PyErr_Occurred()) {
        return;
    }
    try {
        throw;
    } catch (const std::invalid_argument& ex) {
        PyErr_SetString(PyExc_ValueError, ex.what());
    } catch (const std::out_of_range& ex) {
        PyErr_SetString(PyExc_IndexError, ex.what());
    } catch (const std::bad_alloc&) {
        PyErr_NoMemory();
    } catch (const std::exception& ex) {
        PyErr_SetString(PyExc_RuntimeError, ex.what());
    } catch (...) {
        PyErr_SetString(PyExc_RuntimeError, "unknown C++ exception");
    }
}

// ===========================================================================
// Capsule POSetWrap — equivalente Python dell'ExternalPtr R
// ===========================================================================

/// Nome (tag) della capsule: identifica in modo sicuro i POSetWrap.
inline constexpr const char* kPOSetCapsuleName = "poseticDataAnalysis.POSetWrap";

/**
 * @brief Distruttore registrato sulla capsule: libera il POSetWrap C++ quando
 * il garbage collector di Python raccoglie la capsule (equivalente del
 * finalizer R PointerFinalizer<T>).
 */
inline void POSetCapsuleDestructor(PyObject* capsule) {
    void* raw = PyCapsule_GetPointer(capsule, kPOSetCapsuleName);
    if (raw != nullptr) {
        delete static_cast<POSetWrap*>(raw);
    }
}

/**
 * @brief Impacchetta un unique_ptr<POSetWrap> in una capsule Python,
 * cedendo la ownership al GC di Python solo a wrapping riuscito
 * (analogo di RCreate::WrapExternalPtr).
 * @return Nuovo riferimento alla capsule, oppure nullptr con errore Python
 * impostato (in tal caso l'oggetto C++ viene liberato dallo unique_ptr).
 */
inline PyObject* MakePOSetCapsule(std::unique_ptr<POSetWrap> ptr) {
    if (!ptr) {
        Py_RETURN_NONE;
    }
    PyObject* capsule = PyCapsule_New(ptr.get(), kPOSetCapsuleName,
                                      POSetCapsuleDestructor);
    if (capsule == nullptr) {
        return nullptr;  // lo unique_ptr libera l'oggetto: nessun leak
    }
    ptr.release();  // ora la ownership e' della capsule
    return capsule;
}

/**
 * @brief Estrae e valida il POSetWrap* da una capsule ricevuta da Python.
 * @details Equivalente di RConvert::ToPOSetWrap/ToWrapPtr: protegge da
 * dereferenziazioni di puntatori non validi. Una capsule con tag errato o
 * indirizzo nullo produce un'eccezione C++ chiara invece di un crash.
 * @throws std::invalid_argument Se non e' una capsule POSetWrap valida.
 */
[[nodiscard]] inline POSetWrap* GetPOSetWrap(PyObject* capsule) {
    if (!PyCapsule_IsValid(capsule, kPOSetCapsuleName)) {
        // Ripulisce l'errore impostato da PyCapsule_IsValid, se presente.
        PyErr_Clear();
        throw std::invalid_argument(
            "expected a valid POSet handle (was the object already freed, or is "
            "it the wrong type?)");
    }
    void* raw = PyCapsule_GetPointer(capsule, kPOSetCapsuleName);
    if (raw == nullptr) {
        PyErr_Clear();
        throw std::runtime_error("POSet handle is null (already freed).");
    }
    return static_cast<POSetWrap*>(raw);
}

// ===========================================================================
// Capsule generiche per gli altri handle C++ (generatori, valutatori)
// ===========================================================================

/// Nomi (tag) delle capsule per i tipi di handle della fase valutazione.
inline constexpr const char* kLEGCapsuleName =
    "poseticDataAnalysis.LinearGeneratorWrap";
inline constexpr const char* kBDMRPCapsuleName =
    "poseticDataAnalysis.BubleyDyerMRPGenerator";
inline constexpr const char* kBDEvalCapsuleName =
    "poseticDataAnalysis.BubleyDyerEvaluationGenerator";

/**
 * @brief Distruttore generico per capsule che possiedono un oggetto C++ di tipo T.
 * @details Recupera il nome della capsule per estrarre il puntatore in modo
 * type-safe, poi lo distrugge (equivalente del PointerFinalizer<T> di R).
 */
template <typename T>
inline void TypedCapsuleDestructor(PyObject* capsule) {
    void* raw = PyCapsule_GetPointer(capsule, PyCapsule_GetName(capsule));
    if (raw != nullptr) {
        delete static_cast<T*>(raw);
    }
}

/**
 * @brief Impacchetta uno unique_ptr<T> in una capsule con il tag @p name.
 * @return Nuovo riferimento alla capsule, o nullptr (con errore Python) — in tal
 * caso lo unique_ptr libera l'oggetto.
 */
template <typename T>
[[nodiscard]] inline PyObject* MakeCapsule(std::unique_ptr<T> ptr, const char* name) {
    if (!ptr) {
        Py_RETURN_NONE;
    }
    PyObject* capsule = PyCapsule_New(ptr.get(), name, TypedCapsuleDestructor<T>);
    if (capsule == nullptr) {
        return nullptr;
    }
    ptr.release();
    return capsule;
}

/**
 * @brief Estrae e valida un T* da una capsule con tag @p name.
 * @throws std::invalid_argument Se non e' una capsule valida del tipo atteso.
 */
template <typename T>
[[nodiscard]] inline T* GetCapsule(PyObject* capsule, const char* name,
                                   const char* what) {
    if (!PyCapsule_IsValid(capsule, name)) {
        PyErr_Clear();
        throw std::invalid_argument(
            std::string("expected a valid ") + what +
            " handle (was it already freed, or is it the wrong type?)");
    }
    void* raw = PyCapsule_GetPointer(capsule, name);
    if (raw == nullptr) {
        PyErr_Clear();
        throw std::runtime_error(std::string(what) + " handle is null (already freed).");
    }
    return static_cast<T*>(raw);
}
