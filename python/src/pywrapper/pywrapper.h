/**
 * @file pywrapper.h
 * @brief Prototipi delle PyCFunction registrate nel modulo _core.
 *
 * Ogni funzione e' un entry-point Python (METH_VARARGS) speculare a una entry
 * `extern "C"` del package R. Sono raggruppate per modulo sorgente:
 *  - pywrapper_poset_constructors.cpp
 *  - pywrapper_poset_properties.cpp
 */

#pragma once

#define PY_SSIZE_T_CLEAN
#include <Python.h>

// --- Costruttori (pywrapper_poset_constructors.cpp) ------------------------
PyObject* pyx_build_poset(PyObject* self, PyObject* args);
PyObject* pyx_build_linear_poset(PyObject* self, PyObject* args);
PyObject* pyx_build_product_poset(PyObject* self, PyObject* args);
PyObject* pyx_build_lexicographic_product_poset(PyObject* self, PyObject* args);
PyObject* pyx_build_intersection_poset(PyObject* self, PyObject* args);
PyObject* pyx_build_linear_sum_poset(PyObject* self, PyObject* args);
PyObject* pyx_build_disjoint_sum_poset(PyObject* self, PyObject* args);
PyObject* pyx_build_lifting_poset(PyObject* self, PyObject* args);
PyObject* pyx_build_binary_variable_poset(PyObject* self, PyObject* args);
PyObject* pyx_build_fence_poset(PyObject* self, PyObject* args);
PyObject* pyx_build_crown_poset(PyObject* self, PyObject* args);
PyObject* pyx_build_dual_poset(PyObject* self, PyObject* args);

// --- Proprietà e relazioni (pywrapper_poset_properties.cpp) -----------------
PyObject* pyx_elements(PyObject* self, PyObject* args);
PyObject* pyx_incidence_matrix(PyObject* self, PyObject* args);
PyObject* pyx_order_relation(PyObject* self, PyObject* args);
PyObject* pyx_cover_relation(PyObject* self, PyObject* args);
PyObject* pyx_cover_matrix(PyObject* self, PyObject* args);
PyObject* pyx_is_dominated_by(PyObject* self, PyObject* args);
PyObject* pyx_dominates(PyObject* self, PyObject* args);
PyObject* pyx_is_comparable_with(PyObject* self, PyObject* args);
PyObject* pyx_is_incomparable_with(PyObject* self, PyObject* args);
PyObject* pyx_upset_of(PyObject* self, PyObject* args);
PyObject* pyx_is_upset(PyObject* self, PyObject* args);
PyObject* pyx_downset_of(PyObject* self, PyObject* args);
PyObject* pyx_is_downset(PyObject* self, PyObject* args);
PyObject* pyx_comparability_set_of(PyObject* self, PyObject* args);
PyObject* pyx_incomparability_set_of(PyObject* self, PyObject* args);
PyObject* pyx_maximals(PyObject* self, PyObject* args);
PyObject* pyx_minimals(PyObject* self, PyObject* args);
PyObject* pyx_is_maximal(PyObject* self, PyObject* args);
PyObject* pyx_is_minimal(PyObject* self, PyObject* args);
PyObject* pyx_meet(PyObject* self, PyObject* args);
PyObject* pyx_join(PyObject* self, PyObject* args);
PyObject* pyx_incomparabilities(PyObject* self, PyObject* args);
PyObject* pyx_is_extension_of(PyObject* self, PyObject* args);
PyObject* pyx_size(PyObject* self, PyObject* args);

// --- Generatori di estensioni lineari (pywrapper_leg.cpp) -------------------
PyObject* pyx_build_le_generator(PyObject* self, PyObject* args);
PyObject* pyx_build_bubley_dyer_le_generator(PyObject* self, PyObject* args);
PyObject* pyx_leg_get(PyObject* self, PyObject* args);
PyObject* pyx_leg_bubley_dyer_get(PyObject* self, PyObject* args);

// --- Valutazione (pywrapper_poset_evaluation.cpp) ---------------------------
PyObject* pyx_exact_mrp(PyObject* self, PyObject* args);
PyObject* pyx_build_bubley_dyer_mrp_generator(PyObject* self, PyObject* args);
PyObject* pyx_bubley_dyer_mrp(PyObject* self, PyObject* args);
PyObject* pyx_build_bubley_dyer_evaluation_generator(PyObject* self, PyObject* args);
PyObject* pyx_bubley_dyer_evaluation(PyObject* self, PyObject* args);
PyObject* pyx_exact_evaluation(PyObject* self, PyObject* args);
PyObject* pyx_bls_dominance(PyObject* self, PyObject* args);

// --- Proprietà delle relazioni (pywrapper_relations.cpp) --------------------
PyObject* pyx_is_reflexive(PyObject* self, PyObject* args);
PyObject* pyx_is_symmetric(PyObject* self, PyObject* args);
PyObject* pyx_is_antisymmetric(PyObject* self, PyObject* args);
PyObject* pyx_is_transitive(PyObject* self, PyObject* args);
PyObject* pyx_is_preorder(PyObject* self, PyObject* args);
PyObject* pyx_is_partial_order(PyObject* self, PyObject* args);
PyObject* pyx_reflexive_closure(PyObject* self, PyObject* args);
PyObject* pyx_transitive_closure(PyObject* self, PyObject* args);

// --- Lex / Fuzzy (pywrapper_lex_fuzzy.cpp) ----------------------------------
PyObject* pyx_lex_separation(PyObject* self, PyObject* args);
PyObject* pyx_lex_mrp(PyObject* self, PyObject* args);
PyObject* pyx_fuzzy_separation(PyObject* self, PyObject* args);
PyObject* pyx_fuzzy_inbetweenness(PyObject* self, PyObject* args);

// --- Riduzione dimensionale (pywrapper_dimred.cpp) --------------------------
PyObject* pyx_dimensionality_reduction(PyObject* self, PyObject* args);
PyObject* pyx_bidimensional_poset_representation(PyObject* self, PyObject* args);

// --- First Order Dominance (pywrapper_fod.cpp) ------------------------------
PyObject* pyx_first_order_dominance_analysis(PyObject* self, PyObject* args);
