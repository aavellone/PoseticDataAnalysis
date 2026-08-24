/**
 * @file pymodule.cpp
 * @brief Definizione del modulo di estensione nativo `poseticDataAnalysis._core`.
 *
 * @details Equivalente della tabella di registrazione R (useDynLib / R_CallMethodDef):
 * espone tutte le PyCFunction come metodi del modulo. Il layer Python in
 * poseticDataAnalysis costruisce sopra queste primitive un'API idiomatica.
 */

#include "pywrapper.h"

#define PY_SSIZE_T_CLEAN
#include <Python.h>

// ===========================================================================
// Tabella dei metodi del modulo
// ===========================================================================

static PyMethodDef kCoreMethods[] = {
    // --- Costruttori --------------------------------------------------------
    {"build_poset", pyx_build_poset, METH_VARARGS,
     "build_poset(elements, relations) -> handle. relations: sequenza di coppie "
     "(a, b) con significato a < b."},
    {"build_linear_poset", pyx_build_linear_poset, METH_VARARGS,
     "build_linear_poset(elements) -> handle. Catena totale nell'ordine dato."},
    {"build_product_poset", pyx_build_product_poset, METH_VARARGS,
     "build_product_poset(posets) -> handle. Prodotto diretto (componentwise)."},
    {"build_lexicographic_product_poset", pyx_build_lexicographic_product_poset,
     METH_VARARGS,
     "build_lexicographic_product_poset(posets) -> handle. Prodotto lessicografico."},
    {"build_intersection_poset", pyx_build_intersection_poset, METH_VARARGS,
     "build_intersection_poset(posets) -> handle. Intersezione degli ordini."},
    {"build_linear_sum_poset", pyx_build_linear_sum_poset, METH_VARARGS,
     "build_linear_sum_poset(posets) -> handle. Somma lineare (ordinale)."},
    {"build_disjoint_sum_poset", pyx_build_disjoint_sum_poset, METH_VARARGS,
     "build_disjoint_sum_poset(posets) -> handle. Somma disgiunta (antichain di blocchi)."},
    {"build_lifting_poset", pyx_build_lifting_poset, METH_VARARGS,
     "build_lifting_poset(poset, new_element) -> handle. Aggiunge un minimo."},
    {"build_binary_variable_poset", pyx_build_binary_variable_poset, METH_VARARGS,
     "build_binary_variable_poset(variables) -> handle."},
    {"build_fence_poset", pyx_build_fence_poset, METH_VARARGS,
     "build_fence_poset(elements, orientation) -> handle. POSet a zig-zag."},
    {"build_crown_poset", pyx_build_crown_poset, METH_VARARGS,
     "build_crown_poset(elements_1, elements_2) -> handle."},
    {"build_dual_poset", pyx_build_dual_poset, METH_VARARGS,
     "build_dual_poset(poset) -> handle. Ordine invertito."},

    // --- Struttura ----------------------------------------------------------
    {"size", pyx_size, METH_VARARGS, "size(poset) -> int."},
    {"elements", pyx_elements, METH_VARARGS, "elements(poset) -> list[str]."},
    {"incidence_matrix", pyx_incidence_matrix, METH_VARARGS,
     "incidence_matrix(poset) -> list[list[int]] (matrice di adiacenza)."},
    {"cover_matrix", pyx_cover_matrix, METH_VARARGS,
     "cover_matrix(poset) -> list[list[int]] (diagramma di Hasse)."},
    {"order_relation", pyx_order_relation, METH_VARARGS,
     "order_relation(poset) -> list[tuple[str, str]] (u < v)."},
    {"cover_relation", pyx_cover_relation, METH_VARARGS,
     "cover_relation(poset) -> list[tuple[str, str]] (archi di copertura)."},
    {"incomparabilities", pyx_incomparabilities, METH_VARARGS,
     "incomparabilities(poset) -> list[tuple[str, str]]."},

    // --- Confronti pairwise -------------------------------------------------
    {"is_dominated_by", pyx_is_dominated_by, METH_VARARGS,
     "is_dominated_by(poset, a, b) -> list[bool] (a_k <= b_k)."},
    {"dominates", pyx_dominates, METH_VARARGS,
     "dominates(poset, a, b) -> list[bool] (a_k >= b_k)."},
    {"is_comparable_with", pyx_is_comparable_with, METH_VARARGS,
     "is_comparable_with(poset, a, b) -> list[bool]."},
    {"is_incomparable_with", pyx_is_incomparable_with, METH_VARARGS,
     "is_incomparable_with(poset, a, b) -> list[bool]."},

    // --- Upset / Downset ----------------------------------------------------
    {"upset_of", pyx_upset_of, METH_VARARGS, "upset_of(poset, elements) -> list[str]."},
    {"downset_of", pyx_downset_of, METH_VARARGS, "downset_of(poset, elements) -> list[str]."},
    {"is_upset", pyx_is_upset, METH_VARARGS, "is_upset(poset, elements) -> bool."},
    {"is_downset", pyx_is_downset, METH_VARARGS, "is_downset(poset, elements) -> bool."},

    // --- (In)comparabilità di un elemento -----------------------------------
    {"comparability_set_of", pyx_comparability_set_of, METH_VARARGS,
     "comparability_set_of(poset, element) -> list[str]."},
    {"incomparability_set_of", pyx_incomparability_set_of, METH_VARARGS,
     "incomparability_set_of(poset, element) -> list[str]."},

    // --- Estremali ----------------------------------------------------------
    {"maximals", pyx_maximals, METH_VARARGS, "maximals(poset) -> list[str]."},
    {"minimals", pyx_minimals, METH_VARARGS, "minimals(poset) -> list[str]."},
    {"is_maximal", pyx_is_maximal, METH_VARARGS, "is_maximal(poset, element) -> bool."},
    {"is_minimal", pyx_is_minimal, METH_VARARGS, "is_minimal(poset, element) -> bool."},
    {"meet", pyx_meet, METH_VARARGS, "meet(poset, elements) -> str | None."},
    {"join", pyx_join, METH_VARARGS, "join(poset, elements) -> str | None."},

    // --- Estensione ---------------------------------------------------------
    {"is_extension_of", pyx_is_extension_of, METH_VARARGS,
     "is_extension_of(poset_1, poset_2) -> bool."},

    // --- Generatori di estensioni lineari -----------------------------------
    {"build_le_generator", pyx_build_le_generator, METH_VARARGS,
     "build_le_generator(poset) -> handle (esatto, Tree-of-Ideals)."},
    {"build_bubley_dyer_le_generator", pyx_build_bubley_dyer_le_generator, METH_VARARGS,
     "build_bubley_dyer_le_generator(poset, seed) -> (handle, seed_str)."},
    {"leg_get", pyx_leg_get, METH_VARARGS,
     "leg_get(gen, from_start, quante, output_every) -> list[list[str]]."},
    {"leg_bubley_dyer_get", pyx_leg_bubley_dyer_get, METH_VARARGS,
     "leg_bubley_dyer_get(gen, from_start, quante, error, output_every) -> list[list[str]]."},

    // --- Valutazione --------------------------------------------------------
    {"exact_mrp", pyx_exact_mrp, METH_VARARGS,
     "exact_mrp(poset, output_every) -> {mrp, elements, n}."},
    {"build_bubley_dyer_mrp_generator", pyx_build_bubley_dyer_mrp_generator, METH_VARARGS,
     "build_bubley_dyer_mrp_generator(poset, seed) -> (handle, seed_str)."},
    {"bubley_dyer_mrp", pyx_bubley_dyer_mrp, METH_VARARGS,
     "bubley_dyer_mrp(gen, quante, error, output_every) -> {mrp, elements, n}."},
    {"build_bubley_dyer_evaluation_generator", pyx_build_bubley_dyer_evaluation_generator,
     METH_VARARGS,
     "build_bubley_dyer_evaluation_generator(poset, functions, seed) -> (handle, seed_str)."},
    {"bubley_dyer_evaluation", pyx_bubley_dyer_evaluation, METH_VARARGS,
     "bubley_dyer_evaluation(gen, quante, error, output_every) -> {results, n}."},
    {"exact_evaluation", pyx_exact_evaluation, METH_VARARGS,
     "exact_evaluation(poset, functions, output_every) -> {results, n}."},
    {"bls_dominance", pyx_bls_dominance, METH_VARARGS,
     "bls_dominance(poset, relative) -> {matrix, elements}."},

    // --- Proprietà delle relazioni ------------------------------------------
    {"is_reflexive", pyx_is_reflexive, METH_VARARGS, "is_reflexive(set, rel) -> bool."},
    {"is_symmetric", pyx_is_symmetric, METH_VARARGS, "is_symmetric(rel) -> bool."},
    {"is_antisymmetric", pyx_is_antisymmetric, METH_VARARGS, "is_antisymmetric(rel) -> bool."},
    {"is_transitive", pyx_is_transitive, METH_VARARGS, "is_transitive(rel) -> bool."},
    {"is_preorder", pyx_is_preorder, METH_VARARGS, "is_preorder(set, rel) -> bool."},
    {"is_partial_order", pyx_is_partial_order, METH_VARARGS, "is_partial_order(set, rel) -> bool."},
    {"reflexive_closure", pyx_reflexive_closure, METH_VARARGS,
     "reflexive_closure(set, rel) -> list[tuple[str, str]]."},
    {"transitive_closure", pyx_transitive_closure, METH_VARARGS,
     "transitive_closure(rel) -> list[tuple[str, str]]."},

    // --- Lex / Fuzzy --------------------------------------------------------
    {"lex_separation", pyx_lex_separation, METH_VARARGS,
     "lex_separation(modalita) -> {symmetric, asymmetricLower, asymmetricUpper, vertical, horizontal, labels}."},
    {"lex_mrp", pyx_lex_mrp, METH_VARARGS, "lex_mrp(modalita) -> {mrp, labels}."},
    {"fuzzy_separation", pyx_fuzzy_separation, METH_VARARGS,
     "fuzzy_separation(dominance, elements, norm, quali) -> {<type>: matrix, elements}."},
    {"fuzzy_inbetweenness", pyx_fuzzy_inbetweenness, METH_VARARGS,
     "fuzzy_inbetweenness(dominance, elements, norm, quali) -> {<type>: 3d array, elements}."},

    // --- Riduzione dimensionale ---------------------------------------------
    {"dimensionality_reduction", pyx_dimensionality_reduction, METH_VARARGS,
     "dimensionality_reduction(profile, weights, loss, lpom_strategy, output_every, thread_percentage)."},
    {"bidimensional_poset_representation", pyx_bidimensional_poset_representation, METH_VARARGS,
     "bidimensional_poset_representation(profile, weights, loss, lpom_strategy, variable_priority)."},

    // --- First Order Dominance ----------------------------------------------
    {"first_order_dominance_analysis", pyx_first_order_dominance_analysis, METH_VARARGS,
     "first_order_dominance_analysis(posets, freq_matrix, row_labels, col_labels, metrics, "
     "subpopulation_count, total_bins, count, seed, output_every, sep, linear_extensions, n_threads)."},

    {nullptr, nullptr, 0, nullptr}  // sentinella
};

// ===========================================================================
// Definizione e inizializzazione del modulo
// ===========================================================================

static struct PyModuleDef kCoreModule = {
    PyModuleDef_HEAD_INIT,
    "poseticDataAnalysis._core",
    "Native C++20 core for posetic data analysis (CPython C-API binding).",
    -1,
    kCoreMethods,
    nullptr, nullptr, nullptr, nullptr};

PyMODINIT_FUNC PyInit__core(void) {
    return PyModule_Create(&kCoreModule);
}
