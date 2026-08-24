#pragma once

// Dipendenze del progetto
#include "tensor.h"
#include "poset_wrapper.h"
#include "linear_generator_wrapper.h"
#include "my_exception.h"
#include "rwrapper.h"

// Standard Library C++
#include <algorithm>       // Per std::copy
#include <array>           // Per std::array
#include <charconv>        // Per std::from_chars (parsing dei seed a 64 bit)
#include <cstdint>         // Per std::uint8_t, std::uint64_t
#include <limits>          // Per std::numeric_limits
#include <list>            // Per std::list
#include <memory>          // Per std::unique_ptr
#include <optional>        // Per std::optional
#include <stdexcept>       // Per std::invalid_argument, std::runtime_error
#include <string>          // Per std::string
#include <string_view>     // Per std::string_view
#include <unordered_map>   // Per std::unordered_map
#include <utility>         // Per std::pair, std::move
#include <vector>          // Per std::vector

// Intestazioni R
#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>

/**
 * ============================================================================
 * NAMESPACE: RConvert
 * @brief Utility HPC per l'estrazione dati da R (SEXP) verso C++20.
 * @details Converte in modo sicuro i tipi R (vettori, matrici, puntatori) in
 * strutture C++ native. Non acquisisce memory ownership, ma esegue copie
 * o restituisce puntatori osservatori (`const T*`) in O(1).
 * ============================================================================
 */
namespace RConvert {

/**
 * @brief Converte un vettore di stringhe R in uno std::vector<std::string>.
 */
[[nodiscard]] inline std::vector<std::string> ToStringVector(SEXP string_vec_r) {
    if (!Rf_isString(string_vec_r)) {
        throw std::invalid_argument("RConvert::ToStringVector - L'argomento non è un vettore di stringhe.");
    }
    
    const R_len_t len = Rf_length(string_vec_r);
    std::vector<std::string> result;
    result.reserve(static_cast<std::size_t>(len));
    
    for (R_len_t i = 0; i < len; ++i) {
        result.emplace_back(CHAR(STRING_ELT(string_vec_r, i)));
    }
    return result;
}

/**
 * @brief Estrae una lista di archi (comparabilità) da una matrice Nx2 stringhe in R.
 */
[[nodiscard]] inline std::vector<std::pair<std::string, std::string>> ToPairVector(SEXP matrix_r) {
    if (!Rf_isMatrix(matrix_r) || TYPEOF(matrix_r) != STRSXP) {
        throw std::invalid_argument("RConvert::ToPairVector - Richiesta una matrice di stringhe.");
    }
    
    const R_len_t rows = Rf_nrows(matrix_r);
    const R_len_t cols = Rf_ncols(matrix_r);
    
    if (cols != 2) {
        throw std::invalid_argument("RConvert::ToPairVector - La matrice deve avere esattamente 2 colonne.");
    }
    
    std::vector<std::pair<std::string, std::string>> result;
    result.reserve(static_cast<std::size_t>(rows));
    
    for (R_len_t i = 0; i < rows; ++i) {
        std::string from = CHAR(STRING_ELT(matrix_r, i));
        std::string to = CHAR(STRING_ELT(matrix_r, i + rows)); // R usa il column-major order
        result.emplace_back(std::move(from), std::move(to));
    }
    return result;
}

/**
 * @brief Estrae e valida il raw pointer di un oggetto C++ gestito da un ExternalPtr R.
 * (ZERO OWNERSHIP)
 * @details Un EXTPTRSXP con indirizzo nullo e' cio' che R produce quando un external
 * pointer viene deserializzato (workspace .RData ricaricato da RStudio, saveRDS,
 * invio a processi paralleli): l'involucro R sopravvive, l'oggetto C++ no.
 * Senza questo controllo il chiamante dereferenzierebbe nullptr facendo abortire
 * l'INTERO processo R (segfault), con perdita del lavoro non salvato: ogni punto
 * d'ingresso che riceve un ExternalPtr deve passare da qui.
 * @tparam T Tipo dell'oggetto C++ puntato (puo' essere const-qualificato).
 * @param ext_ptr_r External pointer ricevuto da R.
 * @param object_name Nome leggibile dell'oggetto, inserito nel messaggio d'errore.
 * @return Il puntatore, garantito non nullo.
 * @throws std::invalid_argument Se l'argomento non e' un ExternalPtr.
 * @throws std::runtime_error Se il puntatore e' nullo.
 * @warning I messaggi non devono contenere ':' oltre a quelli del prefisso
 * 'RConvert::ToWrapPtr': diversi wrapper R mostrano all'utente solo il frammento
 * successivo all'ultimo ':' (strsplit + ultimo elemento).
 */
template <typename T>
[[nodiscard]] inline T* ToWrapPtr(SEXP ext_ptr_r, const char* object_name) {
    if (TYPEOF(ext_ptr_r) != EXTPTRSXP) {
        throw std::invalid_argument(
            std::string("RConvert::ToWrapPtr - L'argomento non e' un ExternalPtr (")
            + object_name + ").");
    }

    void* raw_ptr = R_ExternalPtrAddr(ext_ptr_r);
    if (raw_ptr == nullptr) {
        throw std::runtime_error(
            std::string("RConvert::ToWrapPtr - ") + object_name +
            " non e' piu' valido (puntatore esterno nullo). Gli oggetti C++ non"
            " sopravvivono al salvataggio e al ricaricamento della sessione R"
            " (.RData, saveRDS) ne' all'invio a processi paralleli, e vanno"
            " ricostruiti con la corrispondente funzione del package.");
    }

    return static_cast<T*>(raw_ptr);
}

/**
 * @brief Estrae il raw pointer osservatore dal wrapper R. (ZERO OWNERSHIP)
 * @details Ideale per il passaggio di const reference alle factory C++.
 * Delega a ToWrapPtr, che valida tipo e indirizzo del puntatore.
 */
[[nodiscard]] inline const POSetWrap* ToPOSetWrap(SEXP ext_ptr_r) {
    return ToWrapPtr<const POSetWrap>(ext_ptr_r, "Il POSet");
}


/**
 * @brief Converts R elements and comparabilities into a C++ Adjacency Matrix.
 * * @details Extracts nodes and edges directly from R SEXP objects. Uses C++20
 * std::string_view for zero-allocation hash map lookups, maximizing HPC performance.
 * It natively exploits R's column-major matrix memory layout.
 * * @param elements_r R character vector containing the node names.
 * @param comparabilities_r R character matrix (N x 2) representing directed edges.
 * @param out_elements Reference to a vector to populate with node names.
 * @return Tensor<std::uint8_t, 2> The resulting adjacency matrix.
 * @throws std::invalid_argument If the inputs do not have the expected R types.
 * @throws MyException If a duplicate node name or an unknown edge endpoint is found.
 */
[[nodiscard]] inline Tensor<std::uint8_t, 2> ToAdjacencyMatrix(
                                                               SEXP elements_r,
                                                               SEXP comparabilities_r,
                                                               std::vector<std::string>& out_elements) {
    
    // 0. Validazione dei tipi R (evita segfault con input errati)
    if (!Rf_isString(elements_r)) {
        throw std::invalid_argument("RConvert::ToAdjacencyMatrix - 'elements_r' deve essere un vettore di stringhe.");
    }
    if (!Rf_isMatrix(comparabilities_r) || TYPEOF(comparabilities_r) != STRSXP ||
        Rf_ncols(comparabilities_r) != 2) {
        throw std::invalid_argument("RConvert::ToAdjacencyMatrix - 'comparabilities_r' deve essere una matrice di stringhe con 2 colonne.");
    }
    
    // 1. Estrazione nodi diretta da R
    const R_len_t num_nodes = Rf_length(elements_r);
    const auto sz_nodes = static_cast<std::size_t>(num_nodes);
    
    out_elements.clear();
    out_elements.reserve(sz_nodes);
    
    // std::string_view azzera le allocazioni heap per le chiavi della mappa!
    std::unordered_map<std::string_view, std::size_t> name_to_id;
    name_to_id.reserve(sz_nodes);
    
    for (R_len_t i = 0; i < num_nodes; ++i) {
        const char* node_name = CHAR(STRING_ELT(elements_r, i));
        
        // emplace_back alloca la std::string una sola volta per l'output
        out_elements.emplace_back(node_name);
        
        // emplace usa lo string_view senza copiare la memoria della stringa
        const auto [it, inserted] = name_to_id.emplace(node_name, static_cast<std::size_t>(i));
        if (!inserted) {
            throw MyException(std::string("ToAdjacencyMatrix: duplicate element name '") + node_name + "'.");
        }
    }
    
    // 2. Creazione della matrice inizializzata a 0
    Tensor<std::uint8_t, 2> adj(std::array<std::uint64_t, 2>{sz_nodes, sz_nodes}, 0);
    
    // 3. Estrazione degli archi bypassando vettori intermedi
    const R_len_t num_edges = Rf_nrows(comparabilities_r);
    
    for (R_len_t i = 0; i < num_edges; ++i) {
        // Accesso O(1) diretto alla memoria column-major di R.
        const char* from_str = CHAR(STRING_ELT(comparabilities_r, i));
        const char* to_str   = CHAR(STRING_ELT(comparabilities_r, i + num_edges));
        
        // .find() ora accetta const char* (convertito in string_view implicitamente) senza allocare
        auto it1 = name_to_id.find(from_str);
        auto it2 = name_to_id.find(to_str);
        
        if (it1 != name_to_id.end() && it2 != name_to_id.end()) {
            adj(it1->second, it2->second) = 1;
        } else {
            throw MyException("ToAdjacencyMatrix: edge refers to an unknown node.");
        }
    }
    
    return adj;
}

/**
 * @brief Infers nodes and builds an Adjacency Matrix directly from an edge list.
 * @details Extracts unique nodes on the fly from an R character matrix (N x 2) of edges.
 * Utilizes C++20 std::string_view to eliminate heap allocations for map keys, and
 * try_emplace to guarantee a single hash computation per node observation.
 * @param comparabilities_r R character matrix (N x 2) representing directed edges.
 * @param out_elements Reference to a vector that will be populated with the inferred unique node names.
 * @return Tensor<std::uint8_t, 2> The resulting adjacency matrix.
 * @throws std::invalid_argument If the input is not an N x 2 character matrix.
 */
[[nodiscard]] inline Tensor<std::uint8_t, 2> ToAdjacencyMatrixFromEdges(
                                                                        SEXP comparabilities_r,
                                                                        std::vector<std::string>& out_elements) {
    
    // Validazione del tipo R (evita segfault con input errati)
    if (!Rf_isMatrix(comparabilities_r) || TYPEOF(comparabilities_r) != STRSXP ||
        Rf_ncols(comparabilities_r) != 2) {
        throw std::invalid_argument("RConvert::ToAdjacencyMatrixFromEdges - 'comparabilities_r' deve essere una matrice di stringhe con 2 colonne.");
    }
    
    out_elements.clear();
    
    // HPC C++20: Zero-allocation string views for map keys.
    std::unordered_map<std::string_view, std::size_t> name_to_id;
    
    const R_len_t num_edges = Rf_nrows(comparabilities_r);
    
    // Euristicamente, pre-allochiamo spazio per evitare il rehashing della mappa.
    // Nel caso peggiore (grafo disconnesso con solo archi isolati) avremo num_edges * 2 nodi,
    // ma tipicamente il numero di nodi è <= num_edges.
    name_to_id.reserve(static_cast<std::size_t>(num_edges));
    
    // Passata 1: Identifica i nodi unici con inserimento a lookup singolo (try_emplace)
    for (R_len_t i = 0; i < num_edges; ++i) {
        const char* from_str = CHAR(STRING_ELT(comparabilities_r, i));
        const char* to_str   = CHAR(STRING_ELT(comparabilities_r, i + num_edges));
        
        // try_emplace tenta l'inserimento. Ritorna una coppia:
        // .first è l'iteratore all'elemento (nuovo o già esistente)
        // .second è un booleano (true se è stato appena inserito, false se esisteva già)
        auto [it_from, inserted_from] = name_to_id.try_emplace(from_str, out_elements.size());
        if (inserted_from) {
            out_elements.emplace_back(from_str);
        }
        
        auto [it_to, inserted_to] = name_to_id.try_emplace(to_str, out_elements.size());
        if (inserted_to) {
            out_elements.emplace_back(to_str);
        }
    }
    
    // Passata 2: Ora che sappiamo la dimensione esatta, allochiamo la matrice HPC
    const std::size_t num_nodes = out_elements.size();
    Tensor<std::uint8_t, 2> adj(std::array<std::uint64_t, 2>{num_nodes, num_nodes}, 0);
    
    for (R_len_t i = 0; i < num_edges; ++i) {
        const char* from_str = CHAR(STRING_ELT(comparabilities_r, i));
        const char* to_str   = CHAR(STRING_ELT(comparabilities_r, i + num_edges));
        
        // Usiamo .find() invece di operator[] per non innescare controlli di inserimento inutili,
        // dato che siamo assolutamente certi che i nodi siano già nella mappa.
        std::size_t from_id = name_to_id.find(from_str)->second;
        std::size_t to_id   = name_to_id.find(to_str)->second;
        
        adj(from_id, to_id) = 1;
    }
    
    return adj;
}

/**
 * @brief Estrae una lista di puntatori POSetWrap da una lista R (VECSXP).
 * @param list_r SEXP di tipo lista contenente ExternalPtr.
 * @return std::vector<const POSetWrap*> Vettore di puntatori osservatori.
 * @throw std::invalid_argument se l'input non è una lista.
 */
[[nodiscard]] inline std::vector<const POSetWrap*> ToPOSetWrapVector(SEXP list_r) {
    if (TYPEOF(list_r) != VECSXP) {
        throw std::invalid_argument("RConvert::ToPOSetWrapVector - L'argomento deve essere una lista.");
    }
    const R_len_t len = Rf_length(list_r);
    std::vector<const POSetWrap*> result;
    result.reserve(static_cast<std::size_t>(len));
    
    for (R_len_t i = 0; i < len; ++i) {
        // Riutilizza ToPOSetWrap per ogni elemento della lista
        result.push_back(ToPOSetWrap(VECTOR_ELT(list_r, i)));
    }
    return result;
}

/**
 * @brief Extracts a single boolean from an R logical scalar.
 * @throws std::invalid_argument If the input is not a logical vector of
 * length >= 1, or if the value is NA.
 */
inline bool ToBool(SEXP sexp) {
    if (!Rf_isLogical(sexp) || Rf_length(sexp) < 1) {
        throw std::invalid_argument("RConvert::ToBool - Expected a logical vector of length >= 1.");
    }
    const int val = LOGICAL(sexp)[0];
    if (val == NA_LOGICAL) {
        throw std::invalid_argument("RConvert::ToBool - NA is not allowed.");
    }
    return val != 0;
}

/**
 * @brief Extracts a single integer from an R integer or double scalar.
 * @details Accepts both INTSXP and REALSXP because R numeric literals
 * (e.g. `100`) are doubles by default. Rejects NA and out-of-range values.
 * @throws std::invalid_argument On wrong type, empty vector, NA, or overflow.
 */
inline int ToInt(SEXP sexp) {
    if (Rf_length(sexp) < 1) {
        throw std::invalid_argument("RConvert::ToInt - Expected a numeric vector of length >= 1.");
    }
    if (Rf_isInteger(sexp)) {
        const int val = INTEGER(sexp)[0];
        if (val == NA_INTEGER) {
            throw std::invalid_argument("RConvert::ToInt - NA is not allowed.");
        }
        return val;
    }
    if (Rf_isReal(sexp)) {
        const double val = REAL(sexp)[0];
        if (ISNAN(val) ||
            val < static_cast<double>(std::numeric_limits<int>::min()) ||
            val > static_cast<double>(std::numeric_limits<int>::max())) {
            throw std::invalid_argument("RConvert::ToInt - Value is NA or out of range for int.");
        }
        return static_cast<int>(val);
    }
    throw std::invalid_argument("RConvert::ToInt - Expected an integer or numeric scalar.");
}

/**
 * @brief Extracts a single double from an R double or integer scalar.
 * @throws std::invalid_argument On wrong type or empty vector.
 */
inline double ToDouble(SEXP sexp) {
    if (Rf_length(sexp) < 1) {
        throw std::invalid_argument("RConvert::ToDouble - Expected a numeric vector of length >= 1.");
    }
    if (Rf_isReal(sexp)) {
        return REAL(sexp)[0];
    }
    if (Rf_isInteger(sexp)) {
        const int val = INTEGER(sexp)[0];
        if (val == NA_INTEGER) {
            throw std::invalid_argument("RConvert::ToDouble - NA is not allowed.");
        }
        return static_cast<double>(val);
    }
    throw std::invalid_argument("RConvert::ToDouble - Expected a numeric scalar.");
}

/**
 * @brief Extracts an unsigned 64-bit integer from an R SEXP, returning an optional value.
 *
 * @details Safely checks if the provided R object contains any elements.
 * Accepts both INTSXP and REALSXP (R numeric literals are doubles by default).
 * NULL, empty vectors and NA yield std::nullopt. Negative or out-of-range
 * values throw, so they can never silently overflow into a huge uint64_t.
 *
 * @param sexp The R object (expected to be an integer/numeric vector or NULL).
 * @return std::optional<std::uint64_t> The extracted value if valid, otherwise std::nullopt.
 * @throws std::invalid_argument On negative values, out-of-range values, or wrong type.
 */
inline std::optional<std::uint64_t> ToOptionalUInt(SEXP sexp) {
    if (Rf_isNull(sexp) || Rf_length(sexp) == 0) {
        return std::nullopt;
    }
    if (Rf_isInteger(sexp)) {
        const int val = INTEGER(sexp)[0];
        if (val == NA_INTEGER) {
            return std::nullopt;
        }
        if (val < 0) {
            throw std::invalid_argument("RConvert::ToOptionalUInt - Negative values are not allowed.");
        }
        return static_cast<std::uint64_t>(val);
    }
    if (Rf_isReal(sexp)) {
        const double val = REAL(sexp)[0];
        if (ISNAN(val)) {
            return std::nullopt;
        }
        if (val < 0.0 ||
            val >= static_cast<double>(std::numeric_limits<std::uint64_t>::max())) {
            throw std::invalid_argument("RConvert::ToOptionalUInt - Value out of range for uint64.");
        }
        return static_cast<std::uint64_t>(val);
    }
    throw std::invalid_argument("RConvert::ToOptionalUInt - Expected an integer or numeric scalar (or NULL).");
}

/**
 * @brief Converte l'elemento i-esimo di 'seed_r' in un seme a 64 bit.
 * @details I semi del motore C++ sono std::uint64_t, un intervallo che R non sa
 * rappresentare: 'integer' e' a 32 bit e 'double' e' esatto solo fino a 2^53.
 * Il canale ufficiale e' percio' la STRINGA di cifre decimali, che attraversa il
 * confine R/C++ senza perdita di precisione. I wrapper R convertono in stringa
 * qualunque seme numerico prima della .Call. I tipi numerici restano accettati
 * qui per le .Call diverte, ma con i limiti di rappresentazione di R.
 * @param seed_r Vettore character (canale preferito), integer o numeric.
 * @param i Indice dell'elemento da convertire.
 * @return Il seme come std::uint64_t.
 * @throws std::invalid_argument Per tipo errato, NA, formato non numerico o
 * valore fuori dall'intervallo di uint64.
 * @warning I messaggi non devono contenere ':' oltre a quelli del prefisso.
 */
[[nodiscard]] inline std::uint64_t ToSeedAt(SEXP seed_r, R_xlen_t i) {
    switch (TYPEOF(seed_r)) {
        case STRSXP: {
            SEXP elem_r = STRING_ELT(seed_r, i);
            if (elem_r == NA_STRING) {
                throw std::invalid_argument("RConvert::ToSeedAt - NA is not a valid seed.");
            }
            const char* first = CHAR(elem_r);
            const char* last  = first + Rf_xlength(elem_r);

            // from_chars e' strict by design: niente spazi, segni o notazione
            // scientifica, solo cifre decimali. Nessuna dipendenza dal locale.
            std::uint64_t value = 0;
            const auto res = std::from_chars(first, last, value);
            if (res.ec == std::errc::result_out_of_range) {
                throw std::invalid_argument(
                    "RConvert::ToSeedAt - Seed out of range for a 64-bit unsigned integer.");
            }
            if (res.ec != std::errc{} || res.ptr != last) {
                throw std::invalid_argument(
                    "RConvert::ToSeedAt - Seed must be a string of decimal digits.");
            }
            return value;
        }
        case INTSXP: {
            const int value = INTEGER(seed_r)[i];
            if (value == NA_INTEGER || value < 0) {
                throw std::invalid_argument(
                    "RConvert::ToSeedAt - Seeds must be non-negative and not NA.");
            }
            return static_cast<std::uint64_t>(value);
        }
        case REALSXP: {
            const double value = REAL(seed_r)[i];
            if (ISNAN(value) || value < 0.0) {
                throw std::invalid_argument(
                    "RConvert::ToSeedAt - Seeds must be non-negative and not NA.");
            }
            if (value >= static_cast<double>(std::numeric_limits<std::uint64_t>::max())) {
                throw std::invalid_argument(
                    "RConvert::ToSeedAt - Seed out of range for a 64-bit unsigned integer.");
            }
            return static_cast<std::uint64_t>(value);
        }
        default:
            throw std::invalid_argument(
                "RConvert::ToSeedAt - 'seed' must be a character, integer or numeric vector.");
    }
}

/**
 * @brief Estrae un seme scalare opzionale a 64 bit.
 * @details NULL o vettore vuoto significano "nessun seme" (il chiamante ne
 * genera uno casuale). Vedi ToSeedAt per i tipi accettati.
 */
[[nodiscard]] inline std::optional<std::uint64_t> ToOptionalSeed(SEXP seed_r) {
    if (Rf_isNull(seed_r) || Rf_xlength(seed_r) == 0) {
        return std::nullopt;
    }
    return ToSeedAt(seed_r, 0);
}

/**
 * @brief Converte una matrice character R in un vettore di LinearExtension.
 * @details Ogni COLONNA della matrice rappresenta una estensione lineare come
 * sequenza di nomi di elementi del poset (posizione 0 = prima riga). La scelta
 * delle colonne sfrutta il layout column-major di R: la lettura di ogni
 * estensione è contigua in memoria. La compatibilità con l'ordine parziale
 * NON viene verificata qui: se ne occupa il costruttore di LEGBubleyDyer.
 * @param matrix_r Matrice character n x K (n = poset.size(), K = numero di LE).
 * @param poset Poset di riferimento per la mappatura nome -> id.
 * @throws std::invalid_argument per tipo/dimensioni errati o NA.
 * @throws MyException se un nome non appartiene al poset (da GetElementId).
 */
[[nodiscard]] inline std::vector<LinearExtension> ToLinearExtensions(SEXP matrix_r, const POSet& poset) {
    if (!Rf_isMatrix(matrix_r) || !Rf_isString(matrix_r)) {
        throw std::invalid_argument("RConvert::ToLinearExtensions - Expected a character matrix "
                                    "(each column is one linear extension of element names).");
    }
    const std::size_t n_rows = static_cast<std::size_t>(Rf_nrows(matrix_r));
    const std::size_t n_cols = static_cast<std::size_t>(Rf_ncols(matrix_r));
    const std::uint64_t n = poset.size();

    if (n_rows != n) {
        throw std::invalid_argument("RConvert::ToLinearExtensions - The matrix must have one row "
                                    "per poset element (rows = poset size, columns = extensions).");
    }
    if (n_cols == 0) {
        throw std::invalid_argument("RConvert::ToLinearExtensions - The matrix must have at least one column.");
    }

    std::vector<LinearExtension> les;
    les.reserve(n_cols);
    for (std::size_t c = 0; c < n_cols; ++c) {
        LinearExtension le(n);
        for (std::size_t r = 0; r < n_rows; ++r) {
            SEXP name_r = STRING_ELT(matrix_r, static_cast<R_xlen_t>(c * n_rows + r));
            if (name_r == NA_STRING) {
                throw std::invalid_argument("RConvert::ToLinearExtensions - NA is not allowed.");
            }
            le.Set(r, poset.GetElementId(CHAR(name_r)));
        }
        les.push_back(std::move(le));
    }
    return les;
}

/**
 * @brief Costruisce i seed per K catene MCMC indipendenti.
 * @details 'seed_r' può essere:
 * - NULL: seed base casuale (da Random::GLOBAL), seed per-catena derivati;
 * - uno scalare: seed base, seed per-catena derivati con SplitMix64(base + k)
 *   (decorrelazione riproducibile: stesso base + stesso K => stesse catene);
 * - un vettore di lunghezza K: seed espliciti, uno per catena.
 */
[[nodiscard]] inline std::vector<std::uint64_t> ToChainSeeds(SEXP seed_r, std::uint64_t n_chains) {
    constexpr auto kSplitMix64 = [](std::uint64_t x) noexcept {
        x += 0x9E3779B97F4A7C15ULL;
        x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
        x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
        return x ^ (x >> 31);
    };

    std::vector<std::uint64_t> seeds;
    seeds.reserve(n_chains);

    const R_xlen_t len = Rf_isNull(seed_r) ? 0 : Rf_xlength(seed_r);

    // Vettore di seed espliciti (uno per catena).
    if (len > 1) {
        if (static_cast<std::uint64_t>(len) != n_chains) {
            throw std::invalid_argument("RConvert::ToChainSeeds - 'seed' must be NULL, a scalar, "
                                        "or a vector with exactly one seed per chain.");
        }
        for (R_xlen_t i = 0; i < len; ++i) {
            seeds.push_back(ToSeedAt(seed_r, i));
        }
        return seeds;
    }

    // Scalare o NULL: seed base + derivazione riproducibile per catena.
    const std::uint64_t base = (len == 0)
        ? Random::GLOBAL.RndNextInt(0, std::numeric_limits<std::uint64_t>::max())
        : ToSeedAt(seed_r, 0);

    for (std::uint64_t k = 0; k < n_chains; ++k) {
        seeds.push_back(kSplitMix64(base + k));
    }
    return seeds;
}

/**
 * @brief Extracts a double scalar from an R SEXP, returning an optional value.
 *
 * @details Safely checks if the provided R object contains any elements.
 * Accepts both REALSXP and INTSXP. NULL, empty vectors and NA yield
 * std::nullopt. By utilizing std::optional instead of heap-allocated
 * pointers, this function guarantees zero-overhead memory management and
 * clear semantics for missing or empty parameters.
 *
 * @param sexp The R object (expected to be a numeric vector or NULL).
 * @return std::optional<double> The extracted double if valid, otherwise std::nullopt.
 * @throws std::invalid_argument On wrong type.
 */
inline std::optional<double> ToOptionalDouble(SEXP sexp) {
    if (Rf_isNull(sexp) || Rf_length(sexp) == 0) {
        return std::nullopt;
    }
    if (Rf_isReal(sexp)) {
        const double val = REAL(sexp)[0];
        if (!ISNA(val)) {
            return val;
        }
        return std::nullopt;
    }
    if (Rf_isInteger(sexp)) {
        const int val = INTEGER(sexp)[0];
        if (val != NA_INTEGER) {
            return static_cast<double>(val);
        }
        return std::nullopt;
    }
    throw std::invalid_argument("RConvert::ToOptionalDouble - Expected a numeric scalar (or NULL).");
}

/**
 * @brief Parses an R list of strings / R functions into parallel vectors.
 *
 * Each element of @p functions_r is either:
 * - A string  → stored in @p internal_funcs by name.
 * - A function → stored in @p external_funcs; key set to `"RFunction"`.
 *
 * @param[in]  functions_r     R list (strings and/or R functions).
 * @param[out] internal_funcs  Internal C++ function name per entry.
 * @param[out] external_funcs  `SEXP` per entry (`R_NilValue` for string entries).
 * @throws std::invalid_argument If an entry is neither a function nor a non-empty character vector.
 */
inline void ToFunctions(SEXP functions_r,
                        std::vector<std::string>& internal_funcs,
                        std::vector<SEXP>&         external_funcs)
{
    const std::size_t n = static_cast<std::size_t>(Rf_length(functions_r));
    internal_funcs.assign(n, "");
    external_funcs.assign(n, R_NilValue);
    
    for (std::size_t k = 0; k < n; ++k) {
        SEXP func = VECTOR_ELT(functions_r, static_cast<R_xlen_t>(k));
        if (Rf_isFunction(func)) {
            external_funcs[k] = func;
            internal_funcs[k] = "RFunction";
        } else if (Rf_isString(func) && Rf_length(func) > 0) {
            internal_funcs[k] = R_CHAR(STRING_ELT(func, 0));
        } else {
            throw std::invalid_argument("RConvert::ToFunctions - Each entry must be a function or a non-empty string.");
        }
    }
}

/**
 * @brief Converts an R numeric matrix (REALSXP or INTSXP) into a C++ Tensor<double, 2>.
 *
 * @details Safely extracts dimensions and verifies the input type to prevent
 * segmentation faults. It reads the R column-major array into the C++
 * row-major matrix (row-outer, column-inner loops: strided reads from R,
 * sequential contiguous writes to the C++ buffer).
 * Integer matrices are converted element-wise (NA_integer_ becomes NA_real_).
 * Returns by value to leverage Return Value Optimization (RVO) or move semantics,
 * avoiding unnecessary smart pointer allocations.
 *
 * @param sexp The R object expected to be a numeric (double or integer) matrix.
 * @return Tensor<double, 2> The parsed C++ matrix returned by value.
 * @throws std::invalid_argument If the input is not a numeric matrix.
 */
[[nodiscard]] inline Tensor<double, 2> ToDoubleMatrix(SEXP sexp) {
    if (!Rf_isMatrix(sexp) ||
        (TYPEOF(sexp) != REALSXP && TYPEOF(sexp) != INTSXP)) {
        throw std::invalid_argument("ToDoubleMatrix: Input must be a numeric (double or integer) matrix.");
    }
    
    SEXP dim_attr = Rf_getAttrib(sexp, R_DimSymbol);
    const auto nrow = static_cast<std::size_t>(INTEGER(dim_attr)[0]);
    const auto ncol = static_cast<std::size_t>(INTEGER(dim_attr)[1]);
    
    Tensor<double, 2> result(std::array<std::uint64_t, 2>{nrow, ncol}, kUninitialized);
    
    if (TYPEOF(sexp) == REALSXP) {
        const double* mat_ptr = REAL(sexp);
        for (std::size_t row = 0; row < nrow; ++row) {
            for (std::size_t col = 0; col < ncol; ++col) {
                // Lettura da R (con salti) -> Scrittura su C++ (sequenziale e contigua)
                result(row, col) = mat_ptr[row + nrow * col];
            }
        }
    } else {
        const int* mat_ptr = INTEGER(sexp);
        for (std::size_t row = 0; row < nrow; ++row) {
            for (std::size_t col = 0; col < ncol; ++col) {
                const int v = mat_ptr[row + nrow * col];
                result(row, col) = (v == NA_INTEGER) ? NA_REAL : static_cast<double>(v);
            }
        }
    }
    
    return result;
}

/**
 * @brief Converts an R numeric into a Tensor<double, 1>.
 * @details Accepts both REALSXP and INTSXP (NA_integer_ becomes NA_real_).
 * NULL yields an empty tensor.
 * @param vec_r The R SEXP object to convert.
 * @return Tensor<double, 1> The extracted C++ vector.
 * @throws std::invalid_argument If the provided SEXP is not a valid numeric-like type.
 */
[[nodiscard]] inline Tensor<double, 1> ToDouble1DTensor(SEXP vec_r) {
    if (Rf_isNull(vec_r)) [[unlikely]] {
        return {};
    }
    
    const R_xlen_t len = Rf_xlength(vec_r);
    Tensor<double, 1> res({static_cast<std::size_t>(len)}, kUninitialized);
    
    if (Rf_isReal(vec_r)) [[likely]] {
        const double* data = REAL(vec_r);
        std::copy(data, data + len, res.begin());
    }
    else if (Rf_isInteger(vec_r)) {
        const int* data = INTEGER(vec_r);
        for (R_xlen_t i = 0; i < len; ++i) {
            res(static_cast<std::size_t>(i)) =
            (data[i] == NA_INTEGER) ? NA_REAL : static_cast<double>(data[i]);
        }
    }
    else [[unlikely]] {
        throw std::invalid_argument("RConvert::ToDouble1DTensor: The provided object is not a numeric vector (expected REALSXP or INTSXP).");
    }
    
    return res;
}

}


/**
 * ============================================================================
 * NAMESPACE: RCreate
 * @brief Utility HPC per il rimpacchettamento di tipi C++ verso R (SEXP).
 * @details Le funzioni trasferiscono la proprietà della memoria dal C++ al
 * Garbage Collector di R registrando appositi finalizzatori.
 * ============================================================================
 */
namespace RCreate {

/**
 * @brief Generic finalizer for C++ objects managed by R external pointers.
 * @details Invoked automatically by R's Garbage Collector when the EXTPTRSXP
 * goes out of scope in R. Safely deletes the underlying C++ object and clears
 * the R pointer to prevent double-free corruption.
 * @tparam T The type of the underlying C++ object.
 * @param ptr_r The R external pointer (EXTPTRSXP).
 */
template <typename T>
inline void PointerFinalizer(SEXP ptr_r) noexcept {
    if (TYPEOF(ptr_r) == EXTPTRSXP) {
        T* obj = static_cast<T*>(R_ExternalPtrAddr(ptr_r));
        if (obj != nullptr) {
            delete obj;                  // Distrugge l'oggetto C++ (richiama il distruttore di T)
            R_ClearExternalPtr(ptr_r);   // Setta il payload dell'EXTPTRSXP a NULL
        }
    }
}

/**
 * @brief Safely wraps any C++ object managed by a std::unique_ptr into an R ExternalPtr.
 * * @details Uses C++ templates to genericize the transfer of ownership from C++
 * to R's garbage collector. It creates the EXTPTRSXP and registers the correctly
 * typed finalizer BEFORE releasing the unique_ptr, so ownership is handed over
 * to R's GC only once the wrapping has fully succeeded.
 * * @tparam T The type of the C++ object being wrapped.
 * @param guard The active RProtectGuard for memory safety.
 * @param ptr The unique_ptr holding the object instance.
 * @return SEXP The R external pointer.
 */
template <typename T>
[[nodiscard]] inline SEXP WrapExternalPtr(RProtectGuard& guard, std::unique_ptr<T> ptr) {
    if (!ptr) return R_NilValue;
    
    // Protegge e crea il puntatore esterno di R. Il unique_ptr mantiene la
    // ownership finché la creazione non è andata a buon fine.
    SEXP ext_ptr = guard.Protect(R_MakeExternalPtr(ptr.get(), R_NilValue, R_NilValue));
    
    // Registrazione del distruttore custom per evitare Memory Leaks
    R_RegisterCFinalizerEx(ext_ptr, PointerFinalizer<T>, TRUE);
    
    // Solo ora il C++ cede la proprietà del puntatore al sistema di GC di R
    ptr.release();
    
    return ext_ptr;
}

/**
 * @brief Converts a C++ boolean into an R logical scalar.
 * @details Uses R's native `Rf_ScalarLogical` for optimized, single-step
 * allocation and initialization of a length-1 LGLSXP vector.
 * @param guard The active RProtectGuard for GC memory safety.
 * @param value The C++ boolean value to convert.
 * @return SEXP The resulting R logical scalar.
 */
[[nodiscard]] inline SEXP FromBool(RProtectGuard& guard, bool value) {
    // Rf_ScalarLogical è più veloce e idiomatico rispetto ad allocare
    // e assegnare manualmente tramite LOGICAL()[0].
    return guard.Protect(Rf_ScalarLogical(value ? TRUE : FALSE));
}

/**
 * @brief Converts a C++ integer into an R integer scalar.
 * @details Uses R's native `Rf_ScalarInteger` for optimized, single-step
 * allocation and initialization of a length-1 INTSXP vector.
 * @param guard The active RProtectGuard for GC memory safety.
 * @param value The C++ integer value to convert.
 * @return SEXP The resulting R integer scalar.
 */
[[nodiscard]] inline SEXP FromInt(RProtectGuard& guard, int value) {
    return guard.Protect(Rf_ScalarInteger(value));
}

/**
 * @brief Converts a C++ double into an R real scalar.
 * @details Uses R's native `Rf_ScalarReal` for optimized, single-step
 * allocation and initialization of a length-1 REALSXP vector.
 * @param guard The active RProtectGuard for GC memory safety.
 * @param value The C++ double value to convert.
 * @return SEXP The resulting R real scalar.
 */
[[nodiscard]] inline SEXP FromDouble(RProtectGuard& guard, double value) {
    return guard.Protect(Rf_ScalarReal(value));
}

/**
 * @brief Converts a C++ std::string into an R character scalar.
 * @details Uses R's native `Rf_mkString` function to optimally allocate
 * and initialize a length-1 STRSXP (character vector) in a single step,
 * completely avoiding intermediate manual allocations or SET_STRING_ELT overhead.
 * @param guard The active RProtectGuard for GC memory safety.
 * @param value The C++ std::string value to convert.
 * @return SEXP The resulting R character scalar.
 */
[[nodiscard]] inline SEXP FromString(RProtectGuard& guard, const std::string& value) {
    // Rf_mkString alloca uno STRSXP di lunghezza 1 e lo popola
    return guard.Protect(Rf_mkString(value.c_str()));
}

/**
 * @brief Converts a C++ list of string pairs into an R character matrix.
 * @details Allocates an (N x 2) R character matrix (STRSXP) and populates it
 * utilizing R's native column-major layout. It explicitly enforces UTF-8
 * encoding to guarantee safe cross-platform text rendering in R.
 * @param guard The active RProtectGuard for GC memory safety.
 * @param pairs The C++ list of string pairs to convert.
 * @return SEXP The resulting N x 2 R character matrix.
 */
[[nodiscard]] inline SEXP FromStringPairList(
                                             RProtectGuard& guard,
                                             const std::list<std::pair<std::string, std::string>>& pairs) {
    
    // Usiamo R_len_t, il tipo ufficiale di R per dimensioni e indici
    const R_len_t n = static_cast<R_len_t>(pairs.size());
    
    // Alloca una matrice di stringhe N x 2
    SEXP result = guard.Protect(Rf_allocMatrix(STRSXP, n, 2));
    
    R_len_t pos = 0;
    for (const auto& p : pairs) {
        // Rf_mkCharCE con CE_UTF8 protegge da corruzioni di encoding (es. accenti)
        SET_STRING_ELT(result, pos,     Rf_mkCharCE(p.first.c_str(), CE_UTF8));
        SET_STRING_ELT(result, pos + n, Rf_mkCharCE(p.second.c_str(), CE_UTF8));
        ++pos;
    }
    
    return result;
}

/**
 * @brief Converts a C++ Tensor<double, 2> into an R matrix (REALSXP).
 *
 * @details This function is heavily optimized for High-Performance Computing (HPC).
 * It leverages Column-Major memory access (standard in R and Fortran) to
 * maximize CPU cache spatial locality. It also extracts raw pointers outside
 * of loops to prevent macro expansion overhead.
 *
 * @param guard A reference to the RProtectGuard to handle R's Garbage Collector.
 * @param values The C++ matrix (Tensor<double, 2>) passed by const reference to avoid copies.
 * @return SEXP The allocated and populated R matrix.
 */
inline SEXP FromDoubleMatrix(RProtectGuard& guard, const Tensor<double, 2>& values) {
    const int nrow = static_cast<int>(values.Extent(0));
    const int ncol = static_cast<int>(values.Extent(1));
    
    SEXP matrix = guard.Protect(Rf_allocMatrix(REALSXP, nrow, ncol));
    
    double* mat_ptr = REAL(matrix);
    
    for (int col = 0; col < ncol; ++col) {
        for (int row = 0; row < nrow; ++row) {
            mat_ptr[row + nrow * col] =
            values(static_cast<std::uint64_t>(row),
                   static_cast<std::uint64_t>(col));
        }
    }
    
    return matrix;
}

/**
 * @brief Converts a 3D C++ vector of doubles into an R 3D array.
 *
 * @details Safely extracts dimensions and allocates a 3D R array (REALSXP).
 * It uses a highly optimized loop structure that respects R's column-major
 * memory layout (r changes fastest, then c, then d). This ensures sequential
 * contiguous memory writes to the R array, preventing cache misses and
 * avoiding redundant index math.
 *
 * @param guard A reference to the RProtectGuard to manage R's Garbage Collector.
 * @param values The 3D nested C++ std::vector to convert.
 * @return SEXP The resulting 3D R array.
 */
inline SEXP From3DDoubleVector(RProtectGuard& guard,
                               const std::vector<std::vector<std::vector<double>>>& values) {
    // 1. Safety Checks per evitare Segfault su vettori vuoti
    int nrow = static_cast<int>(values.size());
    if (nrow == 0) {
        return guard.Protect(Rf_alloc3DArray(REALSXP, 0, 0, 0));
    }
    
    int ncol = static_cast<int>(values[0].size());
    if (ncol == 0) {
        return guard.Protect(Rf_alloc3DArray(REALSXP, nrow, 0, 0));
    }
    
    int ndepth = static_cast<int>(values[0][0].size());
    
    // 2. Allocazione
    SEXP result = guard.Protect(Rf_alloc3DArray(REALSXP, nrow, ncol, ndepth));
    
    // 3. Estrazione Raw Pointer per le massime prestazioni
    double* res_ptr = REAL(result);
    
    // 4. Cicli ottimizzati: L'ordine D -> C -> R rispecchia esattamente
    // la disposizione in memoria di R (Column-Major).
    std::size_t idx = 0;
    for (int d = 0; d < ndepth; ++d) {
        for (int c = 0; c < ncol; ++c) {
            for (int r = 0; r < nrow; ++r) {
                // Poiché scriviamo sequenzialmente, possiamo usare un semplice idx++
                // al posto di ricalcolare (r + c * nrow + ncol * nrow * d)
                res_ptr[idx++] = values[static_cast<std::size_t>(r)]
                [static_cast<std::size_t>(c)]
                [static_cast<std::size_t>(d)];
            }
        }
    }
    
    return result;
}


/**
 * @brief Attaches row and column names to an R matrix.
 *
 * @details In R, the 'dimnames' attribute of a matrix must be a list (VECSXP)
 * of length 2. The first element contains the row names (STRSXP) and the
 * second element contains the column names (STRSXP). This inline function
 * safely builds this list and attaches it to the target matrix.
 *
 * @param guard A reference to the RProtectGuard to manage R's Garbage Collector.
 * @param matrix_r The target R matrix (SEXP) to which names will be attached.
 * @param row_names An R character vector (STRSXP) containing the row names.
 * @param col_names An R character vector (STRSXP) containing the column names.
 */
inline void AttachDimNames(RProtectGuard& guard, SEXP matrix_r, SEXP row_names, SEXP col_names) {
    // Allocate a generic vector (list) of length 2 to hold the names.
    SEXP dimnames = guard.Protect(Rf_allocVector(VECSXP, 2));
    
    // Assign row names to the first slot (index 0)
    SET_VECTOR_ELT(dimnames, 0, row_names);
    
    // Assign column names to the second slot (index 1)
    SET_VECTOR_ELT(dimnames, 1, col_names);
    
    // Attach the newly created list as the 'dimnames' attribute of the matrix
    Rf_setAttrib(matrix_r, R_DimNamesSymbol, dimnames);
}

/**
 * @brief Safely assigns a value and its corresponding name to an R list.
 *
 * @details This utility function encapsulates the repetitive two-step process
 * of populating a named list in the R C API. It ensures that the value
 * and its name are always assigned to the exact same index, preventing
 * common copy-paste alignment bugs.
 *
 * @param list The target R list (VECSXP) being populated.
 * @param names The R character vector (STRSXP) storing the names of the list elements.
 * @param index The 0-based index where the element should be inserted.
 * @param value The R object (SEXP) to insert into the list.
 * @param name The C-string representing the name to assign to this element.
 */
inline void SetListElement(SEXP list, SEXP names, int index, SEXP value, const char* name) {
    SET_VECTOR_ELT(list, index, value);
    SET_STRING_ELT(names, index, Rf_mkChar(name));
}

/**
 * @brief Creates an empty R named list and its parallel names vector.
 *
 * Typical usage:
 * @code
 *   auto [list, names] = RCreate::NamedList(guard, 3);
 *   RCreate::SetListElement(list, names, 0, some_sexp, "key");
 *   Rf_setAttrib(list, R_NamesSymbol, names);
 * @endcode
 *
 * @param guard  Active guard.
 * @param size   Number of elements.
 * @note Call Rf_setAttrib(list, R_NamesSymbol, names) after populating.
 * @return `{VECSXP, STRSXP}` pair ready to be populated.
 */
inline std::pair<SEXP, SEXP> NamedList(RProtectGuard& guard, int size) {
    SEXP list = guard.Protect(Rf_allocVector(VECSXP, size));
    SEXP names = guard.Protect(Rf_allocVector(STRSXP, size));
    return {list, names};
}

/**
 * @brief Converte un seme a 64 bit nella sua rappresentazione decimale per R.
 * @details I semi non entrano negli interi di R (32 bit) e superano la
 * precisione esatta dei double (2^53): la stringa e' l'unico formato che li
 * restituisce senza perdita e che puo' essere rimesso in ingresso tale e quale.
 */
[[nodiscard]] inline SEXP FromSeed(RProtectGuard& guard, std::uint64_t seed) {
    return FromString(guard, std::to_string(seed));
}

/**
 * @brief Converte un vettore di semi a 64 bit in un vettore character di R.
 */
[[nodiscard]] inline SEXP FromSeeds(RProtectGuard& guard, const std::vector<std::uint64_t>& seeds) {
    SEXP result_r = guard.Protect(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(seeds.size())));
    for (std::size_t i = 0; i < seeds.size(); ++i) {
        const std::string text = std::to_string(seeds[i]);
        SET_STRING_ELT(result_r, static_cast<R_xlen_t>(i),
                       Rf_mkCharLenCE(text.data(), static_cast<int>(text.size()), CE_UTF8));
    }
    return result_r;
}

/**
 * @brief Impacchetta un generatore e il seme effettivamente usato per R.
 * @details Restituisce list(ptr = <externalptr>, seed = "<cifre>"). Il seme
 * viene sempre riportato, sia quando l'utente lo ha fornito sia quando e' stato
 * estratto a caso: reimmetterlo riproduce esattamente la stessa sessione.
 * @param ptr_r External pointer del generatore (gia' protetto se necessario).
 * @param seed Seme effettivamente usato.
 */
[[nodiscard]] inline SEXP GeneratorWithSeed(RProtectGuard& guard, SEXP ptr_r, std::uint64_t seed) {
    auto [list_r, names_r] = NamedList(guard, 2);
    SetListElement(list_r, names_r, 0, ptr_r, "ptr");
    SetListElement(list_r, names_r, 1, FromSeed(guard, seed), "seed");
    Rf_setAttrib(list_r, R_NamesSymbol, names_r);
    return list_r;
}

}
