/**
 * @file mainTestR.cpp
 * @brief Eseguibile standalone (R embedded) per testare/debuggare la funzione
 * C `FirstOrderDominanceAnalysis` senza passare dal package R.
 *
 * Test corrente: valutazione FOD multi-catena su un poset "corposo".
 *  - Poset: 40 elementi, DAG casuale riproducibile (set.seed in R).
 *  - Metriche: Dominance, MannWhitneyDominance e MRP (quest'ultima produce
 *    la matrice elemento x elemento del poset, solo binMatrix).
 *  - Run A (riferimento): percorso classico, catena Bubley-Dyer singola,
 *    4 * N estensioni lineari.
 *  - Run B: percorso multi-catena, 4 catene * N estensioni (stesso lavoro
 *    totale), punti di partenza diversi (LE generate in R con Kahn
 *    randomizzato) e un seed per catena.
 *  Confronto: tempi wall e max|diff| tra le matrici risultato (MRP inclusa
 *  nel check di determinismo B vs C).
 *
 * COMPILAZIONE (Linux / macOS):
 *   export R_HOME=$(R RHOME)
 *   g++ -std=c++20 -I"$R_HOME/include" mainTestR.cpp <oggetti .o del package> \
 *       -L"$R_HOME/lib" -lR -o fod_main
 *
 * ESECUZIONE:
 *   R_HOME deve essere impostata nell'ambiente (export R_HOME=$(R RHOME)),
 *   altrimenti Rf_initEmbeddedR fallisce.
 */

// Gli header C++ standard vanno inclusi PRIMA di quelli di R: le macro di R
// (length, error, ...) altrimenti corrompono gli header libc++.
#include <chrono>
#include <cstdio>

#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <R_ext/Parse.h>

// ---------------------------------------------------------------------------
// Funzioni C del package (linkate insieme a questo main)
// ---------------------------------------------------------------------------
extern "C" {
SEXP BuildPOSet(SEXP elements, SEXP comparabilities);
SEXP FirstOrderDominanceAnalysis(SEXP poset_r, SEXP freq_matrix_r,
                                 SEXP subpopulation_count_r, SEXP metrics_r,
                                 SEXP total_bins_r, SEXP count_r,
                                 SEXP seed_r, SEXP output_interval_in_sec_r,
                                 SEXP sep_r, SEXP linear_extensions_r,
                                 SEXP n_threads_r);
}

namespace {

/**
 * @brief Esegue una stringa di codice R (parse + eval nel GlobalEnv).
 * @return Il valore dell'ultima espressione valutata (NON protetto:
 * proteggere subito con PROTECT se va conservato).
 * @note In caso di errore di parsing o valutazione chiama Rf_error.
 */
SEXP EvalRString(const char* code) {
    ParseStatus status = PARSE_NULL;
    SEXP code_r = PROTECT(Rf_mkString(code));
    SEXP exprs = PROTECT(R_ParseVector(code_r, -1, &status, R_NilValue));
    if (status != PARSE_OK) {
        UNPROTECT(2);
        Rf_error("EvalRString: parse error in \"%s\".", code);
    }

    SEXP result = R_NilValue;
    for (R_xlen_t i = 0; i < Rf_xlength(exprs); ++i) {
        int error_occurred = 0;
        result = R_tryEval(VECTOR_ELT(exprs, i), R_GlobalEnv, &error_occurred);
        if (error_occurred) {
            UNPROTECT(2);
            Rf_error("EvalRString: evaluation error in \"%s\".", code);
        }
    }

    UNPROTECT(2);
    return result;
}

// ---------------------------------------------------------------------------
// Wrapper per R_ToplevelExec.
//
// Le funzioni del package segnalano gli errori con Rf_error (longjmp): se
// chiamate direttamente da main, senza un contesto R attivo, il longjmp
// abbatterebbe il processo. R_ToplevelExec stabilisce un contesto top-level:
// in caso di errore ritorna FALSE invece di terminare.
// ---------------------------------------------------------------------------

struct BuildPOSetCall {
    SEXP elements = nullptr;
    SEXP comparabilities = nullptr;
    SEXP result = nullptr;
};

void RunBuildPOSet(void* data) {
    auto* call = static_cast<BuildPOSetCall*>(data);
    call->result = BuildPOSet(call->elements, call->comparabilities);
}

struct FODCall {
    SEXP args[11] = {nullptr};
    SEXP result = nullptr;
};

void RunFOD(void* data) {
    auto* call = static_cast<FODCall*>(data);
    call->result = FirstOrderDominanceAnalysis(
                                               call->args[0], call->args[1], call->args[2],
                                               call->args[3], call->args[4], call->args[5],
                                               call->args[6], call->args[7], call->args[8],
                                               call->args[9], call->args[10]);
}

}  // namespace

int main(int argc, char* argv[]) {
    (void)argc;
    (void)argv;

    // -----------------------------------------------------------------------
    // 0. Inizializzazione dell'interprete R embedded (richiede R_HOME)
    // -----------------------------------------------------------------------
    char* r_argv[] = {
        const_cast<char*>("R"),
        const_cast<char*>("--vanilla"),
        const_cast<char*>("--silent")
    };
    Rf_initEmbeddedR(sizeof(r_argv) / sizeof(r_argv[0]), r_argv);

    int n_protected = 0;

    // -----------------------------------------------------------------------
    // 1. Poset "corposo": 40 elementi, DAG casuale riproducibile.
    //    Le relazioni vi < vj sono generate solo per i < j: aciclicita'
    //    garantita per costruzione. La chiusura transitiva la fa BuildPOSet.
    // -----------------------------------------------------------------------
    EvalRString(
        "set.seed(123)\n"
        "n_el <- 40L\n"
        "elems <- paste0('v', seq_len(n_el))\n"
        "edge_list <- list()\n"
        "for (i in 1:(n_el - 1)) for (j in (i + 1):n_el) {\n"
        "  if (runif(1) < 0.10) edge_list[[length(edge_list) + 1L]] <- c(elems[i], elems[j])\n"
        "}\n"
        "edges <- do.call(rbind, edge_list)\n");

    SEXP elements = PROTECT(EvalRString("elems"));
    ++n_protected;
    SEXP edges = PROTECT(EvalRString("edges"));
    ++n_protected;

    std::printf("Poset: 40 elementi, %d relazioni di copertura generate.\n",
                Rf_nrows(edges));

    BuildPOSetCall build_call;
    build_call.elements = elements;
    build_call.comparabilities = edges;

    if (R_ToplevelExec(RunBuildPOSet, &build_call) == FALSE) {
        REprintf("Errore R durante BuildPOSet: impossibile costruire il POSet.\n");
        UNPROTECT(n_protected);
        Rf_endEmbeddedR(0);
        return 1;
    }
    SEXP poset_ptr = PROTECT(build_call.result);
    ++n_protected;

    // -----------------------------------------------------------------------
    // 2. Matrice delle frequenze: tutti i 40 elementi come valori osservati,
    //    3 gruppi, colonne normalizzate a somma 1.
    // -----------------------------------------------------------------------
    SEXP freq_matrix = PROTECT(EvalRString(
        "freq <- matrix(runif(n_el * 3), nrow = n_el,\n"
        "               dimnames = list(elems, c('GruppoA', 'GruppoB', 'GruppoC')))\n"
        "freq <- sweep(freq, 2, colSums(freq), '/')\n"
        "freq"));
    ++n_protected;

    // -----------------------------------------------------------------------
    // 3. Parametri comuni della chiamata
    // -----------------------------------------------------------------------
    // metrics: (Dominance, MannWhitney, Inferential, MRP). La MRP produce la
    // matrice elemento x elemento del poset: nel risultato ha solo binMatrix
    // (fodClosed/fodMatrix sono NULL, nessuna analisi FOD gruppo x gruppo).
    SEXP metrics = PROTECT(Rf_allocVector(LGLSXP, 4));
    ++n_protected;
    LOGICAL(metrics)[0] = TRUE;
    LOGICAL(metrics)[1] = TRUE;
    LOGICAL(metrics)[2] = FALSE;
    LOGICAL(metrics)[3] = TRUE;

    SEXP total_bins = PROTECT(Rf_ScalarInteger(10));
    ++n_protected;

    SEXP sep = PROTECT(Rf_mkString("_"));
    ++n_protected;

    SEXP subpopulation_count = R_NilValue;
    SEXP output_interval = PROTECT(Rf_ScalarInteger(5));
    ++n_protected;

    constexpr int kNChains = 4;
    constexpr int kLePerChain = 100000;

    // -----------------------------------------------------------------------
    // 4. RUN A (riferimento): catena singola classica, 4 * N estensioni.
    // -----------------------------------------------------------------------
    SEXP count_seq = PROTECT(Rf_ScalarInteger(kNChains * kLePerChain));
    ++n_protected;
    SEXP seed_seq = PROTECT(Rf_ScalarInteger(42));
    ++n_protected;

    FODCall fod_seq;
    fod_seq.args[0] = poset_ptr;
    fod_seq.args[1] = freq_matrix;
    fod_seq.args[2] = subpopulation_count;
    fod_seq.args[3] = metrics;
    fod_seq.args[4] = total_bins;
    fod_seq.args[5] = count_seq;
    fod_seq.args[6] = seed_seq;
    fod_seq.args[7] = output_interval;
    fod_seq.args[8] = sep;
    fod_seq.args[9] = R_NilValue;   // linear_extensions: NULL -> percorso classico
    fod_seq.args[10] = R_NilValue;  // n_threads: non usato nel percorso classico

    std::printf("\n--- RUN A: catena singola, %d LE ---\n", kNChains * kLePerChain);
    const auto t0 = std::chrono::steady_clock::now();
    if (R_ToplevelExec(RunFOD, &fod_seq) == FALSE) {
        REprintf("Errore R durante FirstOrderDominanceAnalysis (catena singola).\n");
        UNPROTECT(n_protected);
        Rf_endEmbeddedR(0);
        return 1;
    }
    const auto t1 = std::chrono::steady_clock::now();
    SEXP result_seq = PROTECT(fod_seq.result);
    ++n_protected;

    // -----------------------------------------------------------------------
    // 5. Matrice delle estensioni lineari di partenza: 4 colonne, generate in
    //    R con un ordinamento topologico randomizzato (Kahn). Ogni colonna e'
    //    una LE valida del poset per costruzione.
    // -----------------------------------------------------------------------
    SEXP le_matrix = PROTECT(EvalRString(
        "rand_le <- function(elems, edges) {\n"
        "  adj <- setNames(vector('list', length(elems)), elems)\n"
        "  indeg <- setNames(integer(length(elems)), elems)\n"
        "  for (r in seq_len(nrow(edges))) {\n"
        "    a <- edges[r, 1]; b <- edges[r, 2]\n"
        "    adj[[a]] <- c(adj[[a]], b)\n"
        "    indeg[b] <- indeg[b] + 1L\n"
        "  }\n"
        "  out <- character(0)\n"
        "  avail <- names(indeg)[indeg == 0L]\n"
        "  while (length(avail) > 0) {\n"
        "    x <- if (length(avail) == 1L) avail else sample(avail, 1)\n"
        "    out <- c(out, x)\n"
        "    avail <- setdiff(avail, x)\n"
        "    for (y in adj[[x]]) {\n"
        "      indeg[y] <- indeg[y] - 1L\n"
        "      if (indeg[y] == 0L) avail <- c(avail, y)\n"
        "    }\n"
        "  }\n"
        "  out\n"
        "}\n"
        "set.seed(456)\n"
        "le_matrix <- sapply(1:4, function(k) rand_le(elems, edges))\n"
        "le_matrix"));
    ++n_protected;

    SEXP chain_seeds = PROTECT(EvalRString("c(101L, 202L, 303L, 404L)"));
    ++n_protected;
    SEXP count_chain = PROTECT(Rf_ScalarInteger(kLePerChain));
    ++n_protected;

    // -----------------------------------------------------------------------
    // 6. RUN B: multi-catena, 4 catene * N estensioni (stesso lavoro totale),
    //    n_threads automatico (NULL -> hardware_concurrency worker; quelli in
    //    eccesso rispetto alle catene trovano la sorgente esaurita ed escono).
    // -----------------------------------------------------------------------
    FODCall fod_mc = fod_seq;
    fod_mc.args[5] = count_chain;   // N per catena
    fod_mc.args[6] = chain_seeds;   // un seed per catena
    fod_mc.args[9] = le_matrix;     // una colonna = una LE di partenza
    fod_mc.args[10] = R_NilValue;   // n_threads automatico
    fod_mc.result = nullptr;

    std::printf("\n--- RUN B: %d catene x %d LE, n_threads = auto ---\n",
                kNChains, kLePerChain);
    const auto t2 = std::chrono::steady_clock::now();
    if (R_ToplevelExec(RunFOD, &fod_mc) == FALSE) {
        REprintf("Errore R durante FirstOrderDominanceAnalysis (multi-catena).\n");
        UNPROTECT(n_protected);
        Rf_endEmbeddedR(0);
        return 1;
    }
    const auto t3 = std::chrono::steady_clock::now();
    SEXP result_mc = PROTECT(fod_mc.result);
    ++n_protected;

    // -----------------------------------------------------------------------
    // 7. RUN C: stesse 4 catene e stessi seed, ma n_threads = 2.
    //    I risultati devono essere IDENTICI a RUN B (dipendono solo da
    //    catene e seed, non dal numero di worker); cambia solo il tempo.
    // -----------------------------------------------------------------------
    SEXP n_threads_2 = PROTECT(Rf_ScalarInteger(2));
    ++n_protected;

    FODCall fod_mc2 = fod_mc;
    fod_mc2.args[10] = n_threads_2;
    fod_mc2.result = nullptr;

    std::printf("\n--- RUN C: %d catene x %d LE, n_threads = 2 ---\n",
                kNChains, kLePerChain);
    const auto t4 = std::chrono::steady_clock::now();
    if (R_ToplevelExec(RunFOD, &fod_mc2) == FALSE) {
        REprintf("Errore R durante FirstOrderDominanceAnalysis (multi-catena, 2 thread).\n");
        UNPROTECT(n_protected);
        Rf_endEmbeddedR(0);
        return 1;
    }
    const auto t5 = std::chrono::steady_clock::now();
    SEXP result_mc2 = PROTECT(fod_mc2.result);
    ++n_protected;

    // -----------------------------------------------------------------------
    // 8. Confronto risultati e tempi
    // -----------------------------------------------------------------------
    Rf_defineVar(Rf_install("ris_seq"), result_seq, R_GlobalEnv);
    Rf_defineVar(Rf_install("ris_mc"), result_mc, R_GlobalEnv);
    Rf_defineVar(Rf_install("ris_mc2"), result_mc2, R_GlobalEnv);

    EvalRString(
        "cat('\\nRUN A: LEType =', ris_seq$LEType,"
        "    '- le_count =', attr(ris_seq, 'le_count'), '\\n')\n"
        "cat('RUN B: LEType =', ris_mc$LEType,"
        "    '- le_count =', attr(ris_mc, 'le_count'), '\\n')\n"
        "cat('RUN C: LEType =', ris_mc2$LEType,"
        "    '- le_count =', attr(ris_mc2, 'le_count'), '\\n')\n"
        "for (m in c('Dominance', 'MannWhitneyDominance')) {\n"
        "  cat(sprintf('%-22s max|binMatrix A-B| = %.6f  max|fodClosed A-B| = %.6f\\n', m,\n"
        "      max(abs(ris_mc[[m]]$binMatrix - ris_seq[[m]]$binMatrix)),\n"
        "      max(abs(ris_mc[[m]]$fodClosed - ris_seq[[m]]$fodClosed))))\n"
        "}\n"
        "stopifnot(is.null(ris_mc$MRP$fodClosed), is.matrix(ris_mc$MRP$binMatrix))\n"
        "cat(sprintf('%-22s max|binMatrix A-B| = %.6f  (fodClosed = NULL, atteso)\\n', 'MRP',\n"
        "    max(abs(ris_mc$MRP$binMatrix - ris_seq$MRP$binMatrix))))\n"
        "diff_bc <- max(sapply(c('Dominance', 'MannWhitneyDominance', 'MRP'), function(m)\n"
        "  max(abs(ris_mc[[m]]$binMatrix - ris_mc2[[m]]$binMatrix))))\n"
        "cat(sprintf('determinismo B vs C (atteso 0, MRP inclusa): max|diff| = %g -> %s\\n',\n"
        "    diff_bc, if (diff_bc == 0) 'OK' else 'FAIL'))");

    const double sec_seq = std::chrono::duration<double>(t1 - t0).count();
    const double sec_mc  = std::chrono::duration<double>(t3 - t2).count();
    const double sec_mc2 = std::chrono::duration<double>(t5 - t4).count();
    std::printf("\nTempo RUN A (1 catena):            %.3f s\n", sec_seq);
    std::printf("Tempo RUN B (%d catene, auto):      %.3f s  -> speedup %.2fx\n",
                kNChains, sec_mc, sec_seq / sec_mc);
    std::printf("Tempo RUN C (%d catene, 2 thread):  %.3f s  -> speedup %.2fx\n",
                kNChains, sec_mc2, sec_seq / sec_mc2);

    // -----------------------------------------------------------------------
    // 9. Chiusura pulita
    // -----------------------------------------------------------------------
    UNPROTECT(n_protected);
    Rf_endEmbeddedR(0);
    return 0;
}
