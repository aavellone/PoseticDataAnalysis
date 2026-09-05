library(poseticDataAnalysis)

# -----------------------------------------------------------------------------
# 1. Poset corposo: 40 elementi, relazioni vi < vj generate solo per i < j
#    (aciclicita' garantita); la chiusura transitiva la fa il backend.
# -----------------------------------------------------------------------------
set.seed(123)
n_el <- 40L
elems <- paste0("v", seq_len(n_el))

edge_list <- list()
for (i in 1:(n_el - 1)) for (j in (i + 1):n_el) {
  if (runif(1) < 0.10) edge_list[[length(edge_list) + 1L]] <- c(elems[i], elems[j])
}
edges <- do.call(rbind, edge_list)
cat(sprintf("Poset: %d elementi, %d relazioni di copertura generate.\n",
            n_el, nrow(edges)))

poset <- POSet(elements = elems, dom = edges)

# -----------------------------------------------------------------------------
# 2. Matrice delle frequenze: tutti gli elementi come valori osservati,
#    3 gruppi, colonne normalizzate a somma 1.
# -----------------------------------------------------------------------------
freq <- matrix(runif(n_el * 3), nrow = n_el,
               dimnames = list(elems, c("GruppoA", "GruppoB", "GruppoC")))
freq <- sweep(freq, 2, colSums(freq), "/")

# -----------------------------------------------------------------------------
# 3. Parametri comuni
# -----------------------------------------------------------------------------
metrics      <- c("Dominance", "MannWhitneyDominance", "MRP")
total_bins   <- 10L
n_chains     <- 4L
le_per_chain <- 100000L
out_interval <- 5L

# -----------------------------------------------------------------------------
# 4. RUN A (riferimento): catena singola classica, 4 * N estensioni.
# -----------------------------------------------------------------------------
cat(sprintf("\n--- RUN A: catena singola, %d LE ---\n", n_chains * le_per_chain))
t_a <- system.time(
  ris_seq <- FirstOrderDominanceAnalysis(
    posets = poset,
    freq_matrix = freq,
    metrics = metrics,
    total_bins = total_bins,
    count = n_chains * le_per_chain,
    seed = 42L,
    output_interval_in_sec = out_interval
  )
)["elapsed"]

# -----------------------------------------------------------------------------
# 5. Matrice delle estensioni lineari di partenza: 4 colonne generate con un
#    ordinamento topologico randomizzato (Kahn). Ogni colonna e' una LE valida
#    del poset per costruzione.
# -----------------------------------------------------------------------------
rand_le <- function(elems, edges) {
  adj <- setNames(vector("list", length(elems)), elems)
  indeg <- setNames(integer(length(elems)), elems)
  for (r in seq_len(nrow(edges))) {
    a <- edges[r, 1]; b <- edges[r, 2]
    adj[[a]] <- c(adj[[a]], b)
    indeg[b] <- indeg[b] + 1L
  }
  out <- character(0)
  avail <- names(indeg)[indeg == 0L]
  while (length(avail) > 0) {
    x <- if (length(avail) == 1L) avail else sample(avail, 1)
    out <- c(out, x)
    avail <- setdiff(avail, x)
    for (y in adj[[x]]) {
      indeg[y] <- indeg[y] - 1L
      if (indeg[y] == 0L) avail <- c(avail, y)
    }
  }
  out
}

set.seed(456)
le_matrix <- sapply(seq_len(n_chains), function(k) rand_le(elems, edges))
chain_seeds <- c(101L, 202L, 303L, 404L)

# -----------------------------------------------------------------------------
# 6. RUN B: multi-catena, n_threads automatico.
# -----------------------------------------------------------------------------
cat(sprintf("\n--- RUN B: %d catene x %d LE, n_threads = auto ---\n",
            n_chains, le_per_chain))
t_b <- system.time(
  ris_mc <- FirstOrderDominanceAnalysis(
    posets = poset,
    freq_matrix = freq,
    metrics = metrics,
    total_bins = total_bins,
    count = le_per_chain,               # N per catena
    seed = chain_seeds,                 # un seed per catena
    output_interval_in_sec = out_interval,
    linear_extensions = le_matrix       # una colonna = una LE di partenza
  )
)["elapsed"]

# -----------------------------------------------------------------------------
# 7. RUN C: stesse catene e stessi seed, n_threads = 2.
# -----------------------------------------------------------------------------
cat(sprintf("\n--- RUN C: %d catene x %d LE, n_threads = 2 ---\n",
            n_chains, le_per_chain))
t_c <- system.time(
  ris_mc2 <- FirstOrderDominanceAnalysis(
    posets = poset,
    freq_matrix = freq,
    metrics = metrics,
    total_bins = total_bins,
    count = le_per_chain,
    seed = chain_seeds,
    output_interval_in_sec = out_interval,
    linear_extensions = le_matrix,
    n_threads = 2L
  )
)["elapsed"]

# -----------------------------------------------------------------------------
# 8. Confronto risultati e tempi.
#    Nomi del post-processing R: bin = binMatrix C, mintr.delta = fodClosed C.
#    La MRP e' a parte: il suo risultato e' direttamente la matrice
#    elemento x elemento (nessun post-processing FOD).
# -----------------------------------------------------------------------------
cat(sprintf("\nRUN A: LEType = %s - le_count = %g\n",
            ris_seq$LEType, attr(ris_seq, "le_count")))
cat(sprintf("RUN B: LEType = %s - le_count = %g\n",
            ris_mc$LEType, attr(ris_mc, "le_count")))
cat(sprintf("RUN C: LEType = %s - le_count = %g\n",
            ris_mc2$LEType, attr(ris_mc2, "le_count")))

fod_metrics <- setdiff(metrics, "MRP")
for (m in fod_metrics) {
  cat(sprintf("%-22s max|bin A-B| = %.6f  max|mintr.delta A-B| = %.6f\n", m,
              max(abs(ris_mc[[m]]$bin - ris_seq[[m]]$bin)),
              max(abs(ris_mc[[m]]$mintr.delta - ris_seq[[m]]$mintr.delta))))
}
if ("MRP" %in% metrics) {
  cat(sprintf("%-22s max|matrice A-B| = %.6f\n", "MRP",
              max(abs(ris_mc$MRP - ris_seq$MRP))))
}

diff_bc <- max(sapply(fod_metrics, function(m)
  max(abs(ris_mc[[m]]$bin - ris_mc2[[m]]$bin))))
if ("MRP" %in% metrics) {
  diff_bc <- max(diff_bc, max(abs(ris_mc$MRP - ris_mc2$MRP)))
}
determinism_ok <- diff_bc == 0
cat(sprintf("determinismo B vs C (atteso 0, MRP inclusa): max|diff| = %g -> %s\n",
            diff_bc, if (determinism_ok) "OK" else "FAIL"))

cat(sprintf("\nTempo RUN A (1 catena):            %.3f s\n", t_a))
cat(sprintf("Tempo RUN B (%d catene, auto):      %.3f s  -> speedup %.2fx\n",
            n_chains, t_b, t_a / t_b))
cat(sprintf("Tempo RUN C (%d catene, 2 thread):  %.3f s  -> speedup %.2fx\n",
            n_chains, t_c, t_a / t_c))

#quit(status = if (determinism_ok) 0L else 1L)
