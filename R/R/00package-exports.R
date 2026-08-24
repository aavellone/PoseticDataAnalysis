#' @useDynLib poseticDataAnalysis, .registration=TRUE, .fixes="C_"
"_PACKAGE"

# ---------------------------------------------------------------------------
# Conversione dei semi per il backend C++ (funzione interna, non esportata)
#
# I semi del generatore C++ sono interi a 64 bit senza segno, un intervallo che
# R non sa rappresentare: 'integer' e' a 32 bit e 'double' e' esatto solo fino
# a 2^53. Il seme viaggia percio' verso C++ come STRINGA di cifre decimali,
# che attraversa il confine senza perdita di precisione.
#
# L'utente puo' indifferentemente passare un numero (comodo, entro i limiti di
# R) o una stringa di cifre (necessaria per usare tutto l'intervallo a 64 bit,
# per esempio per riprodurre un seme generato dal C++).
#
# Vettorizzata: serve anche ai semi per-catena di FirstOrderDominanceAnalysis.
# ---------------------------------------------------------------------------
.seed_to_character <- function(seed, name = "seed") {
  if (is.null(seed)) {
    return(NULL)
  }

  if (is.character(seed)) {
    if (anyNA(seed) || !all(grepl("^[0-9]+$", seed))) {
      stop(sprintf("'%s' must contain only decimal digits when given as a string.", name),
           call. = FALSE)
    }
    return(seed)
  }

  if (!is.numeric(seed) || anyNA(seed) || !all(is.finite(seed)) ||
      any(seed < 0) || any(seed != trunc(seed))) {
    stop(sprintf(paste0("'%s' must be a non-negative whole number, or a string of ",
                        "decimal digits for values above 2^53."), name),
         call. = FALSE)
  }

  # sprintf() invece di as.character(): quest'ultima usa la notazione
  # scientifica per i valori grandi (es. "1e+10"), che il parser C++ rifiuta.
  sprintf("%.0f", seed)
}
