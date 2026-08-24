#' @title
#' Constructing a Partially Ordered Set from an Dominance Matrix.
#'
#' @description
#' Constructs an object of class `POSet` starting from an dominance matrix.
#' Element labels are extracted from the row/column names of the matrix.
#' If no names are provided, they are automatically generated.
#'
#' @param dominance_matrix A square numeric (0/1) or logical matrix representing
#' the partial order relation. An entry `[i, j] == 1` (or `TRUE`) implies that
#' the element `i` is less than or equal to element `j` (\eqn{i \leq j}).
#'
#' @return
#' An object of S4 class `POSet`.
#'
#' @examples
#' # Example with an unnamed logical matrix
#' inc_mat <- matrix(c(TRUE, FALSE, FALSE,
#'                     TRUE, TRUE,  FALSE,
#'                     TRUE, TRUE,  TRUE),
#'                   nrow = 3, byrow = TRUE)
#' pos <- POSetFromDominance(inc_mat)
#'
#' # Example with a named numeric matrix
#' inc_mat_named <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = TRUE)
#' rownames(inc_mat_named) <- colnames(inc_mat_named) <- c("A", "B")
#' pos_named <- POSetFromDominance(inc_mat_named)
#'
#' @name POSetFromDominance
#' @export
POSetFromDominance <- function(dominance_matrix) {

  # --- 1. Validazione dell'input ---
  # Controlla che l'input sia effettivamente una matrice
  if (!is.matrix(dominance_matrix)) {
    stop("'dominance_matrix' must be a matrix.", call. = FALSE)
  }

  # Una matrice di incidenza di un poset deve essere necessariamente quadrata
  if (nrow(dominance_matrix) != ncol(dominance_matrix)) {
    stop("'dominance_matrix' must be a square matrix.", call. = FALSE)
  }

  # Verifica che la matrice contenga solo valori validi (0/1 oppure logici TRUE/FALSE)
  if (!is.numeric(dominance_matrix) && !is.logical(dominance_matrix)) {
    stop("'dominance_matrix' must be either numeric (0/1) or logical.", call. = FALSE)
  }

  # --- Controllo Antisimmetria ---
  # Creiamo una copia temporanea azzerando la diagonale (perché A <= A è ovvio e simmetrico)
  M_check <- dominance_matrix
  diag(M_check) <- 0 # o FALSE se logica

  # L'operatore '&' tra la matrice e la sua trasposta trova le celle simmetriche.
  # Se 'any' è TRUE, significa che per qualche i != j hai sia i <= j che j <= i.
  if (any(M_check & t(M_check))) {
    stop("The dominance matrix violates the antisymmetry property (cycles detected).", call. = FALSE)
  }


  # --- 2. Gestione delle etichette degli elementi (Ground Set) ---
  n_elements <- nrow(dominance_matrix)

  # Tenta di recuperare i nomi dalle righe; se assenti, prova con le colonne
  labels <- rownames(dominance_matrix)
  if (is.null(labels)) {
    labels <- colnames(dominance_matrix)
  }

  # Se la matrice non ha alcun nome, genera etichette automatiche (es. "e1", "e2", ...)
  if (is.null(labels)) {
    labels <- paste0("e", seq_len(n_elements))
  } else if (any(labels == "") || anyNA(labels)) {
    stop("Row/column names of the dominance matrix cannot be empty or NA.", call. = FALSE)
  }

  # --- 3. Estrazione delle relazioni di dominanza ---
  # Crea due matrici che contengono rispettivamente gli indici di riga e di colonna
  r <- row(dominance_matrix)
  c <- col(dominance_matrix)

  # Trova le coordinate (i, j) dove c'è una relazione di dominanza (valore > 0 o TRUE).
  # Viene esclusa la diagonale principale (r != c) perché la riflessività
  # è una proprietà implicita dei poset e non è necessario dichiarare che a <= a.
  idx <- which(dominance_matrix > 0 & r != c)

  # Se ci sono relazioni di dominanza, costruisce la matrice a due colonne 'dom'
  if (length(idx) > 0) {
    # Mappa gli indici numerici estratti sulle etichette corrispondenti
    from_elems <- labels[r[idx]]
    to_elems   <- labels[c[idx]]
    dom <- cbind(from_elems, to_elems)
  } else {
    # Se il poset è un'anticatena (nessun elemento domina un altro),
    # crea una matrice vuota col formato richiesto dal costruttore POSet
    dom <- matrix(character(0), ncol = 2, nrow = 0)
  }

  # --- 4. Costruzione finale ---
  # Passa gli elementi e la matrice delle dominanze alla funzione costruttrice base
  return(poseticDataAnalysis::POSet(elements = labels, dom = dom))
}
