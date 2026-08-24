#' @title
#' Constructing a Partially Ordered Set from a Cover Matrix.
#'
#' @description
#' Constructs an object of class `POSet` starting from a cover matrix
#' (which represents the covering relations, or the edges of the Hasse diagram).
#' Element labels are extracted from the row/column names of the matrix.
#' If no names are provided, they are automatically generated.
#'
#' @param cover_matrix A square numeric (0/1) or logical matrix representing
#' the covering relation. An entry `[i, j] == 1` (or `TRUE`) implies that
#' the element `i` is covered by element `j` (i.e., $i < j$ and there is no
#' element $k$ such that $i < k < j$).
#'
#' @return
#' An object of S4 class `POSet`.
#'
#' @examples
#' # Example with an unnamed logical cover matrix (Hasse diagram edges)
#' cov_mat <- matrix(c(FALSE, TRUE,  FALSE,
#'                     FALSE, FALSE, TRUE,
#'                     FALSE, FALSE, FALSE),
#'                   nrow = 3, byrow = TRUE)
#' pos <- POSetFromCover(cov_mat)
#'
#' # Example with a named numeric cover matrix
#' cov_mat_named <- matrix(c(0, 1, 0,
#'                           0, 0, 1,
#'                           0, 0, 0), nrow = 3, byrow = TRUE)
#' rownames(cov_mat_named) <- colnames(cov_mat_named) <- c("A", "B", "C")
#' pos_named <- POSetFromCover(cov_mat_named)
#'
#' @name POSetFromCover
#' @export
POSetFromCover <- function(cover_matrix) {

  # --- 1. Validazione dell'input ---
  # Controlla che l'input sia effettivamente una matrice
  if (!is.matrix(cover_matrix)) {
    stop("'cover_matrix' must be a matrix.", call. = FALSE)
  }

  # La matrice di copertura tra gli elementi deve essere quadrata
  if (nrow(cover_matrix) != ncol(cover_matrix)) {
    stop("'cover_matrix' must be a square matrix.", call. = FALSE)
  }

  # Verifica che la matrice contenga solo valori validi (0/1 oppure logici TRUE/FALSE)
  if (!is.numeric(cover_matrix) && !is.logical(cover_matrix)) {
    stop("'cover_matrix' must be either numeric (0/1) or logical.", call. = FALSE)
  }

  # --- Controllo Antisimmetria Diretto (Cicli di lunghezza 2) ---
  M_check <- cover_matrix
  diag(M_check) <- 0

  if (any(M_check & t(M_check))) {
    stop("The cover matrix violates the antisymmetry property (direct loop detected).", call. = FALSE)
  }


  # --- 2. Gestione delle etichette degli elementi (Ground Set) ---
  n_elements <- nrow(cover_matrix)

  # Tenta di recuperare i nomi dalle righe; se assenti, prova con le colonne
  labels <- rownames(cover_matrix)
  if (is.null(labels)) {
    labels <- colnames(cover_matrix)
  }

  # Se la matrice non ha alcun nome, genera etichette automatiche (es. "e1", "e2", ...)
  if (is.null(labels)) {
    labels <- paste0("e", seq_len(n_elements))
  } else if (any(labels == "") || anyNA(labels)) {
    stop("Row/column names of the cover matrix cannot be empty or NA.", call. = FALSE)
  }

  # --- 3. Estrazione delle relazioni di copertura ---
  # Crea due matrici che contengono rispettivamente gli indici di riga e di colonna
  r <- row(cover_matrix)
  c <- col(cover_matrix)

  # Trova le coordinate (i, j) dove c'è un arco di copertura (valore > 0 o TRUE).
  # Nei diagrammi di Hasse (cover matrix) la diagonale è già zero per definizione,
  # ma inseriamo 'r != c' come misura di sicurezza aggiuntiva per garantire l'antiriflessività.
  idx <- which(cover_matrix > 0 & r != c)

  # Se ci sono relazioni (archi nel diagramma di Hasse), costruisce la matrice a due colonne 'dom'
  if (length(idx) > 0) {
    # Mappa gli indici numerici estratti sulle etichette corrispondenti
    from_elems <- labels[r[idx]]
    to_elems   <- labels[c[idx]]
    dom <- cbind(from_elems, to_elems)
  } else {
    # Se il poset è un'anticatena (nessun collegamento, matrice vuota di coperture),
    # crea una matrice vuota col formato richiesto dal costruttore POSet
    dom <- matrix(character(0), ncol = 2, nrow = 0)
  }

  # --- 4. Costruzione finale ---
  # Passa gli elementi e le relazioni base alla funzione costruttrice POSet,
  # la quale ricostruirà l'ordine parziale completo a partire da questi collegamenti diretti.
  return(poseticDataAnalysis::POSet(elements = labels, dom = dom))
}
