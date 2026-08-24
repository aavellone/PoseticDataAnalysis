#' @title
#' Computing transitive closure.
#'
#' @description
#' Computes the transitive closure of the input binary relation.
#'
#' @param rel A two-column character matrix, each row comprising an element (pair) of the binary relation.
#'
#' @return
#' A transitive binary relation, as a two-column character matrix (each row comprises an element (pair)
#' of the transitive closure of the input relation).
#'
#' @examples
#' rel <- matrix(c(
#'   "a", "b",
#'   "c", "b",
#'   "d", "a",
#'   "c", "a",
#'   "a", "a",
#'   "b", "b",
#'   "c", "c",
#'   "d", "d"
#' ), ncol = 2, byrow = TRUE)
#'
#' t.clo<-TransitiveClosure(rel)
#'
#' @name TransitiveClosure
#' @export
TransitiveClosure <- function(rel) {
  if (!is.matrix(rel) || ncol(rel) != 2) {
    stop("'rel' must be a two-column matrix.", call. = FALSE)
  }
  if (nrow(rel) > 0 && !is.character(rel)) {
    stop("'rel' must be a two-column character matrix.", call. = FALSE)
  }

  tryCatch({
    return(.Call(C_TransitiveClosure, rel))
  }, error = function(err) {
    err_msg   <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    clean_msg <- if (length(err_parts) > 1L) trimws(err_parts[length(err_parts)]) else trimws(err_msg)
    stop(clean_msg, call. = FALSE)
  })
}
