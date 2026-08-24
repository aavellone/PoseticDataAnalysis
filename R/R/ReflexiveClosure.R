#' @title
#' Computing reflexive closure.
#'
#' @description
#' Computes the reflexive closure of the input binary relation.
#'
#' @param set A character vector (the names of the elements of the set on which the binary relation is defined).
#'
#' @param rel A two-column character matrix, each row comprising an element (pair) of the relation.
#'
#' @return
#' A reflexive binary relation, as a two-column character matrix
#' (each row comprises an element (pair) of the reflexive closure of the input relation).
#'
#' @examples
#' set<-c("a", "b", "c", "d", "e")
#'
#' rel <- matrix(c(
#'   "a", "b",
#'   "c", "b",
#'   "d", "a",
#'   "c", "a",
#'   "a", "a",
#'   "d", "d"
#' ), ncol = 2, byrow = TRUE)
#'
#' r.clo<-ReflexiveClosure(set, rel)
#'
#' @name ReflexiveClosure
#' @export
ReflexiveClosure <- function(set, rel) {
  if (!is.character(set) || length(set) == 0L) {
    stop("'set' must be a non-empty character vector.", call. = FALSE)
  }
  if (anyNA(set)) {
    stop("'set' must not contain NA values.", call. = FALSE)
  }
  if (anyDuplicated(set)) {
    stop("'set' contains duplicated values.", call. = FALSE)
  }
  if (!is.matrix(rel) || ncol(rel) != 2) {
    stop("'rel' must be a two-column character matrix.", call. = FALSE)
  }
  if (nrow(rel) > 0 && !is.character(rel)) {
    stop("'rel' must be a two-column character matrix.", call. = FALSE)
  }
  if (nrow(rel) > 0 && !all(unique(as.vector(rel)) %in% set)) {
    stop("'rel' contains values not belonging to 'set'.", call. = FALSE)
  }

  tryCatch({
    return(.Call(C_ReflexiveClosure, set, rel))
  }, error = function(err) {
    err_msg   <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    clean_msg <- if (length(err_parts) > 1L) trimws(err_parts[length(err_parts)]) else trimws(err_msg)
    stop(clean_msg, call. = FALSE)
  })
}
