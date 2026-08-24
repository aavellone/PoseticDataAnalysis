#' @title
#' Checking binary relation reflexivity.
#'
#' @description
#' Checks whether the input binary relation is reflexive.
#'
#' @param set A list of strings (the names of the elements of the set, on which the binary relation is defined).
#'
#' @param rel A two-columns character matrix, each row comprising an element (pair) of the binary relation.
#'
#' @return
#' A boolean value.
#'
#' @examples
#' set<-c("a", "b", "c", "d", "e")
#'
#' rel <- matrix(c(
#'   "a", "b",
#'   "c", "b",
#'   "b", "d",
#'   "a", "a",
#'   "b", "b",
#'   "c", "c"
#' ), ncol = 2, byrow = TRUE)
#'
#' chk <- IsReflexive(set, rel)
#'
#' @name IsReflexive
#' @export
IsReflexive <- function(set, rel) {
  if (!is.character(set)) {
    stop("'set' must be a character vector.", call. = FALSE)
  }
  if (length(unique(set)) != length(set)) {
    stop("'set' contains duplicated values.", call. = FALSE)
  }
  if (!is.matrix(rel) || ncol(rel) != 2L) {
    stop("'rel' must be a two-column character matrix.", call. = FALSE)
  }
  if (ncol(rel) * nrow(rel) > 0L && !is.character(rel)) {
    stop("'rel' must be a two-column character matrix.", call. = FALSE)
  }
  if (!all(unique(as.vector(rel)) %in% set)) {
    stop("'rel' contains values not belonging to 'set'.", call. = FALSE)
  }

  tryCatch({
    result <- .Call(C_isReflexive, set, rel)
    return(result)
  }, error = function(err) {
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })
}
