#' @title
#' Checking for pre-ordering (or quasi-ordering).
#'
#' @description
#' Checks whether the input relation is a pre-order (aka, a quasi-order), i.e. if it is reflexive and transitive.
#'
#' @param set A list of character strings (the names of the elements of the set, on which the binary relation is defined).
#'
#' @param rel A two-columns character matrix, each row comprising an element (pair) of the binary relation.
#'
#' @return
#' A boolean value.
#'
#' @examples
#' set<-c("a", "b", "c", "d")
#'
#' rel <- matrix(c(
#'   "a", "b",
#'   "c", "b",
#'   "b", "a",
#'   "d", "a",
#'   "c", "a",
#'   "d", "b",
#'   "a", "a",
#'   "b", "b",
#'   "c", "c",
#'   "d", "d"
#' ), ncol = 2, byrow = TRUE)
#'
#' chk<-IsPreorder(set, rel)
#'
#' @name IsPreorder
#' @export
IsPreorder <- function(set, rel) {
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
    result <- .Call(C_isPreorder, set, rel)
    return(result)
  }, error = function(err) {
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })
}
