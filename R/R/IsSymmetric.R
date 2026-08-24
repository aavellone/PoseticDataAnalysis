#' @title
#' Checking binary relation symmetry.
#'
#' @description
#' Checks whether the input binary relation is symmetric.
#'
#' @param rel A two-columns character matrix, each row comprising an element (pair) of the binary relation.
#'
#' @return
#' A boolean value.
#'
#' @examples
#' rel <- matrix(c(
#'   "a", "b",
#'   "c", "b",
#'   "b", "d",
#'   "b", "a",
#'   "b", "c",
#'   "d", "b"
#' ), ncol = 2, byrow = TRUE)
#'
#' chk<-IsSymmetric(rel)
#'
#' @name IsSymmetric
#' @export
IsSymmetric <- function(rel) {
  if (!is.matrix(rel) || ncol(rel) != 2L) {
    stop("'rel' must be a two-column character matrix.", call. = FALSE)
  }
  if (ncol(rel) * nrow(rel) > 0L && !is.character(rel)) {
    stop("'rel' must be a two-column character matrix.", call. = FALSE)
  }

  tryCatch({
    result <- .Call(C_isSymmetric, rel)
    return(result)
  }, error = function(err) {
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })
}
