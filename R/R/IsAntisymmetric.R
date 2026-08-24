#' @title Checking binary relation antisymmetry.
#'
#' @description
#' Checks whether the input binary relation is antisymmetric.
#'
#' @param rel A two-column character matrix, each row comprising an element (pair) of the binary relation.
#'
#' @return A logical value: TRUE if the relation is antisymmetric, FALSE otherwise.
#'
#' @examples
#' rel <- matrix(c(
#'   "a", "b",
#'   "c", "b",
#'   "b", "d",
#'   "a", "a"
#' ), ncol = 2, byrow = TRUE)
#'
#' chk <- IsAntisymmetric(rel)
#'
#' @name IsAntisymmetric
#' @export
IsAntisymmetric <- function(rel) {
  if (!is.matrix(rel) || !is.character(rel) || ncol(rel) != 2L) {
    stop("'rel' must be a two-column character matrix.", call. = FALSE)
  }

  tryCatch({
    result <- .Call(C_isAntisymmetric, rel)
    return(result)
  }, error = function(err) {
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })
}
