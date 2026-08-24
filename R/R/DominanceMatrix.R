#' @title
#' Computing the dominance matrix.
#'
#' @param poset An object of S4 class `POSet`.
#' Argument `poset` must be created by using any function contained in the package aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...) .
#'
#' @description
#' Computes the dominance matrix of the input poset.
#'
#' @return
#' An \eqn{n\times n} boolean matrix \eqn{Z}, where \eqn{n} is the number of poset elements, with \eqn{Z[i,j]=TRUE},
#' if and only if the j-th poset element weakly dominates (\eqn{\leq}) the i-th element, in the input order relation.
#'
#' @examples
#' elems <- c("a", "b", "c", "d")
#' dom <- matrix(c(
#'   "a", "b",
#'   "c", "b",
#'   "b", "d"
#' ), ncol = 2, byrow = TRUE)
#'
#' pos <- POSet(elements = elems, dom = dom)
#'
#' Z <- DominanceMatrix(pos)
#'
#' @name DominanceMatrix
#' @export
DominanceMatrix <- function(poset) {
  # Input validation
  if (!methods::is(poset, "POSet")) {
    stop("'poset' must be an object of S4 class 'POSet'")
  }

  # Call C++ implementation with error handling
  tryCatch({
    result <- .Call(C_IncidenceMatrix, poset@ptr)
    return(result)

  }, error = function(err) {
    # Extract meaningful error message
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })
}
