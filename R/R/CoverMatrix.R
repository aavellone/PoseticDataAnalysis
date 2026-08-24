#' @title
#' Computing the cover matrix of a poset.
#'
#' @description
#' Computes the cover matrix of the input poset.
#'
#' @param poset An object of S4 class POSet.
#' Argument `poset` must be created by using any function contained in the package aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...) .
#'
#' @return
#' A square boolean matrix \eqn{C} (\eqn{C[i,j]=TRUE} if and only if the \eqn{j}-th element of the input poset covers the \eqn{i}-th).
#'
#' @examples
#' elems <- c("a", "b", "c", "d")
#'
#' dom <- matrix(c(
#'   "a", "b",
#'   "c", "b",
#'   "b", "d"
#' ), ncol = 2, byrow = TRUE)
#'
#' pos <- POSet(elements = elems, dom = dom)
#'
#' C <- CoverMatrix(pos)
#'
#' @name CoverMatrix
#' @export
CoverMatrix <- function(poset) {
  # Input validation
  if (!methods::is(poset, "POSet")) {
    stop("'poset' must be an object of S4 class 'POSet'")
  }

  # Call C++ implementation with error handling
  tryCatch({
    result <- .Call(C_CoverMatrix, poset@ptr)
    return(result)

  }, error = function(err) {
    # Extract meaningful error message
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })
}
