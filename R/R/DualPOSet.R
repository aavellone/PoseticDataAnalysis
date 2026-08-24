#' @title
#' Dual of a poset.
#'
#' @description
#' Computes the dual of the input poset.
#'
#' @param poset An object of S4 class `POSet`.
#' Argument `poset` must be created by using any function contained in the package
#' aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...).
#'
#' @return
#' The dual of the input poset, an object of S4 class `POSet`.
#'
#' @details
#' Let \eqn{P=(V,\leq)} be a poset. Then its dual \eqn{P_d=(V,\leq_d)} is defined
#' by \eqn{a\leq_d b} if and only if \eqn{b\leq a} in \eqn{P}.
#' In other words, the dual of \eqn{P} is obtained by reversing its dominances.
#'
#' @examples
#' elems <- c("a", "b", "c", "d")
#'
#' doms <- matrix(c(
#'   "a", "b",
#'   "c", "b",
#'   "b", "d"
#' ), ncol = 2, byrow = TRUE)
#'
#' pos1 <- POSet(elements = elems, dom = doms)
#'
#' dual <- DualPOSet(pos1)
#'
#' @name DualPOSet
#' @export
DualPOSet <- function(poset) {
  # Input validation
  if (!methods::is(poset, "POSet")) {
    stop("'poset' must be an object of S4 class 'POSet'")
  }

  # Call C++ implementation with error handling
  tryCatch({
    ptr <- .Call(C_BuildDualPOSet, poset@ptr)
    result <- methods::new("POSet", ptr = ptr)
    return(result)

  }, error = function(err) {
    # Extract meaningful error message
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })
}
