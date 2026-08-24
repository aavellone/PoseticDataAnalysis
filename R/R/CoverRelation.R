#' @title
#' Computing the cover relation of a poset.
#'
#' @description
#' Computes the cover relation of the input poset.
#'
#' @param poset An object of S4 class `POSet`.
#' Argument `poset` must be created by using any function contained in the package aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...) .
#'
#' @return
#' A two-column matrix  \eqn{M} of character strings (element \eqn{M[i,2]} covers element \eqn{M[i,1]}).
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
#' M <- CoverRelation(pos)
#'
#' @name CoverRelation
#' @export
CoverRelation <- function(poset) {
  # Input validation
  if (!methods::is(poset, "POSet")) {
    stop("'poset' must be an object of S4 class 'POSet'")
  }

  # Call C++ implementation with error handling
  tryCatch({
    result <- .Call(C_CoverRelation, poset@ptr)
    return(result)

  }, error = function(err) {
    # Extract meaningful error message
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })

}
