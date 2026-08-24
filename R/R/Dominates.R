#' @title
#' Checking whether one element dominates another.
#'
#' @description
#' Given two elements \eqn{a} and \eqn{b} of \eqn{V}, checks whether \eqn{a\leq b} in poset \eqn{(V,\leq)}.
#'
#' @param poset An object of S4 class `POSet`.
#' Argument `poset` must be created by using any function contained in the package aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...) .
#' @param element1 A character string or vector (the name(s) of poset element(s)).
#' @param element2 A character string or vector (the name(s) of poset element(s)).
#'
#' @return
#' A boolean value.
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
#' chk <- Dominates(pos, "a", "d")
#'
#' @name Dominates
#' @export
Dominates <- function(poset, element1, element2) {
  # Input validation
  if (!methods::is(poset, "POSet")) {
    stop("'poset' must be an object of S4 class 'POSet'")
  }

  if (!is.character(element1)) {
    stop("'element1' must be a character string or vector")
  }

  if (!is.character(element2)) {
    stop("'element2' must be a character string or vector")
  }

  if (length(element1) != length(element2)) {
    stop("'element1' and 'element2' must have the same length")
  }

  # Call C++ implementation with error handling
  tryCatch({
    result <- .Call(C_Dominates, poset@ptr, element1, element2)
    return(result)

  }, error = function(err) {
    # Extract meaningful error message
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })

}
