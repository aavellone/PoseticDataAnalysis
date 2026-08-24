#' @title
#' Extracting the comparability set of a poset element.
#'
#' @description
#' Extracts the elements comparable with the input element, in the poset.
#'
#' @param poset An object of S4 class `POSet`.
#' Argument `poset` must be created by using any function contained in the package aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...) .
#'
#' @param element A character string (the name of the input poset element).
#'
#' @return
#' A vector of character strings (the names of the poset elements comparable to the input element).
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
#' cmp <- ComparabilitySetOf(pos, "a")
#'
#' @name ComparabilitySetOf
#' @export
ComparabilitySetOf <- function(poset, element) {
  # Input validation
  if (!methods::is(poset, "POSet")) {
    stop("'poset' must be an object of S4 class 'POSet'")
  }

  if (!is.character(element) || length(element) != 1) {
    stop("'element' must be a single character string")
  }

  # Call C++ implementation with error handling
  tryCatch({
    result <- .Call(C_ComparabilitySetOf, poset@ptr, element)
    return(result)

  }, error = function(err) {
    # Extract meaningful error message
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })
}
