#' @title
#' Computing downsets.
#'
#' @description
#' Computes the downset of a set of elements of the input poset.
#'
#' @param poset An object of S4 class `POSet`.
#' Argument `poset` must be created by using any function contained in the package aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...) .
#'
#' @param elements A vector of character strings (the names of the input elements).
#'
#' @return
#' A vector of character strings (the names of the elements of the downset).
#'
#' @examples
#' elems<- c("a", "b", "c", "d")
#'
#' dom <- matrix(c(
#'   "a", "b",
#'   "c", "b",
#'   "a", "d"
#' ), ncol = 2, byrow = TRUE)
#'
#' pos <- POSet(elements = elems, dom = dom)
#'
#' dwn <- DownsetOf(pos, c("b","d"))
#'
#' @name DownsetOf
#' @export
DownsetOf <- function(poset, elements) {
  # Input validation
  if (!methods::is(poset, "POSet")) {
    stop("'poset' must be an object of S4 class 'POSet'")
  }

  if (!is.character(elements) || length(elements) == 0) {
    stop("'elements' must be a non-empty character vector")
  }

  # Remove duplicates before passing to C++
  elements <- unique(elements)

  # Call C++ implementation with error handling
  tryCatch({
    result <- .Call(C_DownsetOf, poset@ptr, elements)
    return(result)

  }, error = function(err) {
    # Extract meaningful error message
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })
}
