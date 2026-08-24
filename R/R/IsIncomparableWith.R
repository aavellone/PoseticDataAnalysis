#' @title
#' Checking incomparability between two elements of a poset.
#'
#' @description
#' Checks whether two elements \eqn{a} and \eqn{b} of \eqn{V} are incomparable, in the input poset \eqn{(V,\leq)}.
#'
#' @param poset An object of S4 class `POSet`.
#' Argument `poset` must be created by using any function contained in the package aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...) .
#' @param element1 A character string (the name of a poset element).
#' @param element2 A character string (the name of a poset element).
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
#' chk <- IsIncomparableWith(pos, "a", "d")
#'
#' @name IsIncomparableWith
#' @export
IsIncomparableWith <- function(poset, element1, element2) {
  if (!methods::is(poset, "POSet")) {
    stop("'poset' must be of class POSet.", call. = FALSE)
  }
  if (!is.character(element1)) {
    stop("'element1' must be a character vector.", call. = FALSE)
  }
  if (!is.character(element2)) {
    stop("'element2' must be a character vector.", call. = FALSE)
  }
  if (length(element1) != length(element2)) {
    stop("'element1' and 'element2' must have the same length.", call. = FALSE)
  }

  tryCatch({
    result <- .Call(C_IsIncomparableWith, poset@ptr, element1, element2)
    return(result)
  }, error = function(err) {
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })
}
