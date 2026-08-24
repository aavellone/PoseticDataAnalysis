#' @title
#' Computing upsets.
#'
#' @description
#' Computes the upset of a set of elements of the input poset.
#'
#' @param poset An object of S4 class `POSet`.
#' Argument `poset` must be created by using any function contained in the package aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...) .
#' @param elements a vector of character strings (the names of the input elements).
#'
#' @return
#' A vector of character strings (the names of the poset elements in the upset).
#'
#' @examples
#' elems<- c("a", "b", "c", "d")
#'
#' dom <- matrix(c(
#'   "a", "b",
#'   "c", "b",
#'   "b", "d"
#' ), ncol = 2, byrow = TRUE)
#'
#' pos <- POSet(elements = elems, dom = dom)
#'
#' up <- UpsetOf(pos, c("a","c"))
#'
#' @name UpsetOf
#' @export
UpsetOf <- function(poset, elements) {
  if (!methods::is(poset, "POSet")) {
    stop("'poset' must be of class POSet.", call. = FALSE)
  }
  if (!is.character(elements) || length(elements) == 0L) {
    stop("'elements' must be a non-empty character vector.", call. = FALSE)
  }

  # Deduplicate before passing to C++
  elements <- unique(elements)

  tryCatch({
    return(.Call(C_UpsetOf, poset@ptr, elements))
  }, error = function(err) {
    err_msg   <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    clean_msg <- if (length(err_parts) > 1L) trimws(err_parts[length(err_parts)]) else trimws(err_msg)
    stop(clean_msg, call. = FALSE)
  })
}
