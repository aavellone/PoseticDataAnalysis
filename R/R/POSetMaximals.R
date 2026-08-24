#' @title
#' Computing the maximal elements of a poset.
#'
#' @description
#' Computes the maximal elements of the input poset, i.e. those elements being strictly dominated by no other elements.
#'
#' @param poset An object of S4 class `POSet`.
#' Argument `poset` must be created by using any function contained in the package aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...) .
#'
#' @return
#' A vector of character strings (the names of the maximal elements).
#'
#' @examples
#' elems <- c("a", "b", "c", "d")
#'
#' dom <- matrix(c(
#'   "a", "b"
#'   ,
#'   "c", "b",
#'   "b", "d"
#' ), ncol = 2, byrow = TRUE)
#'
#' pos <- POSet(elements = elems, dom = dom)
#'
#' mxs <- POSetMaximals(pos)
#'
#' @name POSetMaximals
#' @export
POSetMaximals <- function(poset) {
  if (!methods::is(poset, "POSet")) {
    stop("'poset' must be of class POSet.", call. = FALSE)
  }
  tryCatch({
    return(.Call(C_Maximal, poset@ptr))
  }, error = function(err) {
    err_msg   <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    clean_msg <- if (length(err_parts) > 1L) trimws(err_parts[length(err_parts)]) else trimws(err_msg)
    stop(clean_msg, call. = FALSE)
  })
}
