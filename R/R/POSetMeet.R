#' @title
#' Computing meet (greatest lower bound).
#'
#' @description
#' The function computes the meet (if existing) of a set of elements, in the input poset.
#'
#' @param poset An object of S4 class `POSet`.
#' Argument `poset` must be created by using any function contained in the package
#' aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...).
#'
#' @param elements A character vector (the names of some poset elements).
#'
#' @return
#' A character string (the name of the meet).
#'
#' @examples
#' elems <- c("a", "b", "c", "d", "e")
#'
#' doms <- matrix(c(
#'   "a", "b",
#'   "c", "b",
#'   "c", "e",
#'   "b", "d",
#'   "a", "d",
#'   "c", "d",
#'   "a", "a",
#'   "b", "b",
#'   "c", "c",
#'   "d", "d"
#' ), ncol = 2, byrow = TRUE)
#'
#' pos <- POSet(elements = elems, dom = doms)
#'
#' glb<-POSetMeet(pos, c("b", "e"))
#'
#' @name POSetMeet
#' @export
POSetMeet <- function(poset, elements) {
  if (!methods::is(poset, "POSet")) {
    stop("'poset' must be of class POSet.", call. = FALSE)
  }
  if (!is.character(elements) || length(elements) == 0L) {
    stop("'elements' must be a non-empty character vector.", call. = FALSE)
  }
  if (anyDuplicated(elements)) {
    stop("'elements' contains duplicated values.", call. = FALSE)
  }
  tryCatch({
    return(.Call(C_Meet, poset@ptr, elements))
  }, error = function(err) {
    err_msg   <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    clean_msg <- if (length(err_parts) > 1L) trimws(err_parts[length(err_parts)]) else trimws(err_msg)
    stop(clean_msg, call. = FALSE)
  })
}
