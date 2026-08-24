#' @title
#' Extracting the order relation of a poset.
#'
#' @description
#' Extracts the order relation from the input poset.
#'
#' @param poset An object of S4 class `POSet`.
#' Argument `poset` must be created by using any function contained in the package aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...) .
#'
#' @return
#' A two-column matrix \eqn{M} of character strings (element \eqn{M[i,2]} dominates element \eqn{M[i,1]}).
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
#' M <- OrderRelation(pos)
#'
#' @name OrderRelation
#' @export
OrderRelation <- function(poset) {
  if (!methods::is(poset, "POSet")) {
    stop("'poset' must be of class POSet.", call. = FALSE)
  }
  tryCatch({
    return(.Call(C_OrderRelation, poset@ptr))
  }, error = function(err) {
    err_msg   <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    clean_msg <- if (length(err_parts) > 1L) trimws(err_parts[length(err_parts)]) else trimws(err_msg)
    stop(clean_msg, call. = FALSE)
  })
}
