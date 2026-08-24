#' @title
#' Checking minimality.
#'
#' @description
#' Checks whether the input element is minimal in the input poset.
#'
#' @param poset An object of S4-class `POSet`.
#' Argument `poset` must be created by using any function contained in the package aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...) .

#' @param element A character string (the name of the input element).
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
#' chk <- IsMinimal(pos, "a")

#'
#' @name IsMinimal
#' @export
IsMinimal <- function(poset, element) {
  if (!methods::is(poset, "POSet")) {
    stop("'poset' must be of class POSet.", call. = FALSE)
  }
  if (!is.character(element)) {
    stop("'element' must be a character string.", call. = FALSE)
  }
  if (length(element) != 1L) {
    stop("'element' must be a single value.", call. = FALSE)
  }

  tryCatch({
    result <- .Call(C_IsMinimal, poset@ptr, element)
    return(result)
  }, error = function(err) {
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })
}
