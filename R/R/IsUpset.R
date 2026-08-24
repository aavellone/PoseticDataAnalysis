#' @title
#' Checking upsets.
#'
#' @description
#' Checks whether the input elements form an upset, in the input poset.
#'
#' @param poset An object of S4 class `POSet`.
#' Argument `poset` must be created by using any function contained in the package aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...) .
#'
#' @param elements A vector of character strings (the names of the input elements).
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
#' chk <- IsUpset(pos, c("a", "b", "c"))
#'
#'
#' @name IsUpset
#' @export
IsUpset <- function(poset, elements) {
  if (!methods::is(poset, "POSet")) {
    stop("'poset' must be of class POSet.", call. = FALSE)
  }
  if (!is.character(elements)) {
    stop("'elements' must be a character vector.", call. = FALSE)
  }

  tryCatch({
    result <- .Call(C_IsUpset, poset@ptr, unique(elements))
    return(result)
  }, error = function(err) {
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })
}
