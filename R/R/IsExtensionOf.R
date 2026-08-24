#' @title
#' Checking poset extensions.
#'
#' @description
#' Checks whether poset1 is an extension of poset2.
#'
#' @param poset1 An object of S4 class 'POSet'.
#' Argument `poset1` must be created by using any function contained in the package aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...) .
#'
#' @param poset2 An object of S4 class 'POSet'.
#' Argument `poset2` must be created by using any function contained in the package aimed at building object of S4cl ass `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...) .
#'
#' @return
#' A boolean value.
#'
#' @examples
#' elems <- c("a", "b", "c", "d")
#'
#' dom1 <- matrix(c(
#'   "a", "b",
#'   "c", "b",
#'   "b", "d"
#' ), ncol = 2, byrow = TRUE)
#'
#' dom2 <- matrix(c(
#'   "a", "b",
#'   "c", "b",
#'   "b", "d",
#'   "a", "c"
#' ), ncol = 2, byrow = TRUE)
#'
#' pos1 <- POSet(elements = elems, dom = dom1)
#' pos2 <- POSet(elements = elems, dom = dom2)
#'
#' chk <- IsExtensionOf(pos1, pos2)
#'
#' @name IsExtensionOf
#' @export
IsExtensionOf <- function(poset1, poset2) {
  if (!methods::is(poset1, "POSet")) {
    stop("'poset1' must be of class POSet.", call. = FALSE)
  }
  if (!methods::is(poset2, "POSet")) {
    stop("'poset2' must be of class POSet.", call. = FALSE)
  }

  tryCatch({
    result <- .Call(C_IsExtensionOf, poset1@ptr, poset2@ptr)
    return(result)
  }, error = function(err) {
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })
}
