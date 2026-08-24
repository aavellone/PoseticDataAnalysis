#' @title
#' Lifting posets.
#'
#' @description
#' Lifts the input poset, i.e. adds a (possibly new) bottom element to it.
#'
#' @param poset An object of S4 class `POSet`.
#' Argument `poset` must be created by using any function contained in the package aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...).
#'
#' @param element A character string (the name of the added bottom).
#'
#' @return
#' The lifted poset, an object of S4 class `POSet`.
#'
#' @examples
#' elems <- c("a", "b", "c", "d")
#'
#' doms <- matrix(c(
#'   "a", "b",
#'   "c", "b",
#'   "b", "d"
#' ), ncol = 2, byrow = TRUE)
#'
#' pos <- POSet(elements = elems, dom = doms)
#'
#' #Lifting
#' lifted.pos <- LiftingPOSet(pos, "bot")
#'
#' @name LiftingPOSet
#' @export
LiftingPOSet <- function(poset, element) {
  if (!methods::is(poset, "POSet")) {
    stop("'poset' must be of class POSet.", call. = FALSE)
  }
  if (!is.character(element) || length(element) != 1L || is.na(element)) {
    stop("'element' must be a single non-NA character string.", call. = FALSE)
  }

  tryCatch({
    ptr <- .Call(C_BuildLiftingPOSet, poset@ptr, element)
    if (is.null(ptr)) {
      stop("The C++ constructor returned a NULL pointer.", call. = FALSE)
    }
    return(methods::new("POSet", ptr = ptr))
  }, error = function(err) {
    err_msg   <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    clean_msg <- if (length(err_parts) > 1L) trimws(err_parts[length(err_parts)]) else trimws(err_msg)
    stop(clean_msg, call. = FALSE)
  })
}
