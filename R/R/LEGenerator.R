#' @title
#' Generator of all the linear extensions of a poset.
#'
#' @description
#' Creates an object of S4 class `LEGenerator` to generate all of the linear extensions of a given poset.
#' Actually, this function does not generate the linear extensions, but just the object that will generate them
#' by using function `LEGet`.
#'
#' @param poset An object of S4 class `POSet` representing the poset whose linear extensions are generated.
#' Argument `poset` must be created by using any function contained in the package aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...) .
#'
#' @return
#' An S4 class object `LEGenerator`.
#'
#' @examples
#' el <- c("a", "b", "c", "d")
#'
#' dom <- matrix(c(
#'   "a", "b",
#'   "c", "b",
#'   "b", "d"
#' ), ncol = 2, byrow = TRUE)
#'
#' pos <- POSet(elements = el, dom = dom)
#'
#' LEgen <- LEGenerator(pos)
#'
#' @name LEGenerator
#' @export
LEGenerator <- function(poset) {
  if (!methods::is(poset, "POSet")) {
    stop("'poset' must be of class POSet.", call. = FALSE)
  }

  tryCatch({
    ptr <- .Call(C_BuildLEGenerator, poset@ptr)
    if (is.null(ptr)) {
      stop("The C++ generator returned a NULL pointer.", call. = FALSE)
    }
    return(methods::new("LEGenerator", ptr = ptr))
  }, error = function(err) {
    err_msg   <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    clean_msg <- if (length(err_parts) > 1L) trimws(err_parts[length(err_parts)]) else trimws(err_msg)
    stop(clean_msg, call. = FALSE)
  })
}
