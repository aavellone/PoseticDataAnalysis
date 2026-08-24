#' @title
#' Computing lexicographic product orders.
#'
#' @description
#' Computes the lexicographic product order of the input partial orders.
#'
#' @param poset1 An object of S4 class `POSet`.
#' Argument `poset1` must be created by using any function contained in the package
#' aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...).
#' @param poset2 An object of S4 class `POSet`.
#' Argument `poset2` must be created by using any function contained in the package
#' aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...).
#' @param ... Optional additional objects of S4 class `POSet`.
#' Optional arguments must be created by using any function contained in the package
#' aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...).
#'
#' @return
#' The lexicographic product poset, an object of S4 class `LexicographicProductPOSet` (subclass of `POSet`).
#'
#' @details
#' Let \eqn{P_1 = (V_1, \leq_1),\cdots, P_k = (V_k, \leq_k)} be \eqn{k} posets.
#' The lexicographic product poset \eqn{P_{lxprd}=(V, \leq_{lxprd})} has ground
#' set the Cartesian product of the input ground sets, with
#' \eqn{(a_1,\ldots,a_k)\leq_{lxprd} (b_1,\ldots,b_k)} if and only if \eqn{a_1\leq_1 b_1},
#' or there exists \eqn{j} such that \eqn{a_i=b_i} for \eqn{i<j} and \eqn{a_j\leq_j b_j}.
#'
#' @examples
#' elems1 <- c("a", "b", "c", "d", "e")
#' elems2 <- c("f", "g", "h")
#'
#' dom1 <- matrix(c(
#'   "a", "b",
#'   "c", "b",
#'   "d", "b"
#' ), ncol = 2, byrow = TRUE)
#'
#' dom2 <- matrix(c(
#'   "g", "f",
#'   "h", "f"
#' ), ncol = 2, byrow = TRUE)
#'
#' pos1 <- POSet(elements = elems1, dom = dom1)
#'
#' pos2 <- POSet(elements = elems2, dom = dom2)
#'
#' #Lexicographic product of pos1 and pos2
#' lex.prod <- LexicographicProductPOSet(pos1, pos2)
#'
#' @name LexicographicProductPOSet
#' @export
LexicographicProductPOSet <- function(poset1, poset2, ...) {
  if (!methods::is(poset1, "POSet")) {
    stop("'poset1' must be of class POSet.", call. = FALSE)
  }
  if (!methods::is(poset2, "POSet")) {
    stop("'poset2' must be of class POSet.", call. = FALSE)
  }
  dots <- list(...)
  for (i in seq_along(dots)) {
    if (!methods::is(dots[[i]], "POSet")) {
      stop(paste0("Argument '...[", i, "]' must be of class POSet."), call. = FALSE)
    }
  }

  # Build the ordered list of C++ pointers
  posets <- c(list(poset1@ptr), list(poset2@ptr), lapply(dots, function(p) p@ptr))

  tryCatch({
    ptr <- .Call(C_BuildLexicographicProductPOSet, posets)
    if (is.null(ptr)) {
      stop("The C++ constructor returned a NULL pointer.", call. = FALSE)
    }
    return(methods::new("LexicographicProductPOSet", ptr = ptr))
  }, error = function(err) {
    err_msg   <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    clean_msg <- if (length(err_parts) > 1L) trimws(err_parts[length(err_parts)]) else trimws(err_msg)
    stop(clean_msg, call. = FALSE)
  })
}
