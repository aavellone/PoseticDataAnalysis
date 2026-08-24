#' @title Computing the intersection of a collection of posets
#'
#' @description
#' Computes the poset \eqn{(V, \leq_{\cap})=(V, \leq_1)\cap\cdots\cap(V,\leq_k)}.
#'
#' @param poset1 An object of S4 class `POSet`.
#' Argument `poset1` must be created by using any function contained in the package aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...).
#' @param poset2 An object of S4 class `POSet`.
#' Argument `poset2` must be created by using any function contained in the package aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...).
#' @param ... Optional additional objects of S4 class `POSet`.
#' Optional arguments must be created by using any function contained in the package aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...).
#'
#' @return
#' The intersection poset, an object of S4 class `POSet`.
#'
#' @details
#' Let \eqn{P_1 = (V, \leq_1),\ldots, P_k = (V, \leq_k)} be \eqn{k} posets on the same set \eqn{V}.
#' The intersection poset \eqn{P_{\cap}=P_1 \cap\cdots\cap P_k} is the poset \eqn{(V, \leq_{\cap})} where
#' \eqn{a\leq_{\cap} b} if and only if \eqn{a\leq_i b} for all \eqn{i=1,\ldots, k}.
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
#'   "b", "c",
#'   "b", "d"
#' ), ncol = 2, byrow = TRUE)
#'
#' pos1 <- POSet(elements = elems, dom = dom1)
#'
#' pos2 <- POSet(elements = elems, dom = dom2)
#'
#' pos_int <- IntersectionPOSet(pos1, pos2)
#'
#'
#' @name IntersectionPOSet
#' @export
IntersectionPOSet <- function(poset1, poset2, ...) {
  # Input validation
  if (!methods::is(poset1, "POSet")) {
    stop("'poset1' must be an object of S4 class 'POSet'", call. = FALSE)
  }
  if (!methods::is(poset2, "POSet")) {
    stop("'poset2' must be an object of S4 class 'POSet'", call. = FALSE)
  }

  additional_posets <- list(...)
  for (i in seq_along(additional_posets)) {
    if (!methods::is(additional_posets[[i]], "POSet")) {
      stop(sprintf("Argument %d in '...' must be an object of S4 class 'POSet'", i + 2L),
           call. = FALSE)
    }
  }
  # Build pointer list in one shot, avoiding repeated c() reallocations
  all_posets <- c(list(poset1, poset2), additional_posets)
  posets <- lapply(all_posets, methods::slot, "ptr")

  tryCatch({
    result <- methods::new("POSet", ptr = .Call(C_BuildIntersectionPOSet, posets))
    return(result)
  }, error = function(err) {
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })
}
