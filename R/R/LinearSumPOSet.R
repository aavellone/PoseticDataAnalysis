#' @title
#' Linear sum of posets.
#'
#' @description
#' Computes the linear sum of the input posets.
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
#' The linear sum poset, an object of S4 class `POSet`.
#'
#' @details
#' Let \eqn{P_1=(V_1,\leq_1),\ldots,P_k=(V_k,\leq_k)} be \eqn{k} posets on disjoint ground sets.
#' Their linear sum is the poset \eqn{P=(V,\lhd)} having as ground set the union of the input ground sets,
#' with \eqn{a\leq b} if and only if \eqn{a\leq_i b} for some \eqn{i}, or \eqn{a\in V_i} and \eqn{b\in V_j},
#' with \eqn{i<j}. In other words, the linear sum is obtained by stacking the input posets from bottom,
#' and making all of the minimal elements of \eqn{P_i} covering all of the maximal elements of \eqn{P_{i-1}} (\eqn{i>1}).
#'
#' @examples
#' elems1 <- c("a", "b", "c", "d")
#' elems2 <- c("e", "f", "g", "h")
#'
#' dom1 <- matrix(c(
#'   "a", "b",
#'   "c", "b",
#'   "b", "d"
#' ), ncol = 2, byrow = TRUE)
#'
#' dom2 <- matrix(c(
#'   "e", "f",
#'   "g", "h",
#'   "h", "f"
#' ), ncol = 2, byrow = TRUE)
#'
#' pos1 <- POSet(elements = elems1, dom = dom1)
#'
#' pos2 <- POSet(elements = elems2, dom = dom2)
#'
#' #Linear sum of pos1 and pos2
#' lin.sum <- LinearSumPOSet(pos1, pos2)
#'
#' @name LinearSumPOSet
#' @export
LinearSumPOSet <- function(poset1, poset2, ...) {
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

  # Build an ordered list of C++ pointers
  posets <- c(list(poset1@ptr), list(poset2@ptr), lapply(dots, function(p) p@ptr))

  tryCatch({
    ptr <- .Call(C_BuildLinearSumPOSet, posets)
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
