#' @title
#' Disjoint sum of posets.
#'
#' @description
#' Computes the disjoint sum of the input posets.
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
#' The disjoint sum poset, an object of S4 class `POSet`.
#'
#' @details
#' Let \eqn{P_1=(V_1,\leq_1),\ldots,P_k=(V_k,\leq_k)} be \eqn{k}
#' posets on disjoint ground sets.
#' Their disjoint sum is the poset \eqn{P=(V,\lhd)} having as ground
#' set the union of the input ground sets, with \eqn{a\leq b} if
#' and only if \eqn{a,b\in V_i} and \eqn{a\leq_i b} for some \eqn{i}.
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
#' dsj.sum <- DisjointSumPOSet(pos1, pos2)
#'
#'
#' @name DisjointSumPOSet
#' @export
DisjointSumPOSet <- function(poset1, poset2, ...) {
  # Input validation
  if (!methods::is(poset1, "POSet")) {
    stop("'poset1' must be an object of S4 class 'POSet'")
  }

  if (!methods::is(poset2, "POSet")) {
    stop("'poset2' must be an object of S4 class 'POSet'")
  }


  # Validate additional posets in ...
  additional_posets <- list(...)
  if (length(additional_posets) > 0) {
    for (i in seq_along(additional_posets)) {
      if (!methods::is(additional_posets[[i]], "POSet")) {
        stop(sprintf("Argument %d in '...' must be an object of S4 class 'POSet'", i + 2))
      }
    }
  }

  # Build poset pointer list
  posets <- c(poset1@ptr, poset2@ptr)
  if (length(additional_posets) > 0) {
    for (p in additional_posets) {
      posets <- c(posets, p@ptr)
    }
  }

  # Call C++ implementation with error handling
  tryCatch({
    ptr <- .Call(C_BuildDisjointSumPOSet, posets)
    result <- methods::new("POSet", ptr = ptr)
    return(result)

  }, error = function(err) {
    # Extract meaningful error message
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })
}
