#' @title
#' Constructing the product of posets.
#'
#' @description
#' Constructs the product poset \eqn{(V, \leq_{prd})}, starting from a collection of posets.
#'
#' @param poset1 An object of S4 class `POSet`.
#' Argument `poset1` must be created by using any function contained in the package
#' aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...) .
#' @param poset2 An object of S4 class `POSet`.
#' Argument `poset2` must be created by using any function contained in the package
#' aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...) .
#' @param ... Optional additional objects of S4 class `POSet`.
#' Optional arguments must be created by using any function contained in the package
#' aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...) .
#'
#' @return
#' The product poset, an object of S4 class `ProductPOSet` (subclass of `POSet`).
#'
#' @details
#' Let \eqn{P_1 = (V_1, \leq_1), ..., P_k = (V_k, \leq_k)} be a collection of posets.
#' The product poset \eqn{P=P_1 \times...\times P_k} is the poset \eqn{(V, \leq_{prd})} where
#' \eqn{V=V_1\times...\times V_k} and given \eqn{(a_1, ..., a_k)\in V}
#' and \eqn{(b_1, ..., b_k)\in V},  \eqn{(a_1, ..., a_k)\leq_{prd} (b_1, ..., b_k)} if and only if
#' \eqn{a_i\leq_i b_i} for all \eqn{i=1, ..., k}.
#'
#' @examples
#' elems1 <- c("a", "b", "c", "d")
#' elems2 <- c("x", "y", "z")
#' elems3 <- c("q", "r")
#'
#' dom <- matrix(c(
#'   "a", "b",
#'   "c", "b",
#'   "b", "d"
#' ), ncol = 2, byrow = TRUE)
#'
#' p1 <- POSet(elements = elems1, dom = dom)
#' p2 <- LinearPOSet(elements = elems2)
#' p3 <- LinearPOSet(elements = elems3)
#'
#' prd12 <- ProductPOSet(p1, p2)
#'
#' prd123 <- ProductPOSet(p1, p2, p3)
#'
#' @name ProductPOSet
#' @export
ProductPOSet <- function(poset1, poset2, ...) {
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
    ptr <- .Call(C_BuildProductPOSet, posets)
    if (is.null(ptr)) {
      stop("The C++ constructor returned a NULL pointer.", call. = FALSE)
    }
    return(methods::new("ProductPOSet", ptr = ptr))
  }, error = function(err) {
    err_msg   <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    clean_msg <- if (length(err_parts) > 1L) trimws(err_parts[length(err_parts)]) else trimws(err_msg)
    stop(clean_msg, call. = FALSE)
  })
}
