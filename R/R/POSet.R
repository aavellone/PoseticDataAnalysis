#' @title
#' Constructing a Partially Ordered Set.
#'
#' @description
#' Constructs an object of class `POSet`, representing a partially ordered set (poset) \eqn{P=(V,\leq)}.
#'
#' @param elements A vector of character strings (the labels of the elements of the ground set \eqn{V}).
#' @param dom Two-columns matrix of element labels, representing the dominances in the order relation
#' \eqn{\leq}. The generic `k-th` row of `dom` contains a pair of elements of \eqn{V}, with
#' `dom[k, 1] `\eqn{ \leq} `dom[k, 2]`.
#'
#' @return
#' An object \eqn{(V, \leq)} of S4 class `POSet`, where \eqn{V} is the ground set and \eqn{\leq} is the partial order relation on it
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
#' @name POSet
#' @export
POSet <- function(elements, dom=matrix(character(0), ncol = 2)) {
  if (!is.character(elements) || length(elements) == 0L) {
    stop("'elements' must be a non-empty character vector.", call. = FALSE)
  }
  if (anyNA(elements)) {
    stop("'elements' must not contain NA values.", call. = FALSE)
  }
  if (anyDuplicated(elements)) {
    stop("'elements' contains duplicated values.", call. = FALSE)
  }
  if (!is.matrix(dom) || ncol(dom) != 2) {
    stop("'dom' must be a two-column matrix.", call. = FALSE)
  }
  if (nrow(dom) > 0 && !is.character(dom)) {
    stop("'dom' must be a two-column character matrix.", call. = FALSE)
  }
  if (nrow(dom) > 0 && !all(unique(as.vector(dom)) %in% elements)) {
    stop("'dom' contains values not belonging to 'elements'.", call. = FALSE)
  }

  tryCatch({
    ptr <- .Call(C_BuildPOSet, elements, dom,
                 PACKAGE = "poseticDataAnalysis")
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
