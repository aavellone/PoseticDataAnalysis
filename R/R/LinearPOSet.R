#' @title
#' Constructing a Linearly Ordered Set.
#'
#' @description
#' Constructs a linearly (or completely, or totally) ordered set \eqn{(V,\leq_{lin})}, starting from set \eqn{V}.
#'
#' @param elements A character string vector containing the labels of the elements of \eqn{V} in ascending
#' order according to \eqn{\leq_{lin}}, i.e. such that `elements[h]` \eqn{\leq_{lin}} `elements[k]` if and only if `h` \eqn{\leq} `k`.
#'
#' @return
#' An object of S4 class `LinearPOSet` (subclass of `POSet`)
#'
#' @examples
#' elems <- c("a", "b", "c", "d")
#' linpos <- LinearPOSet(elements = elems)
#'
#' @name LinearPOSet
#' @export
LinearPOSet <- function(elements) {
  if (!is.character(elements) || length(elements) == 0L) {
    stop("'elements' must be a non-empty character vector.", call. = FALSE)
  }
  if (anyNA(elements)) {
    stop("'elements' must not contain NA values.", call. = FALSE)
  }
  if (anyDuplicated(elements)) {
    stop("'elements' contains duplicated values.", call. = FALSE)
  }

  tryCatch({
    ptr <- .Call(C_BuildLinearPOSet, elements)
    if (is.null(ptr)) {
      stop("The C++ constructor returned a NULL pointer.", call. = FALSE)
    }
    return(methods::new("LinearPOSet", ptr = ptr))
  }, error = function(err) {
    err_msg   <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    clean_msg <- if (length(err_parts) > 1L) trimws(err_parts[length(err_parts)]) else trimws(err_msg)
    stop(clean_msg, call. = FALSE)
  })
}
