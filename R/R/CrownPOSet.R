#' @title
#' Building crowns.
#'
#' @description
#' Builds a crown from two unordered collections of elements, with the same size.
#'
#' @param elements_1 A list of character strings.
#'
#' @param elements_2 A list of character strings.
#'
#' @return
#' A crown, an object of S4 class `POSet`.
#'
#' @details
#' Let \eqn{a_1,\ldots,a_n} and \eqn{b_1,\ldots,b_n} be two disjoint collections of \eqn{n} elements.
#' The "crown" over them is the poset \eqn{P=(V,\lhd)} having \eqn{a_1,\ldots,a_n,b_1,\ldots,b_n}
#' as ground set and where \eqn{(a_i||a_j)}, \eqn{(b_i||b_j)}, \eqn{(a_i||b_i)} and \eqn{a_i\lhd b_j},
#' for each \eqn{i\neq j} (\eqn{||} stands for "incomparable to").
#'
#' @examples
#' elems1<-c("a1", "a2", "a3", "a4", "a5")
#' elems2<-c("b1", "b2", "b3", "b4", "b5")

#' crown<-CrownPOSet(elems1, elems2)
#'
#' @name CrownPOSet
#' @export
CrownPOSet <- function(elements_1, elements_2) {
  # Input validation
  if (!is.character(elements_1)) {
    stop("'elements_1' must be a character vector")
  }

  if (anyDuplicated(elements_1)) {
    stop("'elements_1' contains duplicated values")
  }

  if (!is.character(elements_2)) {
    stop("'elements_2' must be a character vector")
  }

  if (anyDuplicated(elements_2)) {
    stop("'elements_2' contains duplicated values")
  }

  if (length(elements_1) != length(elements_2)) {
    stop("'elements_1' and 'elements_2' must have the same length")
  }

  # Call C++ implementation with error handling
  tryCatch({
    ptr <- .Call(C_BuildCrownPOSet, elements_1, elements_2)
    result <- methods::new("POSet", ptr = ptr)
    return(result)

  }, error = function(err) {
    # Extract meaningful error message
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })
}
