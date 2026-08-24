#' @title
#' Constructing a component-wise poset of binary vectors.
#'
#' @description
#' Constructs a component-wise poset, starting from a collection of binary variables.
#'
#' @param variables A vector of character strings (the names of the input binary variables).
#'
#' @return
#' An object of S4 class `BinarVariablePOSet` (subclass of `POSet`).
#'
#' @details
#' Given \eqn{k} input binary variables, the function produces a poset \eqn{(V,\leq_{cmp})}, where \eqn{V} is the set of \eqn{2^k} binary vectors built from the variables, and \eqn{\leq_{cmp}} is the component-wise order.
#'
#' @examples
#' vrbs <- c("var1", "var2", "var3")
#' binPoset <-  BinaryVariablePOSet(variables = vrbs)
#'
#' @name BinaryVariablePOSet
#' @export
BinaryVariablePOSet <- function(variables) {
  # Input validation
  if (!is.character(variables)) {
    stop("'variables' must be a character vector")
  }

  if (length(variables) == 0) {
    stop("'variables' must contain at least one element")
  }

  # Check for duplicates
  if (anyDuplicated(variables)) {
    stop("'variables' contains duplicated values")
  }

  # Call C++ implementation with error handling
  tryCatch({
    ptr <- .Call(C_BuildBinaryVariablePOSet, variables)
    result <- methods::new("BinaryVariablePOSet", ptr = ptr)
    return(result)

  }, error = function(err) {
    # Extract meaningful error message
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })
}
