#' @title Building fences
#'
#' @description
#' Builds a fence from an unordered collection of elements.
#'
#' @param elements A character vector (the names of the fence elements).
#'
#' @param orientation Either "upFirst" (the first element dominates the second) or
#' "downFirst" (the second element dominates the first). Default is "upFirst".
#'
#' @return
#' A fence, an object of S4 class `POSet`.
#'
#' @examples
#' elems <- c("a", "b", "c", "d", "e")
#' fence <- FencePOSet(elems, orientation = "upFirst")
#'
#' @name FencePOSet
#' @export
FencePOSet <- function(elements, orientation="upFirst") {
  # Input validation
  if (!is.character(elements)) {
    stop("'elements' must be a character vector")
  }

  if (anyDuplicated(elements)) {
    stop("'elements' contains duplicated values")
  }

  if (!is.character(orientation) || length(orientation) != 1 ||
      !(orientation %in% c("upFirst", "downFirst"))) {
    stop("'orientation' must be either \"upFirst\" or \"downFirst\"")
  }

  # Call C++ implementation with error handling
  tryCatch({
    ptr <- .Call(C_BuildFencePOSet, elements, orientation == "upFirst")
    result <- methods::new("POSet", ptr = ptr)
    return(result)

  }, error = function(err) {
    # Extract meaningful error message
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })
}
