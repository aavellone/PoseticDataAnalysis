#' @title Fuzzy in-betweenness array computation
#'
#' @description
#' Starting from a poset dominance matrix, computes in-betweenness arrays using a user-supplied t-norm and t-conorm.
#'
#' @param dom A square numeric matrix representing the dominance degree between pairs of poset elements.
#' Row and column names of `dom` are interpreted as the labels of the poset elements.
#' `dom` can be computed by using functions such as [BLSDominance()], [BubleyDyerMRP()], or [ExactMRP()].
#'
#' @param norm An R function defining the t-norm.
#'
#' @param conorm An R function defining the t-conorm.
#'
#' @param types Character vector specifying the types of in-betweenness to be computed.
#' Possible choices are: \code{"symmetric"}, \code{"asymmetricLower"}, and \code{"asymmetricUpper"}.
#' Multiple choices can be provided (e.g., \code{types = c("symmetric", "asymmetricLower")}).
#' Defaults to \code{"symmetric"}.
#' For details on the definition of symmetric and asymmetric in-betweenness, see Fattore et al. (2024).
#'
#' @return
#' A list of three-dimensional arrays, one for each type of in-betweenness selected by `types`.
#' The array element at position \eqn{[i,j,k]} represents \eqn{finb_{p_i,p_j,p_k}} for symmetric in-betweenness,
#' \eqn{finb_{p_i<p_j<p_k}} for asymmetricLower in-betweenness, and \eqn{finb_{p_k<p_j<p_i}} for asymmetricUpper in-betweenness.
#'
#' @references Fattore, M., De Capitani, L., Avellone, A., and Suardi, A. (2024).
#' A fuzzy posetic toolbox for multi-criteria evaluation on ordinal data systems.
#' Annals of Operations Research, https://doi.org/10.1007/s10479-024-06352-3.
#'
#' @examples
#' el <- c("a", "b", "c", "d")
#'
#' dom_list <- matrix(c(
#'   "a", "b",
#'   "c", "b",
#'   "b", "d"
#' ), ncol = 2, byrow = TRUE)
#'
#' pos <- POSet(elements = el, dom = dom_list)
#'
#' BLS <- BLSDominance(pos)
#'
#' tnorm   <- function(x, y) { x * y }
#'
#' tconorm <- function(x, y) { x + y - x * y }
#'
#' FinB <- FuzzyInBetweenness(BLS, norm = tnorm, conorm = tconorm,
#'                             types = c("symmetric", "asymmetricLower"))
#'
#' @name FuzzyInBetweenness
#' @export
FuzzyInBetweenness <- function(dom, norm, conorm,
                               types = c("symmetric", "asymmetricLower", "asymmetricUpper")) {

  # Input validation
  if (!is.matrix(dom)) {
    stop("'dom' must be a square numeric matrix")
  }

  if (ncol(dom) != nrow(dom)) {
    stop("'dom' must be a square matrix")
  }

  if (!is.function(norm)) {
    stop("'norm' must be an R function")
  }

  if (!is.function(conorm)) {
    stop("'conorm' must be an R function")
  }

  # Default fallback if the user provides no arguments for 'types'
  if (missing(types)) {
    types <- "symmetric"
  } else {
    types <- match.arg(types, several.ok = TRUE)
  }

  # Ensure uniqueness of requested types (just in case the user writes types = c("symmetric", "symmetric"))
  functions_list <- unique(types)

  # Call C++ implementation with error handling
  tryCatch({
    # Passed as a simple string vector readable via STRING_ELT
    result <- .Call(C_FuzzyInBetweenness, dom, norm, conorm, functions_list)

    # Return only the requested types in the correct order
    result <- result[functions_list]
    return(result)

  }, error = function(err) {
    # Extract meaningful error message
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })
}
