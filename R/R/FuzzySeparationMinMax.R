#' @title Fuzzy separation computation with minimum t-norm and maximum t-conorm
#'
#' @description
#' Starting from a poset dominance matrix, computes fuzzy separation matrices using minimum t-norm and maximum t-conorm.
#'
#' @param dom A square numeric matrix representing the dominance degree between pairs of poset elements.
#' Row and column names of `dom` are interpreted as the labels of the poset elements.
#' `dom` can be computed by using functions such as [BLSDominance()], [BubleyDyerMRP()], or [ExactMRP()].
#'
#' @param types Character vector specifying the types of fuzzy separation to be computed.
#' Possible choices are: \code{"symmetric"}, \code{"asymmetricLower"}, \code{"asymmetricUpper"},
#' \code{"vertical"}, and \code{"horizontal"}. Multiple choices can be provided
#' (e.g., \code{types = c("symmetric", "vertical")}). Defaults to \code{"symmetric"}.
#' For details on the definition of symmetric, asymmetric, vertical, and horizontal separations, see Fattore et al. (2024).
#'
#' @return
#' A list of the required fuzzy separation matrices.
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
#' FSep <- FuzzySeparationMinMax(BLS, types = c("symmetric", "vertical"))
#'
#' @name FuzzySeparationMinMax
#' @export
FuzzySeparationMinMax <- function(dom, types = c("symmetric", "asymmetricLower", "asymmetricUpper", "vertical", "horizontal")) {

  # Input validation
  if (!is.matrix(dom)) {
    stop("'dom' must be a square numeric matrix")
  }

  if (ncol(dom) != nrow(dom)) {
    stop("'dom' must be a square matrix")
  }

  # Elegant multiple validation with match.arg()
  # Default fallback if the user provides no arguments for 'types'
  if (missing(types)) {
    types <- "symmetric"
  } else {
    types <- match.arg(types, several.ok = TRUE)
  }

  CORE_TYPES  <- c("symmetric", "asymmetricLower", "asymmetricUpper")

  # Determine which core types are needed for computation
  functions_list <- types

  # "vertical" requires asymmetricLower and asymmetricUpper
  if ("vertical" %in% types) {
    functions_list <- c(functions_list, "asymmetricLower", "asymmetricUpper")
  }

  # "horizontal" requires symmetric, asymmetricLower, and asymmetricUpper
  if ("horizontal" %in% types) {
    functions_list <- c(functions_list, "symmetric", "asymmetricLower", "asymmetricUpper")
  }

  # Keep only core types (deduplicate and remove derived types)
  functions_list <- intersect(functions_list, CORE_TYPES)

  # Specify minimum t-norm and maximum t-conorm
  norm <- "minimum"
  conorm <- ""

  # Call C++ implementation with error handling
  tryCatch({
    # functions_list is passed as a simple character vector (STRSXP)
    result <- .Call(C_FuzzySeparation, dom, norm, conorm, functions_list)

    # Compute derived separation matrices
    if ("vertical" %in% types) {
      result[["vertical"]] <- abs(result[["asymmetricLower"]] - result[["asymmetricUpper"]])
    }

    if ("horizontal" %in% types) {
      result[["horizontal"]] <- result[["symmetric"]] - abs(result[["asymmetricLower"]] - result[["asymmetricUpper"]])
    }

    # Return only the requested types in the correct order
    result <- result[types]
    return(result)

  }, error = function(err) {
    # Extract meaningful error message
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })
}
