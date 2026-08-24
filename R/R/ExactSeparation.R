#' @title Exact separation matrices computation
#'
#' @description
#' Computes exact separation matrices by evaluating the average separation over all the linear extensions of the input poset.
#' The linear extensions are generated according to the algorithm given in Habib et al. (2001).
#'
#' @param poset Object of S4 class `POSet` representing the poset whose separation matrix is computed.
#' Argument `poset` must be created by using any function contained in the package aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...).
#'
#' @param output_every_sec Integer specifying a time interval (in seconds).
#' By specifying this argument, during the execution of `ExactSeparation`, a message reporting the number of linear extensions
#' progressively generated is printed on the R-Console, every `output_every_sec` seconds.
#'
#' @param types Character vector specifying the types of separation to be computed.
#' Possible choices are: \code{"symmetric"}, \code{"asymmetricLower"}, \code{"asymmetricUpper"},
#' \code{"vertical"}, and \code{"horizontal"}. Multiple choices can be provided
#' (e.g., \code{types = c("symmetric", "vertical")}). Defaults to \code{"symmetric"}.
#'
#' @details The symmetric separation associated to two elements \eqn{a} and \eqn{b} of the input poset is the average
#' absolute difference between the positions of \eqn{a} and \eqn{b} observed over all linear extensions (whose elements are arranged in ascending order):
#'
#' \eqn{Sep_{ab}=\frac{1}{n}\sum_{i=1}^{n}|Pos_{l_i}(a)-Pos_{l_i}(b)|},
#'
#' where \eqn{n} is the number of linear extensions of the input poset;
#' \eqn{l_i} represents a single linear extension and \eqn{Pos_{l_i}(\cdot)} stands for the position of element \eqn{\cdot}
#' into the sequence of poset elements arranged in increasing order according to \eqn{l_i}.
#'
#' Asymmetric lower and upper separations are defined as:
#' \eqn{Sep_{a < b}=\frac{1}{n}\sum_{i=1}^{n}(Pos_{l_i}(b)-Pos_{l_i}(a))\mathbb{1}(a <_{l_i} b)},
#' \eqn{Sep_{b < a}=\frac{1}{n}\sum_{i=1}^{n}(Pos_{l_i}(a)-Pos_{l_i}(b))\mathbb{1}(b <_{l_i} a)},
#' where \eqn{a\leq_{l_i} b} means that \eqn{a} is lower or equal to \eqn{b} in the linear order defined by linear
#' extension \eqn{l_i} and \eqn{\mathbb{1}} is the indicator function. Note that \eqn{Sep_{ab}=Sep_{a < b}+Sep_{b < a}}.
#'
#' Vertical and horizontal separations (\eqn{vSep} and \eqn{hSep}, respectively) are defined as
#'
#' \eqn{vSep_{ab}=|Sep_{a < b}-Sep_{b < a}|} and \eqn{hSep_{ab}=Sep_{ab}-vSep_{ab}}.
#'
#' For a detailed explanation on why \eqn{vSep} and \eqn{hSep} can be interpreted as vertical and horizontal components
#' of the separation between poset elements, see Fattore et al. (2024).
#'
#' @return A list containing: 1) the required exact separation matrices, according to the `types` parameter;
#' 2) the number of linear extensions generated to compute them.
#'
#' @references Habib, M., Medina, R., Nourine, L., and Steiner, G. (2001). Efficient algorithms on distributive lattices.
#' Discrete Applied Mathematics, 110, 169-187. https://doi.org/10.1016/S0166-218X(00)00258-4.
#'
#' Fattore, M., De Capitani, L., Avellone, A., and Suardi, A. (2024).
#' A fuzzy posetic toolbox for multi-criteria evaluation on ordinal data systems.
#' Annals of Operations Research, https://doi.org/10.1007/s10479-024-06352-3.
#'
#'
#' @examples
#' el <- c("a", "b", "c", "d")
#'
#' dom <- matrix(c(
#'   "a", "b",
#'   "c", "b",
#'   "b", "d"
#' ), ncol = 2, byrow = TRUE)
#'
#' pos <- POSet(elements = el, dom = dom)
#'
#' SEP_matrices <- ExactSeparation(pos, output_every_sec = 5,
#'                                 types = c("symmetric", "asymmetricUpper", "vertical"))
#'
#' @name ExactSeparation
#' @export
ExactSeparation <- function(poset, output_every_sec = NULL,
                            types = c("symmetric", "asymmetricLower", "asymmetricUpper", "vertical", "horizontal")) {

  # Input validation
  if (!methods::is(poset, "POSet")) {
    stop("'poset' must be an object of S4 class 'POSet'")
  }

  if (!is.null(output_every_sec)) {
    if (!is.numeric(output_every_sec) || length(output_every_sec) != 1 ||
        output_every_sec < 1 || output_every_sec != round(output_every_sec)) {
      stop("'output_every_sec' must be a positive integer")
    }
    output_every_sec <- as.integer(output_every_sec)
  }

  # Elegant multiple validation with match.arg()
  # Default fallback if the user provides no arguments for 'types'
  if (missing(types)) {
    types <- "symmetric"
  } else {
    types <- match.arg(types, several.ok = TRUE)
  }

  CORE_TYPES  <- c("symmetric", "asymmetricLower", "asymmetricUpper")

  # Determine which core types are needed for computation (building a flat character vector)
  functions_list <- types

  # "vertical" requires asymmetricLower and asymmetricUpper
  if ("vertical" %in% types) {
    functions_list <- c(functions_list, "asymmetricLower", "asymmetricUpper")
  }

  # "horizontal" requires symmetric, asymmetricLower, and asymmetricUpper
  if ("horizontal" %in% types) {
    functions_list <- c(functions_list, "symmetric", "asymmetricLower", "asymmetricUpper")
  }

  # Keep only core types (intersect automatically deduplicates and returns a character vector)
  functions_list <- intersect(functions_list, CORE_TYPES)

  # Call C++ implementation with error handling
  tryCatch({
    # functions_list is now a pure STRSXP (character vector), safe for STRING_ELT
    result <- .Call(C_ExactSeparation, poset@ptr, output_every_sec, functions_list)

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
