#' @title
#' Approximated separation matrices computation, using the Bubley-Dyer procedure (see Bubley and Dyer, 1999).
#'
#' @description
#' Computes approximated separation matrices, starting from a set of linear extensions
#' sampled according to the Bubley-Dyer procedure.
#'
#' @param generator The approximated separation matrices generator created by function [BuildBubleyDyerSeparationGenerator()].
#'
#' @param n number of linear extensions generated to compute the approximated MRP matrix.
#' See documentation for function [BubleyDyerMRP()] for further information on this argument.
#'
#' @param error A real number in \eqn{(0,1)} representing the "distance" from uniformity of the
#' sampling distribution of the linear extensions.
#' This parameter is used to determine the number of linear extensions to be sampled.
#' If both arguments `n` and `error` are specified by the user, the number of linear extensions actually generated
#' is `n`.
#' See documentation for function `BubleyDyerMRP` for further information on this argument.
#'
#' @param output_every_sec Integer specifying a time interval (in seconds).
#' By specifying this argument, during the execution of `BubleyDyerSeparation`, a message reporting the number of linear extensions
#' progressively generated is printed on the R-Console, every `output_every_sec` seconds.
#'
#' @details See the documentation of [BuildBubleyDyerSeparationGenerator()] for details on how the different types of separations
#' are defined and computed.
#'
#' @return A list containing: 1) the required type of approximated separation matrices, according to the parameter `type` used
#' to build the `generator` ( see[BuildBubleyDyerSeparationGenerator()]); 2) the number of generated linear extensions.
#'
#' @references Bubley, R., Dyer, M. (1999). Faster random generation of linear extensions.
#' Discrete Mathematics, 201, 81-88. https://doi.org/10.1016/S0012-365X(98)00333-1
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
#' BDgen <- BuildBubleyDyerSeparationGenerator(pos, seed = NULL,
#'                   types=c("symmetric", "asymmetricUpper", "vertical"))
#'
#' SEP_matrices <- BubleyDyerSeparation(BDgen, n=10000, output_every_sec=5)
#'
#' @name BubleyDyerSeparation
#' @export
BubleyDyerSeparation <- function(generator, n = NULL, error=NULL, output_every_sec=NULL) {
  # Input validation
  if (!methods::is(generator, "BubleyDyerSeparationGenerator")) {
    stop("'generator' must be an object of S4 class 'BubleyDyerSeparationGenerator'")
  }

  # Check that at least one of n or error is provided
  if (is.null(n) && is.null(error)) {
    stop("Either 'n' or 'error' must be specified")
  }

  if (!is.null(n)) {
    if (!is.numeric(n) || length(n) != 1 || n < 1 || n != round(n)) {
      stop("'n' must be a positive integer")
    }
    n <- as.integer(n)
  }

  if (!is.null(error)) {
    if (!is.numeric(error) || length(error) != 1 || error <= 0 || error >= 1) {
      stop("'error' must be a single numeric value in (0, 1)")
    }
  }

  if (!is.null(output_every_sec)) {
    if (!is.numeric(output_every_sec) || length(output_every_sec) != 1 ||
        output_every_sec < 1 || output_every_sec != round(output_every_sec)) {
      stop("'output_every_sec' must be a positive integer")
    }
    output_every_sec <- as.integer(output_every_sec)
  }

  # Call C++ implementation with error handling
  tryCatch({
    # SeparationTypes: "symmetric", "asymmetricLower", "asymmetricUpper", "vertical", "horizontal"
    # SeparationTypesC:"symmetric", "asymmetricLower", "asymmetricUpper"

    result <- .Call(C_BubleyDyerSeparation, generator@ptr, n, error, output_every_sec)

    # Compute derived separation matrices based on requested types
    if ("vertical" %in% generator@types) {
      vertical <- abs(result[["asymmetricLower"]] - result[["asymmetricUpper"]])
      result[["vertical"]] <- vertical
    }

    if ("horizontal" %in% generator@types) {
      horizontal <- result[["symmetric"]] - abs(result[["asymmetricLower"]] - result[["asymmetricUpper"]])
      result[["horizontal"]] <- horizontal
    }

    # Return only the requested separation types in the correct order
    result <- result[unlist(generator@types)]
    return(result)

  }, error = function(err) {
    # Extract meaningful error message
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })

}
