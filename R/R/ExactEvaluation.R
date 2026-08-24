########################################################################################################################################
########################################################################################################################################
# ExactEvaluation
########################################################################################################################################
########################################################################################################################################
#' @title
#' Computing function mean values on linear extensions
#'
#' @description
#' `ExactEvaluation` computes the mean values of the input functions (defined on linear orders) over
#' the set of linear extensions of the input poset.
#' The linear extensions are generated according to the algorithm given in Habib M, Medina R, Nourine L and Steiner G (2001).
#'
#' @param poset An object of S4 class `POSet` representing the poset from which linear extensions are generated.
#' Argument `poset` must be created by using any function contained in the package aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...) .
#'
#' @param output_every_sec integer specifying a time interval (in seconds).
#' By specifying this argument, during the execution of `BubleyDyerEvaluation`, the number of linear extensions
#' progressively generated is printed on the R-Console, every `output_every_sec`seconds.
#'
#' @param f1 the function whose average value has to be computed.
#' `f1` can be an R-function having as a single argument a linear extension of `poset` and returning a numerical matrix.
#'
#' @param ... Further functions whose averages are to be computed.
#'
#' @return
#' A list of the computed averages, along with the number of linear extensions generated to compute them.
#'
#' @references  Habib M, Medina R, Nourine L and Steiner G (2001). Efficient algorithms on distributive lattices.
#' Discrete Applied Mathematics, 110, 169-187. https://doi.org/10.1016/S0166-218X(00)00258-4.
#'
#' @examples
#' el1 <- c("a", "b", "c", "d")
#' el2 <- c("x", "y")
#' el3 <- c("h", "k")
#' dom <- matrix(c(
#'   "a", "b",
#'   "c", "b",
#'   "b", "d"
#' ), ncol = 2, byrow = TRUE)
#'
#' pos1 <- POSet(elements = el1, dom = dom)
#'
#' pos2 <- LinearPOSet(elements = el2)
#'
#' pos3 <- LinearPOSet(elements = el3)
#'
#' pos <- ProductPOSet(pos1, pos2, pos3)
#'
#' # median_distr computes the frequency distribution of median profile
#'
#' elements <- POSetElements(pos)
#'
#' median_distr <- function(le) {
#'   n <- length(elements)
#'   if (n %% 2 != 0) {
#'     res <- (elements == le[(n + 1) / 2])
#'   } else {
#'     res <- (elements == le[n / 2])
#'   }
#'   res <- as.matrix(res)
#'   rownames(res) <- elements
#'   colnames(res) <- "median_distr"
#'   return (as.matrix(res))
#' }
#'\donttest{
#' res  <- ExactEvaluation(pos, output_every_sec=1, median_distr)
#'}
#'
#' @name ExactEvaluation
#' @export
ExactEvaluation <- function(poset, output_every_sec=NULL, f1, ...) {
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

  # Validate f1
  if (!is.function(f1)) {
    stop("'f1' must be an R function")
  }

  # Build function list and validate additional functions
  functions_list <- list(f1)
  additional_funcs <- list(...)
  if (length(additional_funcs) > 0) {
    for (i in seq_along(additional_funcs)) {
      if (!is.function(additional_funcs[[i]])) {
        stop(sprintf("Argument %d in '...' must be an R function", i + 1))
      }
      functions_list[[length(functions_list) + 1]] <- additional_funcs[[i]]
    }
  }

  # Call C++ implementation with error handling
  tryCatch({
    result <- .Call(C_ExactEvaluation, poset@ptr, output_every_sec, functions_list)
    return(result)

  }, error = function(err) {
    # Extract meaningful error message
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })
}
