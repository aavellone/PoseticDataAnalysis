#' @title
#' Generator for the approximated computation of the mean value of functions over linear
#' extensions.
#'
#' @description
#' `BuildBubleyDyerEvaluationGenerator`
#' creates an object of S4 class `BubleyDyerEvaluationGenerator`, for the estimation of
#' the mean values of the input functions, over linear extensions sampled according to the Bubley-Dyer procedure.
#' Actually, this function does not perform the computation of mean values, but just generates the object that will compute them
#' by using function `BubleyDyerEvaluation`.
#'
#' @param poset An object of S4 class `POSet` representing the poset from which linear extensions are generated.
#' Object `poset` must be created by using any function contained in the package aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...) .
#'
#' @param seed Positive integer to initialize random linear extension generation. Set `seed=NULL` for random initialization.
#' The generator seeds are 64-bit: to use a value above 2^53 (which R cannot
#' represent exactly) pass it as a string of decimal digits, e.g. `"12345678901234567890"`.
#'
#' @param f1 The function whose mean value is to be computed.
#' `f1` must be an R-function having as a single parameter a linear extension of `poset` and returning a numerical matrix.
#'
#' @param ... Further functions whose mean values are to be computed.
#'
#' @return
#' An object of S4-class `BubleyDyerEvaluationGenerator`.
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
#'
#' BDgen <- BuildBubleyDyerEvaluationGenerator(poset = pos, seed = NULL, median_distr)
#'
#' @name BuildBubleyDyerEvaluationGenerator
#' @export
BuildBubleyDyerEvaluationGenerator <- function(poset, seed, f1, ...)  {
  # Input validation
  if (!methods::is(poset, "POSet")) {
    stop("'poset' must be an object of S4 class 'POSet'")
  }

  # Il seme viaggia verso C++ come stringa di cifre decimali (i semi sono a
  # 64 bit, fuori dalla portata degli interi di R). NULL = seme casuale.
  if (!is.null(seed) && length(seed) != 1L) {
    stop("'seed' must be a single value")
  }
  seed <- .seed_to_character(seed)

  # Validate f1
  if (!is.function(f1)) {
    stop("'f1' must be an R function")
  }

  # Build function list starting with f1
  functions_list <- list(f1)

  # Validate and add additional functions
  additional_funcs <- list(...)
  if (length(additional_funcs) > 0) {
    for (i in seq_along(additional_funcs)) {
      func <- additional_funcs[[i]]
      if (!is.function(func)) {
        stop(sprintf("Argument %d in '...' must be an R function", i + 1))
      }
      functions_list[[length(functions_list) + 1]] <- func
    }
  }

  # Call C++ implementation with error handling
  tryCatch({
    # Il backend restituisce list(ptr, seed): vedi LEBubleyDyer.
    res <- .Call(C_BuildBubleyDyerEvaluationGenerator, poset@ptr, seed, functions_list)
    result <- methods::new("BubleyDyerEvaluationGenerator", ptr = res$ptr, seed = res$seed)
    return(result)

  }, error = function(err) {
    # Extract meaningful error message
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })

}
