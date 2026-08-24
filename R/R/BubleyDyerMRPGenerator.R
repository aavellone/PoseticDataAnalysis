########################################################################################################################################
########################################################################################################################################
# BubleyDyerMRPGenerator
########################################################################################################################################
########################################################################################################################################

#' @title
#' Generator of an approximated MRP matrix.
#'
#' @description
#' Creates an object of S4 class `BubleyDyerMRPGenerator` for the computation of an approximated
#' MRP matrix, starting from a set of random linear extensions, sampled according to the Bubley-Dyer procedure.
#' Actually, this function does not compute the MRP matrix, but just the object that will compute it,
#' by using function `BubleyDyerMRP`.
#'
#' @param poset Object of S4 class `POSet` representing the poset whose MRP are computed.
#' Argument `poset` must be created by using any function contained in the package aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...) .
#'
#' @param seed Positive integer to initialize the random linear extension generation.
#' The generator seeds are 64-bit: to use a value above 2^53 (which R cannot
#' represent exactly) pass it as a string of decimal digits, e.g. `"12345678901234567890"`.
#'
#' @return
#' An object of S4 class `BubleyDyerMRPGenerator`.
#' The returned object also carries the 64-bit seed actually used, as a
#' string of decimal digits (slot `seed`): passing it back to `seed`
#' reproduces exactly the same sequence.
#'
#' @examples
##################
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
#' BDgen <- BubleyDyerMRPGenerator(pos)
#'
#' @name BubleyDyerMRPGenerator
#' @export
BubleyDyerMRPGenerator <- function(poset=NULL, seed=NULL) {
  # Input validation
  if (is.null(poset)) {
    stop("'poset' must be provided")
  }

  if (!methods::is(poset, "POSet")) {
    stop("'poset' must be an object of S4 class 'POSet'")
  }

  # Il seme viaggia verso C++ come stringa di cifre decimali (i semi sono a
  # 64 bit, fuori dalla portata degli interi di R). NULL = seme casuale.
  if (!is.null(seed) && length(seed) != 1L) {
    stop("'seed' must be a single value")
  }
  seed <- .seed_to_character(seed)

  # Call C++ implementation with error handling
  tryCatch({
    # Il backend restituisce list(ptr, seed): vedi LEBubleyDyer.
    res <- .Call(C_BuildBubleyDyerMRPGenerator, poset@ptr, seed)
    result <- methods::new("BubleyDyerMRPGenerator", ptr = res$ptr, seed = res$seed)
    return(result)

  }, error = function(err) {
    # Extract meaningful error message
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })

}
