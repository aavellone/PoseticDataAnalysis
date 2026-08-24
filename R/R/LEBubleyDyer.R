#' @title
#' Generator of linear extensions through the Bubley-Dyer procedure.
#'
#' @description
#' Creates an object of S4 class `BubleyDyerGenerator`, needed to sample the linear extensions of a given poset according to the Bubley-Dyer
#' procedure. Actually, this function does not sample the linear extensions, but just generates the object that will sample them
#' by using function `LEGet`.
#'
#' @param poset An object of S4 class `POSet` representing the poset whose linear extensions are generated.
#' Argument `poset` must be created by using any function contained in the package aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...) .
#'
#' @param seed Positive integer to initialize the random linear extension generation. If `NULL` (default),
#' the seed is not set and results will differ across runs.
#' The generator seeds are 64-bit: to use a value above 2^53 (which R cannot
#' represent exactly) pass it as a string of decimal digits, e.g. `"12345678901234567890"`.
#'
#' @return
#' An object of S4 class `BubleyDyerGenerator`.
#' The returned object also carries the 64-bit seed actually used, as a
#' string of decimal digits (slot `seed`): passing it back to `seed`
#' reproduces exactly the same sequence.
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
#' LEgenBD <- LEBubleyDyer(pos)
#'
#' @name LEBubleyDyer
#' @export
LEBubleyDyer <- function(poset, seed=NULL) {
  if (!methods::is(poset, "POSet")) {
    stop("'poset' must be of class POSet.", call. = FALSE)
  }

  # Il seme viaggia verso C++ come stringa di cifre decimali (i semi sono a
  # 64 bit, fuori dalla portata degli interi di R). NULL = seme casuale.
  # NB: qui passava una list(), che il backend rifiutava come tipo non valido
  # (o, se vuota, interpretava silenziosamente come "nessun seme").
  if (!is.null(seed) && length(seed) != 1L) {
    stop("'seed' must be a single value.", call. = FALSE)
  }
  seed_arg <- .seed_to_character(seed)

  tryCatch({
    # Il backend restituisce list(ptr, seed): 'seed' e' il seme effettivamente
    # usato (quello fornito, oppure quello estratto a caso), come stringa.
    res <- .Call(C_BuildBubleyDyerLEGenerator, poset@ptr, seed_arg)
    if (is.null(res$ptr)) {
      stop("The C++ generator returned a NULL pointer.", call. = FALSE)
    }
    return(methods::new("BubleyDyerGenerator", ptr = res$ptr, seed = res$seed))
  }, error = function(err) {
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    clean_msg <- if (length(err_parts) > 1L) trimws(err_parts[length(err_parts)]) else trimws(err_msg)
    stop(clean_msg, call. = FALSE)
  })
}
