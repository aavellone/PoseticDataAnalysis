#' @title
#' Generator of an approximated separation matrix.
#'
#' @description
#' Creates an object of S4 class `BubleyDyerSeparationGenerator` for the computation of approximated
#' separation matrices, starting from a set of random linear extensions, sampled according to the Bubley-Dyer procedure (see Bubley and Dyer, 1999)
#' Actually, this function does not compute the separation matrices, but just the object that will compute them,
#' by using function `BubleyDyerSeparation`.
#'
#' @param poset Object of S4 class `POSet` representing the poset whose separation matrices are to be computed.
#' Argument `poset` must be created by using any function contained in the package aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...) .
#'
#' @param seed Positive integer to initialize the random linear extension generation.
#' The generator seeds are 64-bit: to use a value above 2^53 (which R cannot
#' represent exactly) pass it as a string of decimal digits, e.g. `"12345678901234567890"`.
#'
#' @param types Character vector specifying the types of separation to be computed.
#' Possible choices are: \code{"symmetric"}, \code{"asymmetricLower"}, \code{"asymmetricUpper"},
#' \code{"vertical"}, and \code{"horizontal"}. Multiple choices can be provided
#' (e.g., \code{types = c("symmetric", "vertical")}). Defaults to \code{"symmetric"}.
#'
#' @details The symmetric separation associated to elements \eqn{a} and \eqn{b} in the input poset is the average
#' absolute difference between the positions of \eqn{a} and \eqn{b} observed in the sampled linear extensions (whose elements are arranged in ascending order):
#'
#' \eqn{Sep_{ab}=\frac{1}{n}\sum_{i=1}^{n}|Pos_{l_i}(a)-Pos_{l_i}(b)|},
#'
#' where \eqn{n} is the numbers of sampled linear extensions;
#' \eqn{l_i} represents a sampled linear extension and \eqn{Pos_{l_i}(\cdot)} stands for the position of element \eqn{\cdot}
#' into the sequence of poset elements arranged in increasing order according to \eqn{l_i}.
#'
#' Asymmetric lower and upper separations are defined as:
#' \eqn{Sep_{a < b}=\frac{1}{n}\sum_{i=1}^{n}(Pos_{l_i}(b)-Pos_{l_i}(a))\mathbb{1}(a <_{l_i} b)},
#' \eqn{Sep_{b < a}=\frac{1}{n}\sum_{i=1}^{n}(Pos_{l_i}(a)-Pos_{l_i}(b))\mathbb{1}(b <_{l_i} a)},
#' where \eqn{a\leq_{l_i} b} means that \eqn{a} is lower or equal to \eqn{b} in the linear order defined by linear
#' extension \eqn{l_i} and \eqn{\mathbb{1}} is the indicator function. Note that \eqn{Sep_{ab}=Sep_{a < b}+Sep_{a < b}}.
#'
#' Vertical and horizontal separations (\eqn{vSep} and \eqn{hSep}, respectively) are defined as
#'
#' \eqn{vSep_{ab}=|Sep_{a < b}-Sep_{b < a}|} and \eqn{hSep_{ab}=Sep_{ab}-vSep_{ab}|}.
#'
#' For a detailed explanation on why \eqn{vSep} and \eqn{hSep} can be interpreted as vertical and horizontal components
#' of the separation between two poset elements, see Fattore et. al (2024).
#'
#' @return
#' An object of S4 class `BubleyDyerSeparationGenerator`.
#'
#' @references Bubley, R., Dyer, M. (1999). Faster random generation of linear extensions.
#' Discrete Mathematics, 201, 81-88. https://doi.org/10.1016/S0012-365X(98)00333-1
#'
#' Fattore, M., De Capitani, L., Avellone, A., and Suardi, A. (2024).
#' A fuzzy posetic toolbox for multi-criteria evaluation on ordinal data systems.
#' Annals of Operations Research, https://doi.org/10.1007/s10479-024-06352-3.
#'
#'
#' @examples
#' el <- c("a", "b", "c", "d")
#' dom <- matrix(c("a", "b", "c", "b", "b", "d"), ncol = 2, byrow = TRUE)
#' pos <- POSet(elements = el, dom = dom)
#'
#' # Computes multiple separation types
#' BDgen <- BuildBubleyDyerSeparationGenerator(pos, seed = NULL,
#'               types = c("symmetric", "asymmetricUpper", "vertical"))
#'
#' @name BuildBubleyDyerSeparationGenerator
#' @export
BuildBubleyDyerSeparationGenerator <- function(poset, seed,
                                               types = c("symmetric", "asymmetricLower",
                                                         "asymmetricUpper", "vertical", "horizontal")) {

  # 1. Validazione Input di Base
  if (!methods::is(poset, "POSet")) {
    stop("'poset' must be an object of S4 class 'POSet'")
  }

  # Il seme viaggia verso C++ come stringa di cifre decimali (i semi sono a
  # 64 bit, fuori dalla portata degli interi di R). NULL = seme casuale.
  if (!is.null(seed) && length(seed) != 1L) {
    stop("'seed' must be a single value")
  }
  seed <- .seed_to_character(seed)

  # 2. VALIDAZIONE ELEGANTE MULTIPLA CON match.arg()
  # Se l'utente non specifica nulla, per sicurezza calcoliamo solo "symmetric".
  # Altrimenti match.arg verificherà tutto l'array inserito dall'utente.
  if (missing(types)) {
    types <- "symmetric"
  } else {
    types <- match.arg(types, several.ok = TRUE)
  }

  # 3. Preparazione functions_list per il calcolo
  CORE_TYPES <- c("symmetric", "asymmetricLower", "asymmetricUpper")
  functions_list <- types

  # Se "vertical" è richiesto, servono le due asimmetriche
  if ("vertical" %in% types) {
    functions_list <- c(functions_list, "asymmetricLower", "asymmetricUpper")
  }

  # Se "horizontal" è richiesto, servono la simmetrica e le due asimmetriche
  if ("horizontal" %in% types) {
    functions_list <- c(functions_list, "symmetric", "asymmetricLower", "asymmetricUpper")
  }

  # intersect() scarta automaticamente i duplicati e isola solo i CORE_TYPES
  functions_list <- intersect(functions_list, CORE_TYPES)

  # 4. Chiamata al C++
  tryCatch({
    # Il backend restituisce list(ptr, seed): vedi LEBubleyDyer.
    res <- .Call(C_BuildBubleyDyerSeparationGenerator, poset@ptr, seed, functions_list)

    result <- methods::new("BubleyDyerSeparationGenerator", ptr = res$ptr,
                           types = as.list(types), seed = res$seed)
    return(result)

  }, error = function(err) {
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })
}
