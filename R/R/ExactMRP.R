#' @title
#' Computing Mutual Ranking Probabilities (MRP).
#'
#' @description
#' Computes the MRP matrix of a poset. The MRP associated to \eqn{a\leq b}, with \eqn{a} and \eqn{b} two elements of the input poset, is the share of linear extensions of it where \eqn{b} dominates \eqn{a}. The linear extensions are computed according to the algorithm given in Habib M, Medina R, Nourine L and Steiner G (2001).
#'
#' @param poset An object of S4 class `POSet` representing the poset whose MRP are computed.
#' Argument `poset` must be created by using any function contained in the package aimed at building object of S4 class `POSet`
#' (e.g. [POSet()], [LinearPOSet()], [ProductPOSet()], ...) .
#'
#' @param output_every_sec Integer specifying a time interval (in seconds).
#' By specifying this argument, during the execution of `ExactMRP` a message reporting the number of linear extensions progressively generated
#' is printed on the R-Console, every `output_every_sec` seconds.
#'
#' @return
#' A list of two elements: 1) the MRP matrix and 2) the number of linear extensions generated to compute it.
#'
#' @references Habib M, Medina R, Nourine L and Steiner G (2001). Efficient algorithms on distributive lattices.
#' Discrete Applied Mathematics, 110, 169-187. https://doi.org/10.1016/S0166-218X(00)00258-4.
#'
#' @examples
#'\donttest{
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
#' MRP <- ExactMRP(pos, output_every_sec=1)
#'}
#'
#' @name ExactMRP
#' @export
ExactMRP <- function(poset, output_every_sec=NULL) {
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

  # Call C++ implementation with error handling
  tryCatch({
    result <- .Call(C_ExactMRP, poset@ptr, output_every_sec)
    return(result)

  }, error = function(err) {
    # Extract meaningful error message
    err_msg <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    stop(trimws(err_parts[length(err_parts)]), call. = FALSE)
  })
}
