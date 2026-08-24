#' @title
#' Computing the BLS dominance matrix of a poset
#'
#' @description
#' Computes the dominance matrix of the input poset based on the approximations
#' introduced by Brueggemann, Lerche, and Sørensen (2003). This measure estimates
#' the probability that one element dominates another without the computational burden
#' of generating all linear extensions.
#'
#' @details
#' The function supports two types of computations based on the original paper:
#' \itemize{
#'   \item \code{"absolute"}: Computes the absolute structural dominance
#'   (Equation 13' in the reference). This is the default strategy.
#'   \item \code{"relative"}: Computes the estimated mutual dominance probability
#'   (Equation 13 in the reference). It returns values in the \code{[0, 1]} interval.
#' }
#'
#' @param poset Object of S4 class `POSet`.
#' @param type Character string indicating the type of BLS dominance:
#' \code{"absolute"} (default) or \code{"relative"}.
#'
#' @return
#' A numeric square matrix with row and column names.
#'
#' @references Brueggemann R., Lerche D. B., Sørensen P. B. (2003). First attempts to relate structures of
#' Hasse diagrams with mutual probabilities. NERI Technical Report, No. 479.
#'
#' @examples
#'
#' el <- c("a", "b", "c", "d")
#' dom <- matrix(c(
#'   "a", "b",
#'   "c", "b",
#'   "b", "d"
#' ), ncol = 2, byrow = TRUE)
#' pos <- POSet(elements = el, dom = dom)
#'
#' res <- BLSDominance(pos)
#'
#'
#' @export
BLSDominance <- function(poset, type = c("absolute", "relative")) {
  if (!methods::is(poset, "POSet")) {
    stop("'poset' must be an object of S4 class 'POSet'")
  }

  type <- match.arg(type)
  type_int <- if (type == "absolute") 0L else 1L

  tryCatch({
    # Aggiunto PACKAGE = "poseticDataAnalysis"
    result <- .Call(C_BruggemannLercheSorensenDominance, poset@ptr, type_int)
    return(result)
  }, error = function(e) {
    stop(paste("Error during BLS Dominance computation in C++:", e$message))
  })
}

