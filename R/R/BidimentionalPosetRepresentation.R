#' @title
#' Bidimensional representation of multidimensional ordinal binary data
#'
#' @description
#' Finds an optimal 2D representation of profiles using a pair of lexicographic
#' linear extensions and BLS structural separation.
#'
#' @details
#' The \code{lpomStrategy} parameter allows choosing between two approximation
#' formulas derived from Brüggemann et al. (2003):
#' \itemize{
#'   \item \code{"absolute"}: Equation 13' in the reference (default).
#'   \item \code{"relative"}: Equation 13 in the reference.
#' }
#'
#' @param profile Boolean matrix of unique observed profiles.
#' @param weights Numeric vector of frequencies/weights.
#' @param variablesPriority Integer vector (permutation of 1:k).
#' @param lpomStrategy Character string: \code{"absolute"} (default) or \code{"relative"}.
#'
#' @return a list of 2 elements named `LossVAlue` and `Representation`.
#'
#' `LossVAlue` real number indicating the value of the global error \eqn{L(D^{out}|D^{inp}, p)} corresponding to the representation induced by the chosen `variablesPriority`.
#'
#' `Representation`  a data frame with \eqn{m} values (one value for each observed profile) of 5 variables named `profiles`, `x`, `y`, `weights` and `error`.
#' `$profile` is an integer vector containing the base-10 representation of the \eqn{k}-dimensional Boolean vectors representing observed profiles.
#' `$x` is an integer vector containing the x-coordinates of points representing observed profiles in the bidimensional representation.
#' `$y` is an integer vector containing the y-coordinates of points representing observed profiles in the bidimensional representation.
#' `$weights` is a real vector with the frequencies/weights of each observed profile.
#' `$error` is a real vector with the values of the approximation errors \eqn{L(b|D^{inp}, p)} associated to each observed profile
#' in the bidimensional representation.
#' @references Brueggemann R. et al. (2003). NERI Technical Report, No. 479.
#'
#' @examples
#' #SIMULATING OBSERVED BINARY DATA
#' #number of binary variables
#' k <- 6
#' #building observed profiles matrix
#' profiles <- sapply((0:(2^k-1)) ,function(x){ as.integer(intToBits(x))})
#' profiles <- t(profiles[1:k, ])
#' #building the vector of observation frequencies
#' weights <- sample.int(100, nrow(profiles), replace=TRUE)
#' #Chosing (at random) a variable priority
#' vp <- sample.int(k, k, replace=FALSE)
#' result <- BidimentionalPosetRepresentation(profiles, weights, vp)
#'#'
#' @export
BidimentionalPosetRepresentation <- function(profile, weights, variablesPriority,
                                             lpomStrategy = c("absolute", "relative")) {

  if (!is.matrix(profile)) stop("'profile' must be a matrix")

  k_vars <- ncol(profile)
  m_profiles <- nrow(profile)

  lpomStrategy <- match.arg(lpomStrategy)
  lpom_int <- if (lpomStrategy == "absolute") 0L else 1L

  if (any(weights == 0)) {
    keep_idx <- which(weights != 0)
    profile <- profile[keep_idx, , drop = FALSE]
    weights <- weights[keep_idx]
  }

  storage.mode(profile) <- "integer"
  storage.mode(weights) <- "numeric"
  storage.mode(variablesPriority) <- "integer"
  loss <- "LB"

  tryCatch({
    result <- .Call(C_RunBidimentionalPosetRepresentation,
                    profile,
                    weights,
                    loss,
                    lpom_int,
                    variablesPriority)

    result[["representation"]] <- as.data.frame(result[["representation"]],
                                                stringsAsFactors = FALSE)
    return(result)

  }, error = function(err) {
    stop(sprintf("Internal C++ Error: %s", conditionMessage(err)), call. = FALSE)
  })
}
