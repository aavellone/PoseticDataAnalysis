#' @title Dimensionality reduction of multidimensional ordinal binary data
#'
#' @description Starting from a dataset with \eqn{n} statistical units, scored against \eqn{k}
#' ordinal 0/1-indicators and partially ordered component-wise into a Boolean lattice \eqn{B_k=(\{0,1\}^k,\leq_{cmp})},
#' it finds the bidimensional data representation that optimally preserves the input order relation.
#' The algorithm finding the best bidimensional representation is optimized by using a parallel C++ implementation.
#'
#' @details
#' The bidimensional representation is constructed by approximating the structural separation
#' based on the Local Partial Order Model (LPOM). The \code{lpomStrategy} parameter allows
#' choosing between two approximation formulas derived from Brüggemann et al. (2003):
#' \itemize{
#'   \item \code{"absolute"}: Computes the absolute structural dominance
#'   (referred to as \bold{Equation 13'} in the reference). This is the default strategy.
#'   \item \code{"relative"}: Computes the relative mutual dominance probability
#'   (referred to as \bold{Equation 13} in the reference), normalizing the structural scores
#'   into the \code{[0, 1]} interval.
#' }
#'
#' @param profile Boolean matrix of dimension \eqn{m\times k} of the unique \eqn{m\leq n} different observed profiles.
#' Each observed profile is a row of `profile`. Each observed profile is repeated only once in the matrix `profile`.
#'
#' @param weights real vector of length \eqn{m} with the frequencies/weights of each observed profiles.
#'
#' @param lpomStrategy Character string indicating the LPOM approximation strategy to use.
#' Must be one of \code{"absolute"} (default) or \code{"relative"}.
#'
#' @param output_every_sec Integer specifying a time interval (in seconds).
#' By specifying this argument, a message reporting the number of reversed pairs
#' of lexicographic linear extensions analyzed is printed on the R-Console.
#'
#' @param thread_share real number in the interval \eqn{(0,1]} specifying the share of CPU threads to be involved
#' in the algorithm execution.
#'
#' @return a list of 5 elements named `allLoss`, `variablesPriority`, `bestLossValue`, `bestVariablePriority`, and `bestRepresentation`.
#' (See details in the standard documentation for the output structure).
#'
#' @references Brueggemann R., Lerche D. B., Sørensen P. B. (2003). First attempts to relate structures of
#' Hasse diagrams with mutual probabilities. National Environmental Research Institute,
#' Denmark - NERI Technical Report, No. 479.
#'
#' @examples
#' #SIMULATING OBSERVED BINARY DATA
#' k <- 6
#' profiles <- sapply((0:(2^k-1)) ,function(x){ as.integer(intToBits(x))})
#' profiles <- t(profiles[1:k, ])
#' weights <- sample.int(100, nrow(profiles), replace=TRUE)
#'
#' #FINDING THE OPTIMAL BIDIMENSIONAL REPRESENTATION
#' result <- OptimalBidimensionalEmbedding(profiles, weights, lpomStrategy = "absolute")
#'
#' @name OptimalBidimensionalEmbedding
#' @export
OptimalBidimensionalEmbedding <- function(profile, weights, lpomStrategy = c("absolute", "relative"), output_every_sec=NULL, thread_share=1.0) {
  loss <- "LB"

  if (!is.matrix(profile) || !all(profile %in% c(0, 1))) {
    stop("'profile' must be a 0/1 matrix.", call. = FALSE)
  }
  mode(profile) <- "integer"

  if (!is.numeric(weights) || !is.vector(weights)) {
    stop("'weights' must be a numeric vector.", call. = FALSE)
  }
  if (length(weights) != nrow(profile)) {
    stop("'weights' length must equal the number of rows of 'profile'.", call. = FALSE)
  }
  mode(weights) <- "numeric"

  # --- VALIDAZIONE E CONVERSIONE LPOM STRATEGY ---
  lpomStrategy <- match.arg(lpomStrategy)
  lpom_int <- if (lpomStrategy == "absolute") 0L else 1L

  # Remove zero-weight profiles
  zero_idx <- which(weights == 0)
  if (length(zero_idx) > 0) {
    profile <- profile[-zero_idx, , drop = FALSE]
    weights <- weights[-zero_idx]
  }

  if (!is.null(output_every_sec)) {
    if (!is.numeric(output_every_sec) || length(output_every_sec) != 1L ||
        !is.finite(output_every_sec) || output_every_sec <= 0 || output_every_sec != round(output_every_sec)) {
      stop("'output_every_sec' must be a finite positive integer.", call. = FALSE)
    }
    output_every_sec <- as.integer(output_every_sec)
  }

  if (!is.numeric(thread_share) || length(thread_share) != 1L ||
      !is.finite(thread_share) || thread_share <= 0 || thread_share > 1.0) {
    stop("'thread_share' must be a single finite number in the interval (0, 1].", call. = FALSE)
  }

  tryCatch({
    # ATTENZIONE: Passaggio di lpom_int e aggiunta di PACKAGE
    result <- .Call(C_RunDimensionalityReduction,
                    profile,
                    weights,
                    loss,
                    lpom_int,
                    output_every_sec,
                    thread_share)

    result[["bestRepresentation"]] <- data.frame(result[["bestRepresentation"]])
    return(result)
  }, error = function(err) {
    err_msg   <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    clean_msg <- if (length(err_parts) > 1L) trimws(err_parts[length(err_parts)]) else trimws(err_msg)
    stop(clean_msg, call. = FALSE)
  })
}
