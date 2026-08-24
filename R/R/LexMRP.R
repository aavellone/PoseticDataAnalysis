#' @title
#' MRP matrix computation over the set of lexicographic linear extensions.
#'
#' @description Considering the component-wise poset built starting from \eqn{k} ordinal variables,
#'  computes MRP matrix by analyzing all poset lexicographic linear extensions.
#'
#' @param nvar positive integer specifying the number \eqn{k} of ordinal variables.
#'
#' @param deg parameter specifying the number of degrees of each variable. If all \eqn{k} variables have the same number \eqn{m} of degrees, it can be:
#' 1) the positive integer \eqn{m}. In this case variable degree labels are supposed to be the integers \eqn{1<2<...<m}
#' and columns and rows of the computed MRP matrix are named accordingly to this;
#' 2) a character vector of length \eqn{m} specifying the variable degree labels (in this case columns and rows of the computed MRP matrix are named accordingly to \code{deg}).
#'
#' If the \eqn{k} variables have different number \eqn{(m_1,...,m_k)} of degrees, it can be:
#' 1) a length-\eqn{k} positive integers vector specifying the values of \eqn{m_1,...,m_k}. In this case variable degree labels for the \eqn{j}-th variable are
#' supposed to be the integers \eqn{1<2<...<m_j} and columns and rows of the computed MRP matrix are named accordingly to this;
#' 2) a list of \eqn{k} character vectors. The \eqn{j}-th list element is a character vector of length \eqn{m_j} specifying the degree labels for the \eqn{j}-th variable
#' (in this case columns and rows of the computed MRP matrix are named accordingly to \code{deg}).
#'
#'
#' @return the MRP matrix computed over the set of lexicographic linear extensions.
#'
#' @examples
#'
#' #variables with common number of degrees
#' # default labels for variable degrees
#' nvar <- 3
#' deg  <- 4
#' lMRP <- LexMRP(nvar=nvar, deg=deg)
#'
#' #user supplied variable degree labels
#' nvar <- 3
#' deg  <- c("a","b","c","d")
#' lMRP <- LexMRP(nvar=nvar, deg=deg)
#'
#'
#' #variables with different numbers of degrees
#' # default labels for variable degrees
#' nvar <- 3
#' deg  <- c(4,2,3)
#' lMRP <- LexMRP(nvar=nvar, deg=deg)
#'
#' #user supplied variable degree labels
#' nvar <- 3
#' deg  <- list(c("a","b","c","d"),c("0","1"),c("x","y","z"))
#' lMRP <- LexMRP(nvar=nvar, deg=deg)
#'
#' @name LexMRP
#' @export
LexMRP <- function(nvar, deg) {
  if (!is.numeric(nvar) || length(nvar) != 1L || !is.finite(nvar) || nvar <= 0 || nvar != round(nvar)) {
    stop("'nvar' must be a finite positive integer.", call. = FALSE)
  }

  nvar <- as.integer(nvar)

  valid_deg <- (is.numeric(deg) && length(deg) == 1L && is.finite(deg) && deg > 0 && deg == round(deg)) ||
    (is.character(deg)) ||
    (is.numeric(deg) && length(deg) == nvar && all(is.finite(deg)) && all(deg > 0) && all(deg == round(deg))) ||
    (is.list(deg) && length(deg) == nvar && all(sapply(deg, is.character)))


  if (!valid_deg) {
    stop(paste0(
      "'deg' must be a positive integer, a character vector, ",
      "a positive integer vector of length ", nvar,
      ", or a list of ", nvar, " character vectors."
    ), call. = FALSE)
  }

  # Normalise deg to a list of character vectors
  if (is.numeric(deg) && length(deg) == 1L) {
    deg <- replicate(nvar, as.character(seq_len(deg)), simplify = FALSE)
  } else if (is.character(deg)) {
    deg <- replicate(nvar, deg, simplify = FALSE)
  } else if (is.numeric(deg)) {
    deg <- lapply(deg, function(x) as.character(seq_len(x)))
  } else {
    deg <- as.list(deg)
  }

  tryCatch({
    return(.Call(C_RunLexMRP, deg))
  }, error = function(err) {
    err_msg   <- conditionMessage(err)
    err_parts <- strsplit(err_msg, ":", fixed = TRUE)[[1]]
    clean_msg <- if (length(err_parts) > 1L) trimws(err_parts[length(err_parts)]) else trimws(err_msg)
    stop(clean_msg, call. = FALSE)
  })
}
