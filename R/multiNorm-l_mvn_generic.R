#' Log of the Likelihood of the Multivariate Normal Distribution - Generic
#'
#' Calculates the log of the likelihood function
#' of the multivariate normal distribution
#' for the ith observation.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param x Numeric vector of length `k`.
#'   The ith vector of observations.
#' @param mu Numeric vector of length `k`.
#'   Mean vector.
#' @param sigmacap Numeric matrix of size `k` by `k`.
#'   Covariance matrix.
#'
#' @returns A vector.
#'
#' @examples
#' n <- 5
#' mu <- c(0, 0)
#' sigmacap <- matrix(
#'   data = c(
#'     1, 0.5, 0.5, 1
#'   ),
#'   nrow = 2
#' )
#'
#' xcap <- as.data.frame(
#'   t(
#'     rmvn_chol(
#'       n = n,
#'       mu = mu,
#'       sigmacap = sigmacap
#'     )
#'   )
#' )
#'
#' lapply(
#'   X = xcap,
#'   FUN = l_mvn_generic,
#'   mu = mu,
#'   sigmacap = sigmacap
#' )
#' @export
#' @family Multivariate Normal Distribution Functions
#' @keywords multiNorm
l_mvn_generic <- function(x,
                          mu,
                          sigmacap) {
  stopifnot(
    is.vector(x),
    is.vector(mu),
    is.matrix(sigmacap)
  )
  d <- x - mu
  k <- length(as.vector(d))
  constant <- (k / 2) * log(2 * pi)
  term_01 <- log(det(sigmacap))
  if (k == 1) {
    qcap <- 1 / as.vector(sigmacap)
    term_02 <- d^2 * qcap
  } else {
    qcap <- solve(sigmacap)
    term_02 <- crossprod(d, qcap) %*% d
  }
  return(
    as.vector(
      constant - 0.5 * (
        term_01 + term_02
      )
    )
  )
}
