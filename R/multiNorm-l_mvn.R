#' Log of the Likelihood of the Multivariate Normal Distribution
#'
#' Calculates the log of the likelihood function
#' of the multivariate normal distribution
#' for the ith observation.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param theta Numeric vector.
#'   Parameter vector
#'   \eqn{
#'     \boldsymbol{\theta}
#'     =
#'     \left\{
#'     \boldsymbol{\mu},
#'     \mathrm{vech}
#'     \left(
#'     \boldsymbol{\Sigma}
#'     \right)
#'     \right\}^{\prime}
#' }
#' @param x Numeric vector of length `k`.
#'   The ith vector of observations.
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
#' theta <- c(
#'   mu,
#'   vech(sigmacap)
#' )
#'
#' lapply(
#'   X = xcap,
#'   FUN = l_mvn,
#'   theta = theta
#' )
#' @export
#' @family Multivariate Normal Distribution Functions
#' @keywords multiNorm
l_mvn <- function(theta,
                  x) {
  theta <- theta_helper(theta)
  return(
    l_mvn_generic(
      x = x,
      mu = theta$mu,
      sigmacap = theta$sigmacap
    )
  )
}

# l_mvn2 <- function(theta,
#                   x) {
#  k <- length(x)
#  mu <- theta[1:k]
#  sigmacap <- matrix(
#    theta[(k+1):length(theta)],
#    nrow = k,
#    ncol = k
#  )
#  return(
#    l_mvn_generic(
#      x = x,
#      mu = mu,
#      sigmacap = sigmacap
#    )
#  )
# }
