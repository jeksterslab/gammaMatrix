#' Helper Function to Convert a Vector of Parameters
#' to the Mean Vector and Covariance Matrix
#'
#' @param x Numeric vector.
#'   Vector of parameters
#'   \eqn{
#'     \boldsymbol{\theta}
#'     =
#'     \{
#'       \boldsymbol{\mu} ,
#'       \mathrm{vech} \left( \boldsymbol{\Sigma} \right)
#'     \}^{\prime}
#'   }.
#'
#' @returns A list with elements `mu` and `sigmacap`.
#'
#' @examples
#' mvn_theta_helper(
#'   x = c(0, 0, 1, 0.5, 1)
#' )
#' @export
#' @family Multivariate Normal Distribution Functions
#' @keywords multiNorm
mvn_theta_helper <- function(x) {
  stopifnot(is.vector(x))
  q <- length(x)
  k <- 0.5 * (sqrt(8 * q + 9) - 3)
  # sigmacap is always a matrix even if dimension is 1 x 1
  sigmacap <- matrix(
    data = 0,
    nrow = k,
    ncol = k
  )
  if (nrow(sigmacap) != k) {
    stop("Length of \"x\" is not valid.")
  }
  sigmacap[lower.tri(sigmacap, diag = TRUE)] <- x[(k + 1):q]
  sigmacap[upper.tri(sigmacap)] <- t(sigmacap)[upper.tri(sigmacap)]
  mu <- x[1:k]
  list(
    mu = mu,
    sigmacap = sigmacap
  )
}
