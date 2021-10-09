#' Gradient Vector of the Multivariate Normal Distribution - Generic
#'
#' Calculates gradient vector of the log of the likelihood function
#' of the multivariate normal distribution
#' for the ith observation.
#'
#' @inheritParams l_mvn
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
#'   FUN = grad_l_mvn,
#'   theta = theta
#' )
#' @export
#' @family Multivariate Normal Distribution Functions
#' @keywords multiNorm
grad_l_mvn <- function(x,
                       theta) {
  theta <- mvn_theta_helper(theta)
  return(
    grad_l_mvn_generic(
      x = x,
      mu = theta$mu,
      sigmacap = theta$sigmacap
    )
  )
}
