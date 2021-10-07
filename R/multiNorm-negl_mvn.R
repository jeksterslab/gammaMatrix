#' Negative Log of the Likelihood of the Multivariate Normal Distribution
#'
#' Calculates the negative log of the likelihood function
#' of the multivariate normal distribution
#' for the ith observation.
#'
#' @author Ivan Jacob Agaloos Pesigan
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
#'   FUN = negl_mvn,
#'   theta = theta
#' )
#' @export
#' @family Multivariate Normal Distribution Functions
#' @keywords multiNorm
negl_mvn <- function(theta,
                     x) {
  return(
    -1 * l_mvn(
      theta = theta,
      x = x
    )
  )
}
