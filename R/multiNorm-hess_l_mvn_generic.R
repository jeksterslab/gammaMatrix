#' Hessian Matrix of the Multivariate Normal Distribution - Generic
#'
#' Calculates hessian matrix of the log of the likelihood function
#' of the multivariate normal distribution
#' for the ith observation.
#'
#' @param x Numeric vector of length `k`.
#'   The ith vector of observations.
#' @param mu Numeric vector.
#'   Parameter.
#'   Mean vector
#'   \eqn{\boldsymbol{\mu}}.
#' @param sigmacap Numeric matrix.
#'   Parameter.
#'   Covariance matrix
#'   \eqn{\boldsymbol{\Sigma}}.
#'
#' @returns A matrix.
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
#'   FUN = hess_l_mvn_generic,
#'   mu = mu,
#'   sigmacap = sigmacap
#' )
#' @export
#' @family Multivariate Normal Distribution Functions
#' @keywords multiNorm
hess_l_mvn_generic <- function(x,
                               mu,
                               sigmacap) {
  stopifnot(
    is.vector(x),
    is.vector(mu),
    is.matrix(sigmacap)
  )
  d <- x - mu
  k <- length(as.vector(d))
  q <- k + (k * (k + 1) / 2)
  if (k == 1) {
    sigmacap <- as.vector(sigmacap)
    wrt_mu_mu <- -1 / sigmacap
    wrt_mu_sigmacap <- -1 * (d / sigmacap^2)
    # wrt_sigmacap_sigmacap <- (
    #  1 / (2 * sigmacap^2)
    # ) + (
    #  d^2 / sigmacap^2
    # )
    wrt_sigmacap_sigmacap <- -1 * (
      (d^2 / sigmacap^3) - (0.5 / sigmacap^2)
    )
    hcap <- matrix(
      data = 0,
      nrow = q,
      ncol = q
    )
    hcap[1, 1] <- wrt_mu_mu
    hcap[1, 2] <- wrt_mu_sigmacap
    hcap[2, 1] <- wrt_mu_sigmacap
    hcap[2, 2] <- wrt_sigmacap_sigmacap
  } else {
    qcap <- chol2inv(chol(sigmacap))
    dcap <- dcap(dim(qcap)[1])
    tdcap <- t(dcap)
    wrt_mu_mu <- -1 * qcap
    wrt_mu_sigmacap <- kronecker(
      t(
        -1 * (
          qcap %*% d
        )
      ),
      qcap
    )
    ############################################
    # wrt_mu_sigmacap <- lapply(
    #  X = seq_len(k),
    #  FUN = function(i) {
    #    wrt_mu_sigmacap[i, ]
    #  }
    # )
    # foo <- function(x) {
    #  k <- sqrt(length(x))
    #  x <- matrix(
    #    x,
    #    nrow = k,
    #    ncol = k
    #  )
    #  diags <- diag(x)
    #  x <- x + t(x)
    #  diag(x) <- diags
    #  x <- as.matrix(vech(x))
    #  return(x)
    # }
    # wrt_mu_sigmacap <- lapply(
    #  X = wrt_mu_sigmacap,
    #  FUN = foo
    # )
    wrt_mu_sigmacap <- lapply(
      X = seq_len(k),
      FUN = function(i) {
        tdcap %*% wrt_mu_sigmacap[i, ]
      }
    )
    wrt_mu_sigmacap <- do.call(
      what = "cbind",
      args = wrt_mu_sigmacap
    )
    #######################################
    # page 353 matrix algebra for econometrics
    acap <- (
      qcap %*% (
        2 * tcrossprod(d) - sigmacap
      ) %*% qcap
    )
    wrt_sigmacap_sigmacap <- -1 * (
      0.5 * crossprod(
        dcap,
        kronecker(
          qcap,
          acap
        )
      ) %*% dcap
    )
    hcap <- matrix(
      data = 0,
      nrow = q,
      ncol = q
    )
    hcap[(1:k), (1:k)] <- wrt_mu_mu
    hcap[((k + 1):q), (1:k)] <- wrt_mu_sigmacap
    hcap[(1:k), ((k + 1):q)] <- t(wrt_mu_sigmacap)
    hcap[((k + 1):q), ((k + 1):q)] <- wrt_sigmacap_sigmacap
  }
  return(
    hcap
  )
}
