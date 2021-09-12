#' Generate Data
#' from the Multivariate Normal Distribution
#' Using the Cholesky Decomposition
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param n Positive integer.
#'   `n` variates.
#' @param mu Numeric vector.
#'   Parameter.
#'   Mean vector
#'   \eqn{\boldsymbol{\mu}}.
#' @param sigmacap Numeric matrix.
#'   Parameter.
#'   Covariance matrix
#'   \eqn{\boldsymbol{\Sigma}}.
#' @param varnames Character vector
#'   Optional variable names.
#' @param data_frame Logical.
#'   If `data_frame = TRUE`,
#'   returns a `data.frame`.
#'   If `data_frame = FALSE`,
#'   returns a `matrix`.
#'
#' @returns A matrix (`data_frame = FALSE`) or data.frame (`data_frame = TRUE`).
#'
#' @examples
#' mu <- c(0, 0)
#' sigmacap <- matrix(
#'   data = c(
#'     1, 0.5, 0.5, 1
#'   ),
#'   nrow = 2
#' )
#'
#' rmvn_chol(
#'   n = 5,
#'   mu = mu,
#'   sigmacap = sigmacap
#' )
#' @export
#' @family Multivariate Normal Distribution Functions
#' @keywords multiNorm
rmvn_chol <- function(n,
                      mu,
                      sigmacap,
                      varnames = NULL,
                      data_frame = FALSE) {
  stopifnot(
    is.matrix(sigmacap),
    is.vector(mu),
    length(n) == 1,
    is.logical(data_frame)
  )
  n <- as.integer(n)
  k <- dim(sigmacap)[1]
  stopifnot(
    k == dim(sigmacap)[2],
    k == length(mu),
    sigmacap == t(sigmacap)
  )
  qcap <- chol(sigmacap)
  zcap <- matrix(
    data = stats::rnorm(
      n = n * k
    ),
    nrow = n,
    ncol = k
  )
  ones <- matrix(
    data = 1,
    nrow = n,
    ncol = 1
  )
  output <- (
    zcap %*% qcap + (
      ones %*% mu
    )
  )
  if (!is.null(varnames)) {
    stopifnot(length(varnames) == k)
    colnames(output) <- varnames
  }
  if (data_frame) {
    output <- as.data.frame(output)
  }
  return(output)
}
