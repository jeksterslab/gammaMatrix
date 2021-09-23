#' Asymptotic Covariance Matrix for Ordinary Least Squares Regression - Generic
#'
#' Calculates the asymptotic covariance matrix of the unique elements
#' of the covariance matrix using ordinary least squares estimates.
#'
#' @details
#' # Dependencies
#' * [dcap()]
#' * [vechnames()]
#' * [rmvn_chol()] (test)
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param x Numeric matrix, data frame, or vector.
#'   Values for regressor variables.
#' @param beta Numeric vector.
#'   Partial regression slopes.
#' @param sigmacap Numeric matrix.
#'   Covariance matrix
#'   \eqn{\boldsymbol{\Sigma}}
#'   of
#'   \eqn{\{y, x_1, \cdots, x_p \}^{\prime}}.
#' @param sigmasq Numeric.
#'   Error variance
#'   \eqn{\sigma^2}.
#' @param ke Numeric.
#'   Kurtosis of errors.
#' @inheritParams gammacap_ols
#'
#' @references
#'   Yuan, K.-H., & Chan, W. (2011).
#'   Biases and standard errors of standardized regression coefficients.
#'   Psychometrika,
#'   76, 670–690.
#'   doi:[10.1007/S11336-011-9224-6](https://doi.org/10.1007/S11336-011-9224-6).
#'
#' @returns A matrix.
#'
#' @examples
#' set.seed(42)
#' n <- 1000
#' k <- 2
#' z <- matrix(
#'   data = rnorm(n = n * k), nrow = n, ncol = k
#' )
#' q <- chol(
#'   matrix(
#'     data = c(1.0, 0.5, 0.5, 1.0),
#'     nrow = k, ncol = k
#'   )
#' )
#' x <- as.data.frame(z %*% q)
#' colnames(x) <- c("y", "x")
#' obj <- lm(y ~ x, data = x)
#' beta <- coef(obj)[-1]
#' sigmacap <- cov(x)
#' sigmasq <- summary(obj)$sigma^2
#' ke <- 3
#'
#' gammacap_ols_generic(
#'   x[, -1],
#'   beta = beta,
#'   sigmacap = sigmacap,
#'   sigmasq = sigmasq,
#'   ke = ke
#' )
#' @importFrom methods is
#' @export
#' @family Gamma Matrix Functions
#' @keywords gammaMatrix
gammacap_ols_generic <- function(x,
                                 beta,
                                 sigmacap,
                                 sigmasq,
                                 ke,
                                 type = "gen",
                                 adf_unbiased = TRUE,
                                 bcap = 5000L,
                                 seed = NULL,
                                 yc = FALSE,
                                 names = TRUE,
                                 sep = ".") {
  stopifnot(
    is.vector(beta),
    is.vector(sigmasq),
    length(sigmasq) == 1,
    is.matrix(sigmacap),
    dim(sigmacap)[1] == dim(sigmacap)[2],
    dim(sigmacap)[1] - 1 == length(beta),
    is.vector(ke),
    length(ke) == 1
  )
  # estimates
  p <- length(beta)
  k <- p + 1
  tbeta <- t(beta)
  sigmacapx <- sigmacap[2:k, 2:k, drop = FALSE]
  icap <- diag(p)
  dcap <- dcap(p)
  tdcap <- t(dcap)
  # gammacapx is equal to gamma11
  gammacapx <- gammacap(
    x = x,
    type = type,
    unbiased = adf_unbiased,
    bcap = bcap,
    seed = seed
  )
  # estimate gammacap
  # the calculations below are based on Yuan and Chan page 674
  gammacapx_times_tdcap <- gammacapx %*% tdcap
  beta_kron_icap <- kronecker(beta, icap)
  tbeta_kron_icap <- kronecker(tbeta, icap)
  beta_kron_beta <- kronecker(beta, beta)
  tbeta_kron_tbeta <- kronecker(tbeta, tbeta)
  gamma12 <- gammacapx_times_tdcap %*% beta_kron_icap
  gamma13 <- gammacapx_times_tdcap %*% beta_kron_beta
  gamma22 <- (
    tbeta_kron_icap %*% dcap %*% gamma12
  ) + (
    sigmasq * sigmacapx
  )
  gamma23 <- (
    tbeta_kron_icap %*% dcap %*% gamma13
  ) + (
    2 * (sigmacapx %*% beta) * sigmasq
  )
  gamma33 <- (
    tbeta_kron_tbeta %*% dcap %*% gamma13
  ) + (
    4 * (tbeta %*% sigmacapx %*% beta) * sigmasq
  ) + (
    (ke - 1) * sigmasq^2
  )
  # names
  if (names) {
    varnames <- colnames(sigmacap)
    if (is.null(varnames)) {
      varnames <- paste0("v", seq_len(dim(sigmacap)[2]))
    }
    varnamesy <- varnames[1]
    varnamesx <- varnames[-1]
    varnamesxx <- vechnames(
      varnamesx,
      sep = sep
    )
    varnamesyx <- paste0(
      varnamesy,
      sep,
      varnamesx
    )
    varnamesyy <- paste0(
      varnamesy,
      sep,
      varnamesy
    )
    colnames(gammacapx) <- rownames(gammacapx) <- varnamesxx
    colnames(gamma12) <- varnamesyx
    rownames(gamma12) <- varnamesxx
    colnames(gamma13) <- varnamesyy
    rownames(gamma13) <- varnamesxx
    colnames(gamma22) <- rownames(gamma22) <- varnamesyx
    colnames(gamma23) <- varnamesyy
    rownames(gamma23) <- varnamesyx
    colnames(gamma33) <- rownames(gamma33) <- varnamesyy
  }
  if (yc) {
    # Follow the order from Yuan and Chan page 674
    output <- cbind(
      rbind(
        gammacapx,
        t(gamma12),
        t(gamma13)
      ),
      rbind(
        gamma12,
        gamma22,
        t(gamma23)
      ),
      rbind(
        gamma13,
        gamma23,
        gamma33
      )
    )
  } else {
    # Follow the order of vech of sigmacap
    output <- cbind(
      rbind(
        gamma33,
        gamma23,
        gamma13
      ),
      rbind(
        t(gamma23),
        gamma22,
        gamma12
      ),
      rbind(
        t(gamma13),
        t(gamma12),
        gammacapx
      )
    )
  }
  return(output)
}
