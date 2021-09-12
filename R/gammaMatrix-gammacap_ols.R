#' Asymptotic Covariance Matrix for Ordinary Least Squares Regression
#'
#' Calculates the asymptotic covariance matrix of the unique elements
#' of a sample covariance matrix using ordinary least squares estimates.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param x Object of class `lm`.
#' @param type Character string.
#'   Type of asymptotic covariance matrix of the regressors
#'   If `type = "adf"`,
#'   calculate asymptotic distribution-free covariance matrix.
#'   If `type = "adfnb"`,
#'   calculate nonparametric bootstrapped
#'   asymptotic distribution-free covariance matrix.
#'   If `type = "gen"`,
#'   calculate covariance matrix using the general formula.
#'   If `type = "mvn"`,
#'   calculate covariance matrix with multivariate normal data.
#' @param adf_unbiased Logical.
#'   If `adf_unbiased = TRUE`,
#'   use unbiased asymptotic distribution-free covariance matrix.
#'   If `adf_unbiased = FALSE`,
#'   use consistent asymptotic distribution-free covariance matrix.
#' @param bcap Integer.
#'   Number of bootstrap samples
#'   when `type = "gen"` and `type = "adfnb"`.
#' @param seed Integer.
#'   Random number generation seed
#'   when `type = "gen"` and `type = "adfnb"`.
#' @param ke_unbiased Logical.
#'   If `ke_unbiased = TRUE`,
#'   use estimator
#'   \eqn{
#'     \frac{
#'       \frac{1}{n}
#'       \sum_{i=1}^{n}
#'       \hat{\varepsilon}_{i}^{4}
#'     }{
#'       \left(
#'         \frac{1}{n - p - 1}
#'         \sum_{i=1}^{n}
#'         \hat{\varepsilon}_{i}^{2}
#'       \right)^2
#'     }
#'   }
#'   where the denominator is the square of the unbiased estimator
#'   of the error variance.
#'   If `ke_unbiased = FALSE`,
#'   use estimator
#'   \eqn{
#'     \frac{
#'       \frac{1}{n}
#'       \sum_{i=1}^{n}
#'       \hat{\varepsilon}_{i}^{4}
#'     }{
#'       \left(
#'         \frac{1}{n}
#'         \sum_{i=1}^{n}
#'         \hat{\varepsilon}_{i}^{2}
#'       \right)^2
#'     }
#'   }
#'   where the denominator is the square of the consistent estimator
#'   of the error variance.
#' @param yc Logical.
#'   If `yc = TRUE`, order the output
#'   following Yuan and Chan (2011) page 674.
#'   If `yc = FALSE`, the order of the output
#'   following \eqn{\mathrm{vech} \left( \Sigma \right)}.
#' @inheritParams vech
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
#'
#' gammacap_ols(obj)
#' @importFrom methods is
#' @export
#' @family Gamma Matrix Functions
#' @keywords gammaMatrix
gammacap_ols <- function(x,
                         type = "gen",
                         adf_unbiased = TRUE,
                         bcap = 5000L,
                         seed = NULL,
                         ke_unbiased = FALSE,
                         yc = FALSE,
                         names = TRUE,
                         sep = ".") {
  stopifnot(methods::is(x, "lm"))
  sigmasq <- stats::summary.lm(x)$sigma^2
  dfx <- as.data.frame(x$model[, -1, drop = FALSE])
  beta <- stats::coef(x)[-1]
  varnames <- names(beta)
  beta <- as.vector(beta)
  names(beta) <- varnames
  sigmacap <- as.matrix(stats::cov(x$model))
  # kurtosis of errors
  if (ke_unbiased) {
    ke <- mean(x$residuals^4) / sigmasq^2
  } else {
    # sigmasq should be unbiased estimate of error variance
    # that is, sum of squares errors over n - p - 1
    ke <- mean(x$residuals^4) / mean(x$residuals^2)^2
  }
  return(
    gammacap_ols_generic(
      x = dfx,
      beta = beta,
      sigmacap = sigmacap,
      sigmasq = sigmasq,
      ke = ke,
      type = type,
      adf_unbiased = adf_unbiased,
      bcap = bcap,
      seed = seed,
      yc = yc,
      names = names,
      sep = sep
    )
  )
}
