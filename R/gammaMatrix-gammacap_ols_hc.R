#' Asymptotic Covariance Matrix for Ordinary Least Squares Regression
#' with Sandwich Type Adjustments
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param x Object.
#'   Object of class `lm`.
#' @inheritParams gammacap_ols_hc_qcap
#' @inheritParams gammacap_mvn
#'
#' @references
#'   Dudgeon, P. (2017).
#'   Some improvements in confidence intervals for standardized regression coefficients.
#'   Psychometrika,
#'   82, 928-951.
#'   doi:[10.1007/s11336-017-9563-z](https://doi.org/10.1007/s11336-017-9563-z).
#'
#' @returns A matrix
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
#' gammacap_ols_hc(obj)
#' @export
#' @family Gamma Matrix Functions
#' @keywords gammaMatrix
gammacap_ols_hc <- function(x,
                            type = "hc5",
                            g1 = 1,
                            g2 = 1.5,
                            constant = 0.7,
                            names = TRUE,
                            sep = ".") {
  stopifnot(methods::is(x, "lm"))
  h <- stats::hatvalues(x)
  data <- x$model
  return(
    gammacap_ols_hc_generic(
      x = data,
      h = h,
      type = type,
      g1 = g1,
      g2 = g2,
      constant = constant,
      names = names,
      sep = sep
    )
  )
}
