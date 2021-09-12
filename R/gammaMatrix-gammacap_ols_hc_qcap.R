#' Leverage Adjustment
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param x Object.
#'   Object of class `lm`.
#' @param type Character string.
#'   Correction type.
#'   Possible values are
#'   `"hc0"`,
#'   `"hc1"`,
#'   `"hc2"`,
#'   `"hc3"`,
#'   `"hc4"`,
#'   `"hc4m"`, and
#'   `"hc5"`.
#' @param g1 Numeric.
#'   `g1` value for `type = "hc4m"` or `type = "hc5"`.
#' @param g2 Numeric.
#'   `g2` value for `type = "hc4m"`.
#' @param constant Numeric.
#'   Constant for `type = "hc5"`
#'
#' @references
#'   Dudgeon, P. (2017).
#'   Some improvements in confidence intervals for standardized regression coefficients.
#'   Psychometrika,
#'   82, 928-951.
#'   doi:[10.1007/s11336-017-9563-z](https://doi.org/10.1007/s11336-017-9563-z).
#'
#' @returns A vector.
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
#' h <- hatvalues(obj)
#' head(1 / ((1 - h)^2))
#'
#' head(gammacap_ols_hc_qcap(obj, type = "hc3"))
#' @importFrom methods is
#' @export
#' @family Gamma Matrix Functions
#' @keywords gammaMatrix
gammacap_ols_hc_qcap <- function(x,
                                 type = "hc3",
                                 g1 = 1,
                                 g2 = 1.5,
                                 constant = 0.7) {
  stopifnot(methods::is(x, "lm"))
  k <- ncol(x$model)
  h <- stats::hatvalues(x)
  return(
    gammacap_ols_hc_qcap_generic(
      h = h,
      k = k,
      type = type,
      g1 = g1,
      g2 = g2,
      constant = constant
    )
  )
}
