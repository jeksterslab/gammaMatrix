#' Leverage Adjustment - Generic
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param h Numeric vector.
#'   Leverage values.
#' @param k Positive integer.
#'   `p` number of regressors plus 1.
#' @inheritParams gammacap_ols_hc_qcap
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
#' head(gammacap_ols_hc_qcap_generic(h = h, k = 2))
#' @importFrom methods is
#' @export
#' @family Gamma Matrix Functions
#' @keywords gammaMatrix
gammacap_ols_hc_qcap_generic <- function(h, # leverage
                                         k, # p + 1
                                         type = "hc3",
                                         g1 = 1,
                                         g2 = 1.5,
                                         constant = 0.7) {
  n <- length(h)
  if (type %in% c("hc0", "hc1")) {
    return(
      rep(x = 1, times = n)
    )
  }
  if (type == "hc2") {
    return(
      1 / ((1 - h)^1)
    )
  }
  if (type == "hc3") {
    return(
      1 / ((1 - h)^2)
    )
  }
  if (type == "hc4") {
    delta <- sapply(
      X = h,
      FUN = function(i) {
        min(4, (n * i / k))
      }
    )
    return(
      1 / ((1 - h)^delta)
    )
  }
  if (type == "hc4m") {
    lambda <- sapply(
      X = h,
      FUN = function(i) {
        tmp <- n * i / k
        min(g1, tmp) + min(g2, tmp)
      }
    )
    return(
      1 / ((1 - h)^lambda)
    )
  }
  if (type == "hc5") {
    tmp <- n * constant * max(h) / k
    gamma <- sapply(
      X = h,
      FUN = function(i) {
        min((n * i / k), max(4, tmp))
      }
    )
    return(
      1 / sqrt((1 - h)^gamma)
    )
  }
}
