#' Asymptotic Covariance Matrix for Ordinary Least Squares Regression
#' with Sandwich Type Adjustments - Generic
#'
#' @details
#' # Dependencies
#' * [gammacap_ols_hc_qcap_generic()]
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param x Numeric matrix or data frame.
#'   Data matrix for \eqn{ \{ y, x_1, \cdots, x_p \}}.
#' @inheritParams gammacap_ols
#' @inheritParams gammacap_ols_hc_qcap_generic
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
#' h <- hatvalues(obj)
#'
#' gammacap_ols_hc_generic(x, h = h)
#' @export
#' @family Gamma Matrix Functions
#' @keywords gammaMatrix
gammacap_ols_hc_generic <- function(x,
                                    h = NULL,
                                    type = "hc5",
                                    g1 = 1,
                                    g2 = 1.5,
                                    constant = 0.7,
                                    names = TRUE,
                                    sep = ".") {
  stopifnot(
    is.data.frame(x) || is.matrix(x)
  )
  x <- as.matrix(x)
  d <- scale(
    x = x,
    center = TRUE,
    scale = FALSE
  )
  sigmacap <- stats::cov(x)
  n <- dim(x)[1]
  k <- dim(x)[2]
  if (is.null(h)) {
    design_matrix <- cbind(1, x[, -1, drop = FALSE])
    # leverage
    h <- diag(
      design_matrix %*% solve(
        crossprod(design_matrix)
      ) %*% t(design_matrix)
    )
  }
  qcap <- gammacap_ols_hc_qcap_generic(
    h = h,
    k = k,
    type = type,
    g1 = g1,
    g2 = g2,
    constant = constant
  )
  output <- (
    (
      1 / n
    ) * Reduce(
      f = "+",
      x = lapply(
        X = seq_len(n),
        FUN = function(i) {
          qcap[i] * tcrossprod(
            vech(
              tcrossprod(
                d[i, ]
              ) - sigmacap
            )
          )
        }
      )
    )
  )
  if (names) {
    colnames(output) <- rownames(output) <- gammacapnames(x, sep = sep)
  }
  return(output)
}
