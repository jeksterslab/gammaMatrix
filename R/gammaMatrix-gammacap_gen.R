#' Asymptotic Covariance Matrix (General)
#'
#' Calculates the asymptotic covariance matrix of the unique elements
#' of the covariance matrix.
#'
#' @details
#' # Dependencies
#' * [vech()]
#' * [vechnames()]
#' * [rmvn_chol()] (test)
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param x Numeric matrix, data frame, or vector.
#' @inheritParams vech
#'
#' @references
#'   Yuan, K.-H., & Hayashi, K. (2006).
#'   Standard errors in covariance structure models:
#'   Asymptotics versus bootstrap.
#'   British Journal of Mathematical and Statistical Psychology,
#'   59, 397–417.
#'   doi:[10.1348/000711005X85896](https://doi.org/10.1348/000711005X85896).
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
#' x <- z %*% q
#'
#' gammacap_gen(x)
#' @export
#' @family Gamma Matrix Functions
#' @keywords gammaMatrix
gammacap_gen <- function(x,
                         names = TRUE,
                         sep = ".") {
  stopifnot(
    is.data.frame(x) || is.matrix(x) || is.vector(x)
  )
  if (is.vector(x)) {
    x <- as.matrix(x)
  }
  d <- scale(
    x = x,
    center = TRUE,
    scale = FALSE
  )
  sigmacap <- stats::cov(x)
  n <- dim(x)[1]
  # Yuan and Hayashi (2006) page 496
  # Yuan and Chan (2011) page 674 equation 12
  output <- (
    1 / n
  ) * Reduce(
    f = "+",
    x = lapply(
      X = seq_len(n),
      FUN = function(i) {
        tcrossprod(
          vech(
            tcrossprod(
              d[i, ]
            ) - sigmacap
          )
        )
      }
    )
  )
  if (names) {
    colnames(output) <- rownames(output) <- gammacapnames(x, sep = sep)
  }
  return(output)
}
