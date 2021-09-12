#' Asymptotic Distribution-Free Covariance Matrix
#'
#' Calculates the asymptotic distribution-free (ADF)
#' covariance matrix of the unique elements
#' of a sample covariance matrix.
#'
#' @details
#' # Dependencies
#' * [rmvn_chol()] (test)
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param unbiased Logical.
#'   If `unbiased = TRUE`,
#'   returns unbiased asymptotic distribution-free covariance matrix.
#'   If `unbiased = FALSE`,
#'   returns consistent asymptotic distribution-free covariance matrix.
#' @inheritParams gammacap_gen
#'
#' @references
#'   Browne, M. W. (1984).
#'   Asymptotically distribution-free methods
#'   for the analysis of covariance structures.
#'   British Journal of Mathematical and Statistical Psychology,
#'   37, 62–83.
#'   doi:[10.1111/j.2044-8317.1984.tb00789.x](https://doi.org/10.1111/j.2044-8317.1984.tb00789.x).
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
#' gammacap_adf(x, unbiased = TRUE)
#' gammacap_adf(x, unbiased = FALSE)
#' @export
#' @family Gamma Matrix Functions
#' @keywords gammaMatrix
gammacap_adf <- function(x,
                         unbiased = TRUE,
                         names = TRUE,
                         sep = ".") {
  stopifnot(
    is.data.frame(x) || is.matrix(x) || is.vector(x),
    is.logical(unbiased)
  )
  if (is.vector(x)) {
    x <- as.matrix(x)
  }
  n <- dim(x)[1]
  # ML covariance
  sigmacaptilde <- stats::cov(x) * (n - 1) / n
  d <- scale(
    x = x,
    center = TRUE,
    scale = FALSE
  )
  # Browne (1984) page 71 equation 3.4
  gammacaptilde <- (1 / n) * Reduce(
    f = "+",
    x = lapply(
      X = seq_len(n),
      FUN = function(i) {
        tcrossprod(
          vech(
            tcrossprod(d[i, ])
          )
        )
      }
    )
  ) - tcrossprod(
    vech(sigmacaptilde)
  )
  if (unbiased) {
    output <- gammacap_adfunbiased(
      gammacaptilde = gammacaptilde,
      sigmacaptilde = sigmacaptilde,
      n = n
    )
  } else {
    output <- gammacaptilde
  }
  if (names) {
    colnames(output) <- rownames(output) <- gammacapnames(x, sep = sep)
  }
  return(output)
}

#' @rdname gammacap_adf
#' @param gammacaptilde Numeric matrix.
#'   Consistent estimate of the asymptotic distribution-free
#'   covariance matrix.
#' @param sigmacaptilde Numeric matrix.
#'   Consistent estimate of the sample covariance matrix.
#' @param n Positive integer. Sample size.
#' @export
gammacap_adfunbiased <- function(gammacaptilde,
                                 sigmacaptilde,
                                 n) {
  # Browne (1984) page 72 equation 3.8
  return(
    (
      (
        (n * (n - 1)) / ((n - 2) * (n - 3))
      ) * gammacaptilde
    ) - (
      (
        n / ((n - 2) * (n - 3))
      ) * (
        gammacap_mvn(sigmacap = sigmacaptilde) - (
          (2 / (n - 1)) * tcrossprod(
            vech(sigmacaptilde)
          )
        )
      )
    )
  )
}
