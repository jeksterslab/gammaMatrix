#' Asymptotic Covariance Matrix Assuming Multivariate Normal Data
#'
#' Calculates the covariance matrix of the unique elements
#' of the covariance matrix assuming
#' multivariate normal data.
#'
#' @details
#' # Dependencies
#' * [dcap()]
#' * [kcap()]
#' * [rmvn_chol()] (test)
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param x Numeric matrix, data frame, or vector.
#' @param sigmacap Numeric matrix.
#'   Optional argument.
#'   Sample covariance matrix.
#' @inheritParams vech
#'
#' @references
#'   Browne, M. W., & Arminger, G. (1995).
#'   Specification and estimation of mean-and covariance-structure  models.
#'   Handbook of statistical modeling forthe social and behavioral sciences.
#    Springer US.
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
#' gammacap_mvn(x)
#' gammacap_mvn(sigmacap = cov(x))
#' @export
#' @family Gamma Matrix Functions
#' @keywords gammaMatrix
gammacap_mvn <- function(x,
                         sigmacap = NULL,
                         names = TRUE,
                         sep = ".") {
  if (is.null(sigmacap)) {
    stopifnot(
      is.data.frame(x) || is.matrix(x) || is.vector(x)
    )
    if (is.vector(x)) {
      sigmacap <- as.matrix(stats::var(x))
    } else {
      sigmacap <- stats::cov(x)
    }
  } else {
    stopifnot(is.matrix(sigmacap))
  }
  stopifnot(
    dim(sigmacap)[1] == dim(sigmacap)[2],
    sigmacap == t(sigmacap)
  )
  kcap <- kcap(dim(sigmacap)[1])
  # Browne and Arminger (1995) equation 4.9 on page 189
  output <- 2 * kcap %*% tcrossprod(
    kronecker(
      sigmacap,
      sigmacap
    ),
    kcap
  )
  if (names) {
    colnames(output) <- rownames(output) <- gammacapnames(sigmacap, sep = sep)
  }
  return(output)
}
