#' Asymptotic Covariance Matrix with Adjustment (Variant 1)
#'
#' Calculates the covariance matrix of the unique elements
#' of the covariance matrix with adjustment for nonnormality.
#'
#' @details
#' # Dependencies
#' * [rmvn_chol()] (test)
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @inheritParams gammacap_gen
#'
#' @references
#'   Add appropriate references here...
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
#' gammacap_mvnadj1(x)
#' @export
#' @family Gamma Matrix Functions
#' @keywords gammaMatrix
gammacap_mvnadj1 <- function(x,
                             names = TRUE,
                             sep = ".") {
  mvn <- gammacap_mvn(
    x = x,
    sigmacap = NULL,
    names = names,
    sep = sep
  )
  gen <- gammacap_gen(
    x = x,
    names = names,
    sep = sep
  )
  invmvn <- chol2inv(chol(mvn))
  output <- (
    mvn %*% (
      invmvn %*% gen %*% invmvn
    ) %*% mvn
  )
  if (names) {
    colnames(output) <- rownames(output) <- colnames(gen)
  }
  return(output)
}
