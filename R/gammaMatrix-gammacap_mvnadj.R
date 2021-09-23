#' Asymptotic Covariance Matrix with Adjustment
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
#' gammacap_mvnadj(x)
#' @export
#' @family Gamma Matrix Functions
#' @keywords gammaMatrix
gammacap_mvnadj <- function(x,
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
    varnames <- colnames(output)
  }
  output <- sym_of_vech(vech(output))
  colnames(output) <- rownames(output) <- varnames
  return(output)
}
