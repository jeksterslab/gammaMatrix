#' Nonparametric Bootstrapped Covariance Matrix
#'
#' Calculates the nonparametric bootstrapped
#' asymptotic distribution-free (ADF)
#' covariance matrix of the unique elements
#' of the covariance matrix.
#'
#' @details
#' # Dependencies
#' * [rmvn_chol()] (test)
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param bcap Integer.
#'   Number of bootstrap samples.
#' @param seed Integer.
#'   Random number generation seed.
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
#'   Yung, Y.-F., & Bentler, P. M. (1994).
#'   Bootstrap-corrected ADF test statistics
#'   in covariance structure analysis.
#'   British Journal of Mathematical and Statistical Psychology,
#'   47, 63–84.
#'   doi:[10.1111/j.2044-8317.1994.tb01025.x](https://doi.org/10.1111/j.2044-8317.1994.tb01025.x).
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
#' gammacap_adfnb(x)
#' @export
#' @family Gamma Matrix Functions
#' @keywords gammaMatrix
gammacap_nb <- function(x,
                        bcap = 5000L,
                        seed = NULL,
                        names = TRUE,
                        sep = ".") {
  stopifnot(
    is.data.frame(x) || is.matrix(x) || is.vector(x)
  )
  bcap <- as.integer(bcap)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.vector(x)) {
    x <- as.matrix(x)
  }
  n <- dim(x)[1]
  k <- dim(x)[2]
  if (k == 1) {
    nb <- function(i, x) {
      stats::var(
        x[sample(
          seq_len(dim(x)[1]),
          replace = TRUE
        ), ]
      )
    }
  } else {
    nb <- function(i, x) {
      vech(
        stats::var(
          x[sample(
            seq_len(dim(x)[1]),
            replace = TRUE
          ), ]
        )
      )
    }
  }
  vechsigmacap <- lapply(
    X = seq_len(bcap),
    FUN = nb,
    x = x
  )
  vechsigmacap <- do.call(
    what = "rbind",
    args = vechsigmacap
  )
  output <- stats::var(vechsigmacap) * n
  if (names) {
    colnames(output) <- rownames(output) <- gammacapnames(x, sep = sep)
  }
  return(output)
}
