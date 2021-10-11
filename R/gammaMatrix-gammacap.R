#' Asymptotic Covariance Matrix
#'
#' Calculates the
#' covariance matrix of the unique elements
#' of the covariance matrix.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param x Numeric matrix, data frame, or vector.
#' @param type Character string.
#'   If `type = "adf"`,
#'   calculate asymptotic distribution-free covariance matrix.
#'   If `type = "adfnb"`,
#'   calculate nonparametric bootstrapped
#'   asymptotic distribution-free covariance matrix.
#'   If `type = "gen"`,
#'   calculate covariance matrix using the general formula.
#'   If `type = "mvn"`,
#'   calculate covariance matrix with multivariate normal data.
#'   If `type = "mvnadj1"`,
#'   calculate covariance matrix with adjustment variant 1.
#'   If `type = "mvnadj2"`,
#'   calculate covariance matrix with adjustment variant 2.
#'   If `type = "nb"`,
#'   calculate covariance matrix from
#'   nonparametric bootstrapped
#'   covariances.
#' @param sigmacap Numeric matrix.
#'   The argument is used when `type = "mvn"`.
#'   Optional argument.
#'   Sample covariance matrix.
#' @param unbiased Logical.
#'   The argument is used when `type = "adf"`.
#'   If `unbiased = TRUE`,
#'   returns unbiased asymptotic distribution-free covariance matrix.
#'   If `unbiased = FALSE`,
#'   returns consistent asymptotic distribution-free covariance matrix.
#' @param bcap Integer.
#'   The argument is used when `type = "adfnb" or type = "nb"`.
#'   Number of bootstrap samples.
#' @param missing Logical.
#'   The argument is used when `type = "mvnadj2"`.
#'   If `missing = TRUE`,
#'   the mean vector and the covariance matrix
#'   will be estimated using the EM algorithm.
#'   If `missing = FALSE`,
#'   all missing values will be dropped
#'   and the mean vector and covariance matrix
#'   will be estimated using `x`.
#' @param ml_cov Logical.
#'   The argument is used when `type = "mvnadj2"`.
#'   If `missing = FALSE` and `ml_cov = TRUE`,
#'   use maximum likelihood estimator of the covariance matrix.
#' @param seed Integer.
#'   The argument is used when `type = "adfnb" or type = "nb"`.
#'   Random number generation seed.
#' @inheritParams vech
#'
#' @returns A matrix.
#'
#' @examples
#' set.seed(42)
#' n <- 1000
#' k <- 2
#' q <- chol(
#'   matrix(
#'     data = c(1.0, 0.5, 0.5, 1.0),
#'     nrow = k, ncol = k
#'   )
#' )
#' z <- matrix(
#'   data = rnorm(n = n * k), nrow = n, ncol = k
#' )
#' x <- z %*% q
#'
#' gammacap(x, type = "adf")
#' gammacap(x, type = "adfnb")
#' gammacap(x, type = "gen")
#' gammacap(x, type = "mvn")
#' gammacap(x, type = "mvnadj1")
#' gammacap(x, type = "mvnadj2")
#' gammacap(x, type = "nb")
#' gammacap(sigmacap = cov(x), type = "mvn")
#' @export
#' @family Gamma Matrix Functions
#' @keywords gammaMatrix
gammacap <- function(x,
                     type = "mvn",
                     sigmacap = NULL,
                     unbiased = TRUE,
                     bcap = 1000L,
                     seed = NULL,
                     missing = FALSE,
                     ml_cov = TRUE,
                     names = TRUE,
                     sep = ".") {
  # gammaR
  switch(
    EXPR = type,
    gen = {
      return(
        gammacap_gen(
          x = x,
          names = names,
          sep = sep
        )
      )
    },
    mvn = {
      return(
        gammacap_mvn(
          x = x,
          sigmacap = sigmacap,
          names = names,
          sep = sep
        )
      )
    },
    mvnadj1 = {
      return(
        gammacap_mvnadj1(
          x = x,
          names = names,
          sep = sep
        )
      )
    },
    mvnadj2 = {
      return(
        gammacap_mvnadj2(
          x = x,
          missing = FALSE,
          ml_cov = ml_cov,
          drop_means = TRUE, # only elements for the covariances for uniformity
          names = names,
          sep = sep
        )
      )
    },
    adf = {
      return(
        gammacap_adf(
          x = x,
          unbiased = unbiased,
          names = names,
          sep = sep
        )
      )
    },
    adfnb = {
      return(
        gammacap_adfnb(
          x = x,
          bcap = bcap,
          seed = seed,
          names = names,
          sep = sep
        )
      )
    },
    nb = {
      return(
        gammacap_nb(
          x = x,
          bcap = bcap,
          seed = seed,
          names = names,
          sep = sep
        )
      )
    }
  )
}
