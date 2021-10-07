#' Asymptotic Covariance Matrix with Adjustment (Variant 2)
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
#' @param missing Logical.
#'   If `missing = TRUE`,
#'   the mean vector and the covariance matrix
#'   will be estimated using the EM algorithm.
#'   If `missing = FALSE`,
#'   all missing values will be dropped
#'   and the mean vector and covariance matrix
#'   will be estimated using `x`.
#' @param ml_cov Logical.
#'   If `missing = FALSE` and `ml_cov = TRUE`,
#'   use maximum likelihood estimator of the covariance matrix.
#' @param drop_means Logical.
#'   If `drop_means = TRUE`,
#'   remove the means from the output matrix.
#'   If `drop_means = FALSE`,
#'   the first `ncol(x)` elements of the output matrix
#'   will be names for the `ncol(x)` means.
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
#' gammacap_mvnadj2(x)
#' gammacap_mvnadj2(x, drop_means = FALSE)
#' @export
#' @family Gamma Matrix Functions
#' @keywords gammaMatrix
gammacap_mvnadj2 <- function(x,
                             missing = FALSE,
                             ml_cov = TRUE,
                             drop_means = TRUE,
                             names = TRUE,
                             sep = ".") {
  stopifnot(
    is.data.frame(x) || is.matrix(x),
    missing == FALSE # missing data under construction
  )
  # only supports >= 2 variables
  stopifnot(ncol(x) > 1)
  # check missing values
  n_orig <- nrow(x)
  #  if (n_orig - nrow(x[stats::complete.cases(x), , drop = FALSE]) == 0) {
  #      if (missing) {
  #        message(
  #          paste0(
  #            "There are no missing observations in `x`.\n",
  #            "The function will perform the complete data approach."
  #          )
  #        )
  #      }
  #      missing <- FALSE
  #  } else {
  #      if (!missing) {
  #        message(
  #          paste0(
  #            "The complete data approach will be performed.\n",
  #            "List-wise deletion will be applied to remove rows with missing observations."
  #          )
  #        )
  #      }
  #  }
  if (missing) {
  } else {
    x <- x[stats::complete.cases(x), , drop = FALSE]
    n <- nrow(x)
    k <- ncol(x)
    mu <- colMeans(x)
    sigmacap <- stats::cov(x)
    if (ml_cov) {
      sigmacap <- sigmacap * (n - 1) / n
    }
    gradient_vector <- function(data,
                                theta,
                                q,
                                k) {
      mu <- theta[1:k]
      sigmacap <- sym_of_vech(theta[(k + 1):q])
      return(
        tcrossprod(
          grad_l_mvn_generic(
            x = data,
            mu = mu,
            sigmacap = sigmacap
          )
        )
      )
    }
    hessian_matrix <- function(data,
                               theta,
                               q,
                               k) {
      mu <- theta[1:k]
      sigmacap <- sym_of_vech(theta[(k + 1):q])
      return(
        hess_l_mvn_generic(
          x = data,
          mu = mu,
          sigmacap = sigmacap
        )
      )
    }
    theta <- c(
      mu,
      vech(sigmacap)
    )
    q <- length(theta)
    tx <- as.data.frame(t(x))
    gradient <- lapply(
      X = tx,
      FUN = gradient_vector,
      theta = theta,
      q = q,
      k = k
    )
    b <- (1 / n) * Reduce(
      "+",
      gradient
    )
    hessian <- lapply(
      X = tx,
      FUN = hessian_matrix,
      theta = theta,
      q = q,
      k = k
    )
    a <- (-1 / n) * Reduce(
      "+",
      hessian
    )
    inva <- solve(a)
    output <- inva %*% b %*% inva
    if (names) {
      colnames(output) <- rownames(output) <- gammacapnames(
        x,
        sep = sep,
        mean_structure = TRUE
      )
    }
    if (drop_means) {
      output <- output[(k + 1):q, (k + 1):q, drop = FALSE]
    }
    return(output)
  }
}
