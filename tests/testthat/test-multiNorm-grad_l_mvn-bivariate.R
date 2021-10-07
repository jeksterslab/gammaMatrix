## ---- test-multiNorm-grad_l_mvn-bivariate
tol_i <- 0.05
n_i <- 5
mu_i <- c(0, 0)
sigmacap_i <- matrix(
  data = c(
    1, 0.5, 0.5, 1
  ),
  nrow = 2
)
xcap_i <- as.data.frame(
  t(
    rmvn_chol(
      n = n_i,
      mu = mu_i,
      sigmacap = sigmacap_i
    )
  )
)
theta_i <- c(
  mu_i,
  vech(sigmacap_i)
)
answer_i <- lapply(
  X = xcap_i,
  FUN = function(i, theta) {
    foo <- function(theta, data) {
      l_mvn(
        theta = theta,
        x = data
      )
    }
    numDeriv::grad(
      func = foo,
      x = theta,
      data = i
    )
  },
  theta = theta_i
)
answer_i <- (1 / n_i) * Reduce(
  "+",
  answer_i
)
answer_i <- as.vector(
  answer_i
)
answer_j <- lapply(
  X = xcap_i,
  FUN = function(i, theta) {
    foo <- function(theta, data) {
      negl_mvn(
        theta = theta,
        x = data
      )
    }
    numDeriv::grad(
      func = foo,
      x = theta,
      data = i
    )
  },
  theta = theta_i
)
answer_j <- (1 / n_i) * Reduce(
  "+",
  answer_j
)
answer_j <- as.vector(
  answer_j
)
result_i <- lapply(
  X = xcap_i,
  FUN = function(i, theta) {
    grad_l_mvn(i, theta)
  },
  theta = theta_i
)
result_i <- (1 / n_i) * Reduce(
  "+",
  result_i
)
result_i <- as.vector(
  result_i
)
testthat::test_that("test-multiNorm-grad_l_mvn-bivariate", {
  testthat::expect_true(
    all(
      abs(
        result_i - answer_i
      ) <= tol_i
    )
  )
})
testthat::test_that("test-multiNorm-grad_l_mvn-bivariate negl", {
  testthat::expect_true(
    all(
      abs(
        -1 * answer_i - answer_j
      ) <= tol_i
    )
  )
})
# clean environment
rm(
  tol_i,
  n_i,
  mu_i,
  sigmacap_i,
  xcap_i,
  theta_i,
  answer_i,
  answer_j,
  result_i
)
