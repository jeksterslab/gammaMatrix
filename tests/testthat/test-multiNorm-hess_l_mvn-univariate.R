## ---- test-multiNorm-hess_l_mvn-univariate
tol_i <- 0.05
n_i <- 5
mu_i <- 0
sigmacap_i <- matrix(
  data = 1,
  nrow = 1
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
    numDeriv::hessian(
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
result_i <- lapply(
  X = xcap_i,
  FUN = function(i, theta) {
    hess_l_mvn(i, theta)
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
testthat::test_that("test-multiNorm-hess_l_mvn-univariate", {
  testthat::expect_true(
    all(
      abs(
        result_i - answer_i
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
  result_i
)
