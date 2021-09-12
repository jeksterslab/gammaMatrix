## ---- test-linearAlgebra-vech
x_i <- matrix(
  data = c(
    1.0, 0.5, 0.4,
    0.5, 1.0, 0.6,
    0.4, 0.6, 1.0
  ),
  ncol = 3
)
answer_i <- c(1.0, 0.5, 0.4, 1.0, 0.6, 1.0)
result_i <- vech(x_i, names = FALSE)
testthat::test_that("linearAlgebra-vech 3 by 3", {
  testthat::expect_equal(
    result_i,
    answer_i
  )
})
x_i <- matrix(
  c("a", "b", "b", "c"),
  ncol = 2
)
answer_i <- c("a", "b", "c")
result_i <- vech(x_i, names = FALSE)
testthat::test_that("linearAlgebra-vech 2 by 2", {
  testthat::expect_equal(
    result_i,
    answer_i
  )
})
n_i <- 100
k_i <- sample(
  1:10,
  size = 1
)
mu_i <- rep(
  x = 0,
  times = k_i
)
sigmacap_i <- matrix(
  runif(
    n = 1,
    min = 0,
    max = 1
  ),
  nrow = k_i,
  ncol = k_i
)
diag(sigmacap_i) <- 1
x_i <- matrix(
  data = rnorm(
    n = n_i * k_i
  ),
  nrow = n_i,
  ncol = k_i
) %*% (
  chol(sigmacap_i)
) + (
  matrix(
    data = 1,
    nrow = n_i,
    ncol = 1
  ) %*% mu_i
)
x_i <- cov(x_i)
answer_i <- x_i[lower.tri(x_i, diag = TRUE)]
result_i <- vech(x_i, names = FALSE)
testthat::test_that("linearAlgebra-vech random cov", {
  testthat::expect_equal(
    result_i,
    answer_i
  )
})
# coverage
vech(x_i, names = TRUE)
# clean environment
rm(
  x_i,
  answer_i,
  result_i,
  n_i,
  k_i,
  mu_i,
  sigmacap_i
)
