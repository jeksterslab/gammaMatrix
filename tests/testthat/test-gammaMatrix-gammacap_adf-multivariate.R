## ---- test-gammaMatrix-gammacap_adf-multivariate
tol_i <- 0.5
x_i <- rmvn_chol(
  n = 1000,
  mu = c(0.00, 0.00),
  sigmacap = matrix(
    data = c(1.00, 0.5, 0.5, 1),
    nrow = 2
  )
)
normal_i <- gammacap(x_i, type = "mvn")
testthat::test_that("test-gammaMatrix-gammacap_adf-multivariate unbiased = TRUE", {
  testthat::expect_true(
    all(
      abs(normal_i - gammacap(x_i, type = "adf", unbiased = TRUE)) <= tol_i
    )
  )
})
testthat::test_that("test-gammaMatrix-gammacap_adf-multivariate unbiased = FALSE", {
  testthat::expect_true(
    all(
      abs(normal_i - gammacap(x_i, type = "adf", unbiased = FALSE)) <= tol_i
    )
  )
})
# clean environment
rm(
  x_i,
  tol_i,
  normal_i
)
