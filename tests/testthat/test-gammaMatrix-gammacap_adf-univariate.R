## ---- test-gammaMatrix-gammacap_adf-univariate
tol_i <- 0.5
x_i <- as.vector(scale(rnorm(1000)))
normal_i <- gammacap(x_i, type = "mvn")
testthat::test_that("test univariate adf unbiased = TRUE", {
  testthat::expect_true(
    all(
      abs(normal_i - gammacap(x_i, type = "adf", unbiased = TRUE)) <= tol_i
    )
  )
})
testthat::test_that("test univariate adf unbiased = FALSE", {
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
