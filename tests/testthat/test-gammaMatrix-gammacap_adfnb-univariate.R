## ---- test-gammaMatrix-gammacap_adfnb-univariate
tol_i <- 0.5
x_i <- as.vector(scale(rnorm(1000)))
normal_i <- gammacap(x_i, type = "mvn")
testthat::test_that("test-gammaMatrix-gammacap_adfnb-univariate", {
  testthat::expect_true(
    all(
      abs(normal_i - gammacap(x_i, type = "adfnb")) <= tol_i
    )
  )
})
testthat::test_that("test-gammaMatrix-gammacap_adfnb-univariate seed", {
  testthat::expect_true(
    all(
      abs(normal_i - gammacap(x_i, type = "adfnb", seed = 42)) <= tol_i
    )
  )
})
# clean environment
rm(
  x_i,
  tol_i,
  normal_i
)
