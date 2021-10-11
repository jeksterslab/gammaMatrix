## ---- test-gammaMatrix-gammacap_mvnadj2-univariate
set.seed(42)
tol_i <- 0.5
x_i <- as.vector(scale(rnorm(1000)))
normal_i <- gammacap(x_i, type = "mvn")
testthat::test_that("test-gammaMatrix-gammacap_mvnadj2-univariate", {
  testthat::expect_error(
    # univariate cases not supported
    gammacap(x_i, type = "mvnadj2")
    # all(
    #  abs(normal_i - gammacap(x_i, type = "mvnadj2")) <= tol_i
    # )
  )
})
# clean environment
rm(
  x_i,
  tol_i,
  normal_i
)
