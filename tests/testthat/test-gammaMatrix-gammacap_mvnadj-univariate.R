## ---- test-gammaMatrix-gammacap_mvnadj-univariate
tol_i <- 0.5
x_i <- as.vector(scale(rnorm(1000)))
normal_i <- gammacap(x_i, type = "mvn")
testthat::test_that("test univariate mvnadj", {
  testthat::expect_true(
    all(
      abs(normal_i - gammacap(x_i, type = "mvnadj")) <= tol_i
    )
  )
})
# clean environment
rm(
  x_i,
  tol_i,
  normal_i
)
