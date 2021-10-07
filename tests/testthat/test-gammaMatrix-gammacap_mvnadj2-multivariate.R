## ---- test-gammaMatrix-gammacap_mvnadj2-multivariate
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
sigmacap_i <- cov(x_i) * (nrow(x_i) - 1) / nrow(x_i)
testthat::test_that("test-gammaMatrix-gammacap_mvnadj2-multivariate", {
  testthat::expect_true(
    all(
      abs(normal_i - gammacap(x_i, type = "mvnadj2")) <= tol_i
    )
  )
})
# means
mu_i <- gammacap_mvnadj2(
  x_i,
  ml_cov = TRUE,
  drop_means = FALSE
)
testthat::test_that("test-gammaMatrix-gammacap_mvnadj2-multivariate means", {
  testthat::expect_true(
    all(
      as.vector(
        sigmacap_i
      ) - as.vector(
        mu_i[(1:2), (1:2)]
      ) <= tol_i
    )
  )
})
# missing
x_i[sample(dim(x_i)[1], size = 5), 1] <- NA
x_i[sample(dim(x_i)[1], size = 5), 2] <- NA
# testthat::test_that("test-gammaMatrix-gammacap_mvnadj2-multivariate missing", {
#  testthat::expect_true(
#    all(
#      abs(normal_i - gammacap(x_i, type = "mvnadj2", missing = TRUE)) <= tol_i
#    )
#  )
# })
## means
# mu_i <- gammacap_mvnadj2(
#  x_i,
#  missing = TRUE,
#  ml_cov = TRUE,
#  drop_means = FALSE
# )
# testthat::test_that("test-gammaMatrix-gammacap_mvnadj2-multivariate missing means", {
#  testthat::expect_true(
#    all(
#      as.vector(
#        sigmacap_i
#      ) - as.vector(
#        mu_i[(1:2), (1:2)]
#      ) <= tol_i
#    )
#  )
# })
# clean environment
rm(
  x_i,
  tol_i,
  normal_i,
  mu_i
)
