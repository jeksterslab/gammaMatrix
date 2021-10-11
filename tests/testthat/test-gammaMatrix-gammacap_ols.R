## ---- test-gammaMatrix-gammacap_ols
set.seed(42)
tol_i <- 0.5
x_i <- rmvn_chol(
  n = 1000,
  mu = c(0.00, 0.00),
  sigmacap = matrix(
    data = c(1.00, 0.5, 0.5, 1),
    nrow = 2
  ),
  varnames = c("y", "x"),
  data_frame = TRUE
)
normal_i <- gammacap(x_i, type = "mvn")
obj_i <- lm(
  y ~ x,
  data = x_i
)
testthat::test_that("test-gammaMatrix-gammacap_ols", {
  testthat::expect_true(
    all(
      abs(
        normal_i - gammacap_ols(obj_i)
      ) <= tol_i
    )
  )
})
# coverage
gammacap_ols(obj_i, yc = TRUE)
gammacap_ols(obj_i, yc = FALSE)
gammacap_ols(obj_i, ke_unbiased = TRUE)
gammacap_ols(obj_i, ke_unbiased = FALSE)
# gammacap_ols_generic
gammacap_ols_generic(
  x = as.data.frame(obj_i$model[, -1, drop = FALSE]),
  beta = stats::coef(obj_i)[-1],
  sigmacap = unname(as.matrix(stats::cov(obj_i$model))),
  sigmasq = stats::summary.lm(obj_i)$sigma^2,
  ke = mean(obj_i$residuals^4) / mean(obj_i$residuals^2)^2
)
# clean environment
rm(
  x_i,
  tol_i,
  normal_i,
  obj_i
)
