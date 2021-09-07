## ---- test-gammaMatrix-gammacap_ols_hc
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
testthat::test_that("test ols_hc", {
  testthat::expect_true(
    all(
      abs(
        normal_i - gammacap_ols_hc(obj_i)
      ) <= tol_i
    )
  )
})
# coverage
head(gammacap_ols_hc_qcap(obj_i, type = "hc0"))
head(gammacap_ols_hc_qcap(obj_i, type = "hc1"))
head(gammacap_ols_hc_qcap(obj_i, type = "hc2"))
head(gammacap_ols_hc_qcap(obj_i, type = "hc3"))
head(gammacap_ols_hc_qcap(obj_i, type = "hc4"))
head(gammacap_ols_hc_qcap(obj_i, type = "hc4m"))
head(gammacap_ols_hc_qcap(obj_i, type = "hc5"))
# gammacap_ols_hc_generic
gammacap_ols_hc_generic(
  x = obj_i$model,
  h = NULL
)
# clean environment
rm(
  x_i,
  tol_i,
  normal_i,
  obj_i
)
