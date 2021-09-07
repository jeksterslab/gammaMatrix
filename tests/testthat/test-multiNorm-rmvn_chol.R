## ---- test-multiNorm-rmvn_chol
tol_i <- 0.01
k_i <- sample(x = 2:10, size = 1)
data_i <- rmvn_chol(
  n = 1000000,
  mu = rep(x = 0, times = k_i),
  sigmacap = toeplitz((k_i:1) / k_i)
)
testthat::test_that("means", {
  testthat::expect_true(
    all(
      abs(
        colMeans(data_i) - 0
      ) <= tol_i
    )
  )
})
testthat::test_that("covariances", {
  testthat::expect_true(
    all(
      abs(
        as.vector(
          cov(data_i)
        ) - as.vector(
          toeplitz((k_i:1) / k_i)
        )
      ) <= tol_i
    )
  )
})
# coverage
data_i <- rmvn_chol(
  n = 1000000,
  mu = rep(x = 0, times = k_i),
  sigmacap = toeplitz((k_i:1) / k_i),
  varnames = paste0("x", seq_len(k_i)),
  data_frame = TRUE
)
# clean environment
rm(
  tol_i,
  k_i,
  data_i
)
