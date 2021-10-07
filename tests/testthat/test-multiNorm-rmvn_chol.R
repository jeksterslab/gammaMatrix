## ---- test-multiNorm-rmvn_chol
tol_i <- 0.05
k_i <- sample(x = 2:10, size = 1)
data_i <- rmvn_chol(
  n = 10000,
  mu = rep(x = 0, times = k_i),
  sigmacap = toeplitz((k_i:1) / k_i)
)
testthat::test_that("test-multiNorm-rmvn_chol means", {
  testthat::expect_true(
    all(
      abs(
        round(
          colMeans(data_i),
          digits = 0
        ) - 0
      ) <= tol_i
    )
  )
})
testthat::test_that("test-multiNorm-rmvn_chol covariances", {
  testthat::expect_true(
    all(
      abs(
        round(
          as.vector(
            cov(data_i)
          ),
          digits = 2
        ) - round(
          as.vector(
            toeplitz((k_i:1) / k_i)
          ),
          digits = 2
        )
      ) <= tol_i
    )
  )
})
# coverage
data_i <- rmvn_chol(
  n = 10,
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
