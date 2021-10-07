## ---- test-multiNorm-theta_helper
k_i <- sample(x = 2:10, size = 1)
vech_i <- toeplitz((k_i:1) / k_i)
vech_i <- vech_i[lower.tri(vech_i, diag = TRUE)]
theta_i <- theta_helper(
  c(
    rep(x = 0, times = k_i),
    vech_i
  )
)
testthat::test_that("test-multiNorm-theta_helper means", {
  testthat::expect_equal(
    theta_i$mu,
    rep(x = 0, times = k_i)
  )
})
testthat::test_that("test-multiNorm-theta_helper covariances", {
  testthat::expect_equal(
    theta_i$sigmacap,
    toeplitz((k_i:1) / k_i)
  )
})
# expect_error
testthat::test_that("test-multiNorm-theta_helper error", {
  testthat::expect_error(
    theta_helper(1)
  )
})
# clean environment
rm(
  k_i,
  theta_i,
  vech_i
)
