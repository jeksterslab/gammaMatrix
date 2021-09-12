## ---- test-linearAlgebra-vechnames
testthat::test_that("linearAlgebra-vechnames 1 by 1", {
  testthat::expect_equal(
    "x1.x1",
    vechnames("x1")
  )
})
testthat::test_that("linearAlgebra-vechnames 2 by 2", {
  testthat::expect_equal(
    c("x1.x1", "x1.x2", "x2.x2"),
    vechnames(c("x1", "x2"))
  )
})
testthat::test_that("linearAlgebra-vechnames 3 by 3", {
  testthat::expect_equal(
    c("x1.x1", "x1.x2", "x1.x3", "x2.x2", "x2.x3", "x3.x3"),
    vechnames(c("x1", "x2", "x3"))
  )
})
