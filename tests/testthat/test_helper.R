library(testthat)

tol <- 1e-4

test_that("orth_bspline generates orthogonalized B-splines", {
  set.seed(1234)
  # Test parameters
  total_knots <- 10
  degree <- 3
  boundaries <- c(0, 1)
  predictors <- seq(from = 0, to = 1, length.out = 30)
  knots <- seq(from = 0, to = 1, length.out = total_knots)
  
  # Call the function
  bsplines <- orthgonize_bspline(knots, boundaries, degree, predictors)
  
  # Check if the output is a matrix
  expect_true(is.matrix(bsplines$bsplines))
  
  # Check if the number of rows/columns is equal to the length of predictors
  expect_equal(ncol(bsplines$bsplines), total_knots + degree + 1)
  expect_equal(nrow(bsplines$bsplines), length(predictors))

  expect_lte(sum(bsplines$bsplines) - 8.017837, tol)
  expect_lte(sum(bsplines$z) - 15, tol)
})
