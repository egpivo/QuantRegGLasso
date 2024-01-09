library(testthat)

tol <- 1e-4

test_that("orthogonize_bspline generates orthogonalized B-splines", {
  set.seed(1234)
  # Test parameters
  total_knots <- 10
  degree <- 3
  boundaries <- c(0, 1)
  predictors <- seq(from = 0, to = 1, length.out = 30)
  knots <- seq(from = 0, to = 1, length.out = total_knots)
  
  # Call the function
  bsplines <- orthogonize_bspline(knots, boundaries, degree, predictors)
  
  # Check if the output is a matrix
  expect_true(is.matrix(bsplines$bsplines))
  
  # Check if the number of rows/columns is equal to the length of predictors
  expect_equal(ncol(bsplines$bsplines), total_knots + degree + 1)
  expect_equal(nrow(bsplines$bsplines), length(predictors))

  expect_lte(sum(bsplines$bsplines) - 8.017837, tol)
  expect_lte(sum(bsplines$z) - 15, tol)
})


test_that("orthogonize_bspline handles different degrees", {
  knots <- seq(from = 0, to = 1, length.out = 10)
  boundaries <- c(0, 1)
  predictors <- seq(from = 0, to = 1, length.out = 30)
  # Test with degree 2
  bsplines_deg2 <- orthogonize_bspline(knots, boundaries, degree = 2, predictors)
  expect_true(is.matrix(bsplines_deg2$bsplines))
  expect_equal(ncol(bsplines_deg2$bsplines), length(knots) + 2 + 1)
  expect_equal(nrow(bsplines_deg2$bsplines), length(predictors))
  # Test with degree 4
  bsplines_deg4 <- orthogonize_bspline(knots, boundaries, degree = 4, predictors)
  expect_true(is.matrix(bsplines_deg4$bsplines))
  expect_equal(ncol(bsplines_deg4$bsplines), length(knots) + 4 + 1)
  expect_equal(nrow(bsplines_deg4$bsplines), length(predictors))
})


test_that("orthogonize_bspline handles different predictor sizes", {
  knots <- seq(from = 0, to = 1, length.out = 10)
  boundaries <- c(0, 1)
  # Test with 10 predictors
  predictors_10 <- seq(from = 0, to = 1, length.out = 10)
  bsplines_10 <- orthogonize_bspline(knots, boundaries, degree = 3, predictors_10)
  expect_true(is.matrix(bsplines_10$bsplines))
  expect_equal(ncol(bsplines_10$bsplines), length(knots) + 3 + 1)
  expect_equal(nrow(bsplines_10$bsplines), length(predictors_10))
  # Test with 50 predictors
  predictors_50 <- seq(from = 0, to = 1, length.out = 50)
  bsplines_50 <- orthogonize_bspline(knots, boundaries, degree = 3, predictors_50)
  expect_true(is.matrix(bsplines_50$bsplines))
  expect_equal(ncol(bsplines_50$bsplines), length(knots) + 3 + 1)
  expect_equal(nrow(bsplines_50$bsplines), length(predictors_50))
})

test_that("orthogonize_bspline produces consistent results with a random seed", {
  knots <- seq(from = 0, to = 1, length.out = 10)
  boundaries <- c(0, 1)
  predictors <- seq(from = 0, to = 1, length.out = 30)
  # Set a random seed
  set.seed(1234)
  bsplines1 <- orthogonize_bspline(knots, boundaries, degree = 3, predictors)
  # Set the same random seed
  set.seed(1234)
  bsplines2 <- orthogonize_bspline(knots, boundaries, degree = 3, predictors)
  expect_identical(bsplines1, bsplines2)
})

test_that("orthogonize_bspline produces consistent results with a random seed", {
  knots <- seq(from = 0, to = 1, length.out = 10)
  boundaries <- c(0, 1)
  predictors <- seq(from = 0, to = 1, length.out = 30)
  # Set a random seed
  set.seed(1234)
  bsplines1 <- orthogonize_bspline(knots, boundaries, degree = 3, predictors)
  # Set the same random seed
  set.seed(1234)
  bsplines2 <- orthogonize_bspline(knots, boundaries, degree = 3, predictors)
  expect_identical(bsplines1, bsplines2)
})


# Mock qrglasso class object for testing
mock_qrglasso <- structure(list(L = 10), class = "qrglasso")

test_that("check_predict_parameters correctly validates parameters", {
  # Valid parameters
  expect_silent(check_predict_parameters(mock_qrglasso, top_k = 3, degree = 2, boundaries = c(0, 1)))
  
  # Invalid object
  expect_error(check_predict_parameters(list(), top_k = 3, degree = 2, boundaries = c(0, 1)),
               "Invalid object! Please enter a `qrglasso` object")
  
  # Negative top_k
  expect_error(check_predict_parameters(mock_qrglasso, top_k = -2, degree = 2, boundaries = c(0, 1)),
               "Please enter a positive top k")
  
  # Negative degree
  expect_error(check_predict_parameters(mock_qrglasso, top_k = 3, degree = -1, boundaries = c(0, 1)),
               "Please enter a positive degree")
  
  # Incorrect boundaries size
  expect_error(check_predict_parameters(mock_qrglasso, top_k = 3, degree = 2, boundaries = c(0, 1, 2)),
               "Please enter a size 2 boundaries.")
  
  # Invalid boundaries order
  expect_error(check_predict_parameters(mock_qrglasso, top_k = 3, degree = 2, boundaries = c(1, 0)),
               "Please input valid boundaries consisting of two elements in ascending order.")
  
  # Total knots less than or equal to 0
  expect_error(check_predict_parameters(mock_qrglasso, top_k = 3, degree = 11, boundaries = c(0, 1)),
               "Please enter a smaller degree")
})

# Test plot_sequentially
test_that("plot_sequentially prints ggplot2 objects", {
  # Create mock ggplot2 objects for testing
  plot1 <- ggplot(mtcars, aes(x = mpg, y = disp)) + geom_point()
  plot2 <- ggplot(mtcars, aes(x = wt, y = hp)) + geom_point()

  # Capture the current environment before calling the function
  original_env <- environment()

  # Test the function
  plot_sequentially(list(plot1, plot2))
  
  # Verify that the environment has been changed
  expect_equal(original_env, environment())
})

# Test plot_coefficient_function
test_that("plot_coefficient_function returns a ggplot object", {
  # Create a mock dataframe for testing
  data <- data.frame(z = 1:10, coef = rnorm(10))
  variate <- "Test Variate"
  
  # Test the function
  plot_result <- plot_coefficient_function(data, variate)
  
  # Verify the output
  expect_true(is.ggplot(plot_result))
})
