# test_qrglasso.R

# Load required libraries
library(testthat)
library(QuantRegGLasso)  # Replace with your actual package name

# Test case for qrglasso function
test_that("qrglasso returns expected results", {
  # Create sample data for testing
  set.seed(123)
  n <- 100
  pL <- 10
  Y <- matrix(rnorm(n), n, 1)
  W <- matrix(rnorm(n * pL), n, pL)
  
  # Call the qrglasso function
  result <- qrglasso(Y = Y, W = W, L = 2)
  
  # Perform assertions
  expect_s3_class(result, "qrglasso")
  expect_true(all(names(result) %in% c("gamma", "xi", "phi", "BIC", "lambda", "L", "omega")))
})
