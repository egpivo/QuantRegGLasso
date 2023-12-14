# test_qrglasso.R

# Load required libraries
library(testthat)
library(QuantRegGLasso)  # Replace with your actual package name
tol <- 1e-4

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

test_that("qrglasso with omega", {
  # Generate some example data
  set.seed(123)
  n <- 100
  p <- 100
  L <- 3
  Y <- matrix(rnorm(n), n, 1)
  W <- matrix(runif(n * p, min = 0, max = 5), nrow = n)
  omega <- matrix(runif(p * L), ncol = L)
  
  # Call qrglasso with omega
  result <- qrglasso(
    Y = Y,
    W = W,
    L = L,
    omega = omega,
    tau = 0.7,
    qn = 1.5,
    lambda = c(0.01, 0.1, 1),
    maxit = 500,
    thr = 1e-05
  )
  
  # Perform assertions based on your expectations
  expect_true(is.list(result))
  expect_equal(dim(result$BIC)[1], 3)
  expect_lte(min(result$phi[,3]) + 0.1345752,  tol)
  expect_lte(min(result$xi[,3]) + 1.99916,  tol)
  expect_lte(min(result$gamma[,3]) + 0.1345752, tol)
})
