library(testthat)
library(QuantRegGLasso)
tol <- 1e-4

# Test case for qrglasso function
test_that("qrglasso returns expected results", {
  # Create sample data for testing
  set.seed(123)
  n <- 100
  p <- 2
  L <- 5
  pL <- 2 * 5
  Y <- matrix(rnorm(n), n, 1)
  W <- matrix(rnorm(n * pL), n, pL)
  
  # Call the qrglasso function
  result <- qrglasso(Y = Y, W = W, p)
  
  # Perform assertions
  expect_s3_class(result, "qrglasso")
  expected_fields <- c(
    "gamma", "xi", "phi", "BIC", "lambda", "L",
    "omega", "solver_iterations", "solver_max_error"
  )
  expect_setequal(names(result), expected_fields)
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
    p = p,
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
  expect_equal(dim(result$omega), c(p, 2))
  expect_equal(length(result$solver_iterations), length(result$lambda))
  expect_equal(length(result$solver_max_error), length(result$lambda))
})

test_that("qrglasso validates W shape", {
  set.seed(123)
  n <- 20
  Y <- matrix(rnorm(n), n, 1)
  W <- matrix(rnorm(n * 5), n, 5)
  expect_error(qrglasso(Y = Y, W = W, p = 3), "Number of columns in W must be divisible by p")
})

test_that("lambda sequence is sorted when unsorted provided", {
  set.seed(123)
  n <- 50
  p <- 2
  L <- 4
  Y <- matrix(rnorm(n), n, 1)
  W <- matrix(rnorm(n * p * L), n, p * L)
  expect_warning(
    result <- qrglasso(Y = Y, W = W, p = p, lambda = c(1, 0, 0.5)),
    "lambda sequence is not non-decreasing"
  )
  expect_false(is.unsorted(result$lambda))
  expect_equal(result$lambda, sort(c(1, 0, 0.5)))
})

test_that("awgl handles zero-norm groups", {
  Y <- matrix(0, 4, 1)
  W <- matrix(0, 4, 2)
  res <- QuantRegGLasso:::awgl(
    Y = Y,
    W = W,
    lambda = c(0, 0.5),
    tau = 0.5,
    L = 2,
    qn = 1,
    zeta = 1,
    zetaincre = 1,
    maxit = 10,
    tol = 1e-04
  )
  expect_true(is.finite(res$solver_max_error[2]))
  expect_false(is.nan(res$solver_max_error[2]))
})

test_that("qrglasso validates omega dimensions", {
  set.seed(1)
  n <- 30
  p <- 4
  L <- 3
  Y <- matrix(rnorm(n), n, 1)
  W <- matrix(rnorm(n * p * L), n, p * L)
  bad_omega <- matrix(runif(5 * L), nrow = 5)
  expect_error(
    qrglasso(Y = Y, W = W, p = p, omega = bad_omega),
    "omega must have p rows"
  )
})

# Mock qrglasso class object for testing
mock_qrglasso <- structure(list(
  L = 6,
  gamma = matrix(rnorm(400), nrow = 5),
  BIC = matrix(runif(10), nrow = 5),
  omega = matrix(runif(120), nrow = 6)
), class = "qrglasso")

test_that("predict coefficient functions", {
  # Valid parameters
  expect_silent(predict(mock_qrglasso))
  
  # Invalid object
  expect_error(predict(list()))
  
  # Negative top_k
  expect_error(predict(mock_qrglasso, "BIC", -2))
  
  # Negative degree
  expect_error(predict(mock_qrglasso, "BIC", 3, -1))
  
  # Incorrect boundaries size
  expect_error(predict(mock_qrglasso, "BIC", 3, 2, c(0, 1, 2)))
  
  # Invalid boundaries order
  expect_error(predict(mock_qrglasso, "BIC", 3, 2, c(1, 0)))
  
  # Invalid metric_type
  expect_error(predict(mock_qrglasso, "invalid_metric"))
})
