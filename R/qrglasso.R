#' @title Adaptively Weighted Group Lasso
#'
#' @description
#' The function `qrglasso` performs Adaptively Weighted Group Lasso for semiparametric quantile regression models. It estimates the coefficients of a quantile regression model with adaptively weighted group lasso regularization. The algorithm supports the use of B-spline basis functions to model the relationship between covariates and the response variable. Regularization is applied across different groups of covariates, and an adaptive weighting scheme is employed to enhance variable selection.
#'
#' @param Y A \eqn{n \times 1} data matrix where \eqn{n} is the sample size.
#' @param W A \eqn{n \times (p \times L)} B-spline matrix where \eqn{L} is the number of groups and \eqn{p} is the number of covariates.
#' @param p A numeric indicating the number of covariates.
#' @param omega A \eqn{p \times 1} weight matrix. Default value is NULL.
#' @param tau A numeric quantile of interest. Default value is 0.5.
#' @param qn A numeric bound parameter for HDIC. Default value is 1.
#' @param lambda A sequence of tuning parameters. Default value is NULL.
#' @param maxit The maximum number of iterations. Default value is 1000.
#' @param thr Threshold for convergence. Default value is \eqn{10^{-4}}.
#' @return A list with the following components:
#'   \item{\code{gamma}}{A target estimate.}
#'   \item{\code{xi}}{An auxiliary estimate in the ADMM algorithm.}
#'   \item{\code{phi}}{An auxiliary estimate in the ADMM algorithm.}
#'   \item{\code{BIC}}{A sequence of BIC values with respect to different lambdas.}
#'   \item{\code{lambda}}{A sequence of tuning parameters used in the algorithm.}
#'   \item{\code{L}}{The number of groups.}
#'   \item{\code{omega}}{A \eqn{p \times 1} weight matrix used in the algorithm.}
#' @author Wen-Ting Wang
#' @references
#' Toshio Honda, Ching-Kang Ing, Wei-Ying Wu (2019). Adaptively weighted group Lasso for semiparametric quantile regression models. \emph{Bernoulli} \bold{225} 4B.
#'
#' @export
#' 
#' @examples
#' # Example: One true non-linear covariate function
#' # Define the function g1
#' g1 <- function(x) { 
#'   (3 * sin(2 * pi * x) / (2 - sin(2 * pi * x))) - 0.4641016 
#' }
#' 
#' # Set parameters
#' n <- 100
#' p <- 100
#' err_sd <- 0.1 ** 2
#' tau <- 0.7
#' 
#' # Generate synthetic data
#' set.seed(1234)
#' x <- matrix(runif(n * p, min = 0, max = 1), n, p)
#' error_tau <- rnorm(n, sd = err_sd) - qnorm(tau, sd = err_sd)
#' y <- g1(x[, 1]) + error_tau
#' y <- y - mean(y)
#' 
#' # B-spline parameters
#' total_knots <- 5
#' degree <- 2
#' boundaries <- c(0, 1)
#' xx <- seq(from = 0, to = 1, length.out = total_knots)
#' knots <- xx[2:(total_knots - 1)]
#' 
#' # Create B-spline matrix W
#' L <- total_knots + degree - 1
#' W <- matrix(0, nrow = n, ncol = p * (L - 1))
#' 
#' for (i in 1:n) {
#'   bspline_result <- orthogonize_bspline(knots, boundaries, degree, x[i, ])
#'   W[i, ] <- matrix(t(sqrt(L) * bspline_result$bsplines[, -1]), ncol = p * (L - 1), nrow = 1)
#' }
#' 
#' # Perform quantile regression with group Lasso
#' result <- qrglasso(as.matrix(y), W, p)
#' # BIC Results
#' plot(result)
#' # Prediction
#' estimate = predict(result)
#' plot(estimate)
#' 
#' @export
qrglasso <- function(Y,
                     W,
                     p,
                     omega = NULL,
                     tau = 0.5,
                     qn = 1,
                     lambda = NULL,
                     maxit = 1000,
                     thr = 1e-04) {
  if (is.null(lambda)) {
    nlambda <- 51
    max.lambda <- 10
    lambda <- c(0, exp(seq(log(max.lambda / 1e4), log(max.lambda), length = (nlambda - 1))))
  } else {
    nlambda <- length(lambda)
  }
  
  zeta <- 10
  zetaincre <- 1
  L_star <- dim(W)[2] / p
  
  if (is.null(omega))
    result <- awgl(Y, W, lambda, tau, L_star, qn, zeta, zetaincre, maxit, thr)
  else
    result <- awgl_omega(Y, W, omega, lambda, tau, qn, zeta, zetaincre, maxit, thr)
  
  result$phi[, 1] <- result$gamma[, 1]
  
  obj.cv <- list(
    gamma = result$gamma,
    xi = result$xi,
    phi = result$phi,
    BIC = result$BIC,
    lambda = lambda,
    L = L_star + 1,
    omega = result$omega
  )
  
  class(obj.cv) <- "qrglasso"
  return(obj.cv)
}


#' @title Predict Top-k Coefficient Functions
#'
#' @description Predict the top-k coefficient functions based on a \code{qrglasso} class object.
#'
#' @param qrglasso_object A \code{qrglasso} class object.
#' @param metric_type Character. Metric type for gamma selection, e.g., `BIC`, `BIC-log`. Default is `BIC`.
#' @param top_k Integer. The number of top estimated functions to predict. Default is 5.
#' @param degree Integer. Degree of the piecewise polynomial. Default is 2.
#' @param boundaries Array. Two boundary points for the piecewise polynomial. Default is c(0, 1).
#' @param is_approx Logical. If TRUE, the size of covariate indexes will be 1e6; otherwise, 1e4. Default is FALSE.
#' @seealso \code{\link{qrglasso}}
#' @return A list containing:
#'   \item{\code{coef_functions}}{Matrix. The estimated top-k coefficient functions with dimension (\eqn{m \times k}), where \eqn{m} is the size of \code{z}.}
#'   \item{\code{z}}{Array. Index predictors used in the generation.}
#' @examples
#' set.seed(123)
#' n <- 100
#' p <- 5
#' L <- 5
#' Y <- matrix(rnorm(n), n, 1)
#' W <- matrix(rnorm(n * p * (L - 1)), n, p * (L - 1))
#'
#' # Call qrglasso with default parameters
#' result <- qrglasso(Y = Y, W = W, p = p)
#' estimate <- predict(result) 
#' print(dim(estimate$coef_functions))
#' 
#' @export
predict <- function(qrglasso_object,
                    metric_type = "BIC",
                    top_k = 5,
                    degree = 2,
                    boundaries = c(0, 1),
                    is_approx = FALSE) {
  check_predict_parameters(qrglasso_object, metric_type, top_k, degree, boundaries)
  total_knots <- qrglasso_object$L - degree + 1
  x <- seq(boundaries[1], boundaries[2], length.out = total_knots)
  knots <- x[2:(total_knots - 1)]
  approx_bsplines <- orthogonize_bspline(knots, boundaries, degree, is_approx = is_approx)
  bsplines <- approx_bsplines$bsplines[, -1]
  z <- approx_bsplines$z
  if (metric_type == "BIC") {
    gamma_hat <- qrglasso_object$gamma[, which.min(qrglasso_object$BIC[, 1])] 
  } else {
    gamma_hat <- qrglasso_object$gamma[, which.min(qrglasso_object$BIC[, 2])]
  }
  
  estimate <- bsplines %*% matrix(gamma_hat, nrow = dim(bsplines)[2])
  obj.predict <- list(
    coef_functions = as.matrix(estimate[, 1:min(top_k, dim(estimate)[2])]),
    z = z
  )
  class(obj.predict) <- "qrglasso.predict"
  return(obj.predict)
}

#' @title Display BIC Results from `qrglasso`
#'
#' @description Visualize the HDIC BIC results corresponding to hyperparameters obtained from `qrglasso`.
#'
#' @param x An object of class \code{qrglasso} for the \code{plot} method.
#' @param ... Additional parameters not used directly.
#' @return \code{NULL}
#' @seealso \code{\link{qrglasso}}
#'
#' @export
#' @method plot qrglasso
#' @examples
#' set.seed(123)
#' n <- 100
#' p <- 5
#' L <- 5
#' Y <- matrix(rnorm(n), n, 1)
#' W <- matrix(rnorm(n * p * (L - 1)), n, p * (L - 1))
#'
#' # Call qrglasso with default parameters
#' result <- qrglasso(Y = Y, W = W, p = p)
#' 
#' # Visualize the BIC results
#' plot(result)
#' 
#' @export
plot.qrglasso <- function(x, ...) {
  if (!inherits(x, "qrglasso")) {
    stop("Invalid object! Please enter a `qrglasso` object.")
  }
  originalPar <- par(no.readonly = TRUE)
  result <- list()
  variates <- c("BIC", "BIC-log")
  for (i in 1:2) {
    variate <- variates[i]
    data <- data.frame(lambda = x$lambda, bic = x$BIC[, i])
    result[[variate]] <- plot_bic_result(data, variate)
  }
  plot_sequentially(result)
  par(originalPar)
}

#' @title Display Predicted Coefficient Functions from `qrglasso`
#'
#' @description Visualize the predicted coefficient functions selected by BIC.
#'
#' @param x An object of class \code{qrglasso.predict} for the \code{plot} method.
#' @param ... Additional parameters not used directly.
#' @return \code{NULL}
#' @seealso \code{\link{qrglasso}}
#'
#' @export
#' @method plot qrglasso.predict
#' @examples
#' set.seed(123)
#' n <- 100
#' p <- 5
#' L <- 5
#' Y <- matrix(rnorm(n), n, 1)
#' W <- matrix(rnorm(n * p * (L - 1)), n, p * (L - 1))
#'
#' # Call qrglasso with default parameters
#' result <- qrglasso(Y = Y, W = W, p = p)
#' 
#' # Predict the top-k coefficient functions
#' estimate <- predict(result, top_k = 2)
#' 
#' # Display the predicted coefficient functions
#' plot(estimate)
#'
#' @export
plot.qrglasso.predict <- function(x, ...) {
  if (!inherits(x, "qrglasso.predict")) {
    stop("Invalid object! Please enter a `qrglasso.predict` object.")
  }
  originalPar <- par(no.readonly = TRUE)
  k <- dim(x$coef_functions)[2]
  result <- list()
  for (i in 1:k) {
    variate <- paste0("Coefficient function - g", i)
    data <- data.frame(z = x$z, coefficient = x$coef_functions[, i])
    result[[variate]] <- plot_coefficient_function(data, variate)
  }
  plot_sequentially(result)
  par(originalPar)
}
