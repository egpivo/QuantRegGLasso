#' Adaptively weighted group Lasso
#'
#' @param Y A  \eqn{n \times 1} data matrix
#' @param W A  \eqn{n \times pL}B-spline matrix.
#' @param L The number of groups
#' @param omega A $p x 1$ weight matrix. Default value is NULL.
#' @param tau A numeric quantile of interest. Default value is 0.5.
#' @param qn  A nuneric bound parameter for HDIC.  Default value is 1.
#' @param lambda A sequence of tuning parameters. Default value is NULL.
#' @param maxit The maximum number of iterations. Default value is 1000.
#' @param thr Threshold for convergence. Default value is \eqn{10^{-4}}.
#' @return This function returns a \code{list} including:
#' \itemize{
#'  \item{gamma}{A target estimate}
#'  \item{xi}{An auxiliary estimate in the ADMM algorithm}
#'  \item{phi}{An auxiliary estimate in the ADMM algorithm}
#'  \item{BIC}{A sequence of BIC values w.r.t. different lambdas}
#'  \item{lambda}{A sequence of tuning parameters used in the algorithm.}
#'  \item{omega}{A $p x 1$ weight matrix used in the algorithm.}
#' }
#' @author Wen-Ting Wang
#' @references Toshio Honda, Ching-Kang Ing, Wei-Ying Wu (2019). Adaptively weighted group Lasso for semiparametric quantile regression models. \emph{Bernoulli} \bold{225} 4B.
#' @export
#' @examples
#' # Example: One true non-linear covariate function
#' # Define the function g1
#' g1 <- function(x) { 
#'   (3 * sin(2 * pi * x) / (2 - sin(2 * pi * x))) - 0.4641016 
#' }
#' 
#' # Set parameters
#' n <- 150
#' p <- 150
#' err_sd <- 0.1 ** 2
#' tau <- 0.9
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
#' result <- qrglasso(as.matrix(y), W, L - 1)
#' 
#' 
#' # Extract relevant information from the result
#' approx_bsplines <- orthogonize_bspline(knots, boundaries, degree)
#' bsplines <- approx_bsplines$bsplines[, -1]
#' z <- approx_bsplines$z
#' L_star <- L - 1
#' gamma_hat <- result$gamma[, which.min(result$BIC[, 1])]
#' g1_BIC <- bsplines %*% gamma_hat[1:L_star]
#' g2_BIC <- bsplines %*% gamma_hat[(L_star + 1):(2 * L_star)]
#' 
#' g10 <- bsplines %*% result$gamma[, 1][1:L_star]
#' g20 <- bsplines %*% result$gamma[, 1][(L_star + 1):(2 * L_star)]
#' 
#' # Plotting
#' original_par <- par(no.readonly = TRUE)
#' par(mfrow = c(1, 2))
#' 
#' # Plot for g1(z)
#' plot(z, g1(z), type = 'l', lwd = 4, cex.lab = 2, lty = 2, cex.axis = 1.5, ylim = c(-3, 3))
#' lines(z, g10, col = 2, lwd = 2)
#' lines(z, g1_BIC, col = 3, lwd = 2)
#' legend("topleft", c("True", expression(lambda == 0), expression(lambda[BIC])), col = c(1, 2, 3), 
#'        lty = c(2, 1, 1), lwd = c(3, 3, 3), horiz = TRUE, bty = "n")
#'
#' # Plot for g2(z)
#' plot(z, rep(0, length(z)), type = 'l', lwd = 4, cex.lab = 2, lty = 2, cex.axis = 1.5, 
#'      ylab = "g2(z)", ylim = c(-0.1, 0.1))
#' lines(z, g20, col = 2, lwd = 2)
#' lines(z, g2_BIC, col = 3, lwd = 2)
#' legend("topleft", c("True", expression(lambda == 0), expression(lambda[BIC])), col = c(1, 2, 3, 4), 
#'        lty = c(2, 1, 1, 1), lwd = c(3, 3, 3, 3), horiz = TRUE, bty = "n")
#'        
#' par(original_par)
#' 
qrglasso <-
  function(Y,
           W,
           L,
           omega = NULL,
           tau = 0.5,
           qn = 1,
           lambda = NULL,
           maxit = 1000,
           thr = 1e-04) {
    if (is.null(lambda)) {
      nlambda <- 51
      max.lambda <- 10
      lambda <-
        c(0, exp(seq(
          log(max.lambda / 1e4), log(max.lambda), length = (nlambda - 1)
        )))
    } else {
      nlambda <- length(lambda)
    }
    
    zeta <- 10
    zetaincre <- 1
    
    if (is.null(omega))
      result <-
      awgl(Y, W, lambda, tau, L, qn, zeta, zetaincre, maxit, thr)
    else
      result <-
      awgl_omega(Y, W, omega, lambda, tau, qn, zeta, zetaincre, maxit, thr)
    
    result$phi[, 1] <- result$gamma[, 1]
    
    obj.cv <- list(
      gamma = result$gamma,
      xi = result$xi,
      phi = result$phi,
      BIC = result$BIC,
      lambda = lambda,
      L = L,
      omega = result$omega
    )
    
    class(obj.cv) <- "qrglasso"
    return(obj.cv)
  }


#' @title Predict the coefficient functions
#'
#' @description Predict the top-k coefficient functions
#'
#' @param qrglasso_object An `qrglasso` class object.
#' @param top_k Integer. A matrix of the top K estimated functions. Default is 5.
#' @param degree Integer. Degree of the piecewise polynomial. Default is 2.
#' @param boundaries Array. Two boundary points. Default is c(0, 1).
#' @param is_approx Logical. If TRUE, the size of covariate indexes will be 1e6; otherwise, 1e4. Default is FALSE.
#' @seealso \link{qrglasso}
#' @return A list containing:
#'   \item{coef_functions}{Matrix. Top-k coefficient function estimates with dimenstion (\eqn{m \times k}) where $m$ is size of `z`.}
#'   \item{z}{Array. Index predictors used in generation}
#' @examples
#' set.seed(123)
#' n <- 100
#' p <- 5
#' L <- 5
#' Y <- matrix(rnorm(n), n, 1)
#' W <- matrix(rnorm(n * p * (L - 1)), n, p * (L - 1))
#'
#' # Call qrglasso with default parameters
#' result <- qrglasso(Y = Y, W = W, L = 5)
#' estimate <- predict(result) 
#' print(dim(estimate$coef_functions))
#' 
predict <-
  function(qrglasso_object,
           top_k = 5,
           degree = 2,
           boundaries = c(0, 1),
           is_approx = FALSE) {
    check_predict_parameters(qrglasso_object, top_k, degree, boundaries)
    total_knots <- qrglasso_object$L - degree + 1
    x <- seq(boundaries[1],  boundaries[2], length.out = total_knots)
    knots <- x[2:(total_knots - 1)]
    approx_bsplines <-orthogonize_bspline(knots, boundaries, degree, is_approx=is_approx)
    bsplines <- approx_bsplines$bsplines[,-1]
    z <- approx_bsplines$z
    gamma_hat <- qrglasso_object$gamma[, which.min(qrglasso_object$BIC[, 1])]
    estimate <- bsplines %*%  matrix(gamma_hat, nrow = dim(bsplines)[2])
    obj.predict <- list(
      coef_functions = as.matrix(estimate[, 1:min(top_k, dim(estimate)[2])]),
      z = z
    )
    class(obj.predict) <- "qrglasso.predict"
    return(obj.predict)
  }

#' @title  Display the estimated coefficient functions
#'
#' @description Display the estimated coefficient functions by BIC
#'
#' @param x An qrglasso.predict class object for `plot` method
#' @param ... Not used directly
#' @return `NULL`.
#' @seealso \link{qrglasso}
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
#' result <- qrglasso(Y = Y, W = W, L = 5)
#' estimate <- predict(result, top_k = 2)
#' plot(estimate)
#'
plot.qrglasso.predict <- function(x, ...) {
  if (!inherits(x, "qrglasso.predict")) {
    stop("Invalid object! Please enter a `qrglasso.predict` object")
  }
  originalPar <- par(no.readonly = TRUE)
  k <- dim(x$coef_functions)[2]
  result <- list()
  for (i in 1:k) {
    variate <- paste0("Coefficient function - g", i)
    data <- data.frame(z = x$z, coef = x$coef_functions[,i])
    result[[variate]] <- plot_coefficient_function(data, variate)
  }
  plot_sequentially(result)
  par(originalPar)
}
