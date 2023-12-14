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
#' @export
#' @examples
#' 
#' @return This function returns a \code{list} including: 
#' \itemize{
#'  \item{gamma}{A target estimate}
#'  \item{xi}{An auxiliary estimate in the ADMM algorithm}
#'  \item{phi}{An auxiliary estimate in the ADMM algorithm}
#'  \item{BIC}{A sequence of BIC values w.r.t. different lambdas}
#'  \item{lambda}{A sequence of tuning parameters used in the algorithm.}
#'  \item{omega}{A $p x 1$ weight matrix used in the algorithm.}
#' }
#' 
qrglasso <- function(Y, W, L, omega = NULL, tau = 0.5, qn = 1, lambda = NULL, maxit = 1000, thr = 1e-04){
  if(is.null(lambda)) {
    nlambda <- 51
    max.lambda <- 10
    lambda <- c(0, exp(seq(log(max.lambda / 1e4), log(max.lambda), length = (nlambda - 1))))
  } else {
    nlambda <- length(lambda)
  }
  
  zeta <- 10
  zetaincre <- 1
  
  if(is.null(omega))
    result <- awgl(Y, W, lambda, tau, L, qn, zeta, zetaincre, maxit, thr)
  else
    result <- awgl_omega(Y, W, omega, lambda, tau, qn, zeta, zetaincre, maxit, thr)
  
  result$phi[, 1] <- result$gamma[, 1]
  
  obj.cv <- list(gamma = result$gamma,
                 xi = result$xi,
                 phi = result$phi,
                 BIC = result$BIC,
                 lambda = lambda,
                 omega = result$omega)
  
  class(obj.cv) <- "qrglasso"
  return(obj.cv)
}
