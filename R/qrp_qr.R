#' Adaptively weighted group Lasso
#'
#' @param Y 
#' @param W
#' @param omega
#' @param tau
#' @param qn
#' @param lambda
#' @param maxit
#' @param thr
#' @export
#' @examples
#' 
#' @return This function returns a \code{list} including: 
#' \itemize{
#'  \item gamma
#'  \item xi
#'  \item phi
#'  \item BIC
#'  \item lambda
#'  \item omega
#' }
#' 
grp_qr <- function(Y, W, L, omega = NULL, tau, qn = 1, lambda = NULL, maxit = 1000, thr = 1e-04){
  
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
    result <- qr_rcpp(Y, W, lambda, tau, L, qn, zeta, zetaincre, maxit, thr)
  else
    result <- qr_rcpp_omega(Y, W, omega, lambda, tau, qn,zeta, zetaincre, maxit, thr)
  
  result$phi[, 1] <- result$gamma[, 1]
  
  obj.cv <- list(gamma = result$gamma,
                 xi = result$xi,
                 phi = result$phi,
                 BIC = result$BIC,
                 lambda = lambda,
                 omega = result$omega)
  
  return(obj.cv)
}

#' B-spline orthogonalization
#'
#' @param LLI
#' @param degree
#' @param bkn
#' @param crep_input
#' @param is_approx
#' 
#' @export
#' @examples
#' @return This function returns a \code{list} including: 
#' \itemize{
#'  \item b_function
#'  \item z
#' }
B_orth <- function(LLI, degree, bkn, crep_input = NULL, is_approx = FALSE){
  if(is_approx) {
    if(is.null(crep_input)) {
      nrep <- 1e6
      crep <- 1:nrep
      crep <- crep / nrep
    } else if(length(crep_input) < 1e6) {
      nrep <- 1e6
      crep <- sort(c(crep_input,
                     runif(1e6 - length(crep_input), min = 0, max = 1)))
    } else {
      nrep <- length(crep_input)
      crep <- crep_input
    }
  } else{
    if(is.null(crep_input)){
      nrep <- 1e4
      crep <- 1:nrep
      crep <- crep/nrep
    } else{
      nrep <- length(crep_input)
      crep <- crep_input
    }
  }
  
  bs0 <- bs(crep, knots = ikn,Boundary.knots = bkn, degree = degree, intercept = TRUE)
  LL <- dim(bs0)[2]
  bs1 <- matrix(0, nrow = nrep, ncol = LL)
  a0 <- diag(LL)
  
  sqLL<- 1 / sqrt(LL)
  bs1[, 1] <- sqLL
  bs1[, 2] <- sqrt(12) * sqLL * (crep - 0.5)
  
  for (ii in 3:LL)
    bs1[, ii] <- bs0[, (ii - 1)]

  for (jj in 3:LL){
    ee <- as.numeric(bs1[, jj])
    
    for (ii in 1:(jj - 1))
      a0[jj, ii] <- -LL*(mean(ee * bs1[, ii]))
    
    for (ii in 1:(jj - 1))
      ee <- ee + a0[jj, ii] * as.numeric(bs1[, ii])
    
    nee <- sqrt(mean(ee ** 2))
    bs1[,jj] <- ee / nee * sqLL
    a0[jj,] <- a0[jj, ] * sqLL / nee
  }
  
  if(is.null(crep_input) || is_approx == FALSE) {
    b_function <- bs1
    z <- crep
  } else {
    idx <- sapply(1:length(crep_input), function(i) which(crep == crep_input[i])[1])
    b_function <- bs1[unlist(idx), ]
    z <- crep_input
  }

  return(list(b_function = b_function, z = z))
}

