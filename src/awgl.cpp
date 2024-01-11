// includes from the plugin
// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
using namespace arma;
using namespace Rcpp;
using namespace std;


arma::vec soft1(const arma::colvec v, const double tau) {
  int n = v.n_elem;
  arma::vec ans;
  double b = tau / n;
  double c = (1 - tau) / n;
  ans.zeros(n);
  
  for(int  j = 0; j < n; ++j) {
    if(v[j] > b)
      ans[j] = v[j] - b;
    else if(v[j] < - c)
      ans[j] = v[j] + c;
    else
      ans[j] = 0;
  }
  
  return(ans);
}

arma::vec soft3(const arma::colvec v,
                const arma::mat omega,
                const double lambda,
                const int L) {
  int m, pL = v.n_elem;
  arma::vec ans, zero;
  double temp;
  int p = (int) pL / L;
  
  ans.zeros(pL);
  if(L > 1)
    zero.zeros(L - 1);
  
  for(int j = 0; j < p; j++) {
    m = j * L;
    //constant
    temp = 1 - lambda * omega(j, 0) / norm(v.row(m), 2);
    if(temp > 0)
      ans.row(m) =  temp * v.row(m);
    else
      ans.row(m) = 0;
    //nonconstant
    if(L > 1) {
      temp = 1 - lambda * omega(j, 1) / norm(v.rows(m + 1, m + L - 1), 2);
      if(temp > 0)
        ans.rows(m + 1, m + L - 1) = temp * v.rows(m + 1, m + L - 1);
      else
        ans.rows(m + 1, m + L - 1) = zero;
    }
  }
  
  return(ans);
}

arma::vec scad_derivative(const arma::colvec v, 
                          const double lambda,
                          const double a){
  int pL = v.n_elem;
  arma::vec absv = abs(v);
  arma::vec ans;
  ans.zeros(pL);
  
  for(int i = 0; i < pL; i++) {
    if(absv[i] < a * lambda && absv[i] > 0) {
      if(absv[i] <= lambda)
        ans[i] = lambda;
      else
        ans[i] = (a * lambda - absv[i]) / (a - 1);
      ans[i] *= (v[i] > 0 ? 1:-1);
    }
    else if(v[i] == 0) {
      ans[i] = lambda;
    }
  }
  
  return(ans);
}

arma::mat omega_weight(const arma::colvec v, const int p, const int L) {
  int m;
  arma::mat ans;
  ans.zeros(p, 2);
  
  for(int j = 0; j < p; j++) {
    m = j * L;
    ans(j, 0) = std::abs(v[m]) / sqrt(L);
    ans(j, 1) = norm(v.rows(m + 1 , m + L - 1), 2) / sqrt(L);
  }
  
  return(ans);
}

arma::mat BIC(const arma::mat xi,
              const arma::colvec gamma,
              const double tau,
              const int L,
              const int qn) {
  int n = xi.n_rows, pL = gamma.n_elem, m;
  int Sc = 0, Sv = 0;
  double loss = 0;
  int p = (int) pL / L;
  arma::mat ans;
  ans.zeros(1, 2);
  
  for(int j = 0; j < p; j++){
    m = j * L;
    //constant
    if(gamma[m] != 0) 
      Sc++;
    //nonconstant
    for(int k = 0; k < L - 1; k++)
      if(gamma[m + 1 + k] != 0)
        Sv++;
  }
  
  for(int i = 0; i < n; i++){
    if(xi[i] > 0)
      loss += tau * xi[i];
    else
      loss -= (1 - tau) * xi[i];
  }
  
  ans[0] = loss / n + (Sc + (L - 1) * Sv) * (qn * log(n) / (2.0 * n));
  ans[1] = log(loss / n) + (Sc + (L - 1) * Sv) * (qn * log(n) / (2.0 * n));
  return(ans);
}

void qrcore(const arma::mat Y,
            const arma::mat W,
            const arma::mat Wt,
            const arma::mat Winv,
            const arma::mat Winvt,
            const arma::mat omega,
            arma::mat &gamma,
            arma::mat &xi,
            arma::mat &phi,
            arma::mat &theta1,
            arma::mat &theta2,
            const double lambda,
            const double tau,
            double zeta,
            const double zetaincre,
            const int maxit,
            const double tol) {
  int n = Y.n_rows;
  int p = omega.n_rows;
  int pL = W.n_cols;
  int L = (int) pL / p;
  int iter = 0;
  arma::vec er(3);
  
  arma::mat gammaold = gamma;
  arma::mat xiold = xi;
  arma::mat phiold = phi;
  arma::mat theta1old = theta1;
  arma::mat theta2old = theta2;
  arma::mat temp1, temp2;
  
  for(iter = 0; iter < maxit; iter++) {
    gamma = Winvt * (Y - xiold - theta1old / zeta) + Winv * (phiold -theta2old / zeta);
    xi = soft1(zeta * Y - zeta * W * gamma - theta1old, tau) / zeta;
    phi= soft3(gamma + theta2old / zeta, omega, lambda/zeta, L);
    temp1 = xi - Y + W * gamma;
    temp2 = gamma - phi;
    theta1 = theta1old + zeta*(temp1);
    theta2 = theta2old + zeta*(temp2);
    
    er[0] = arma::norm(temp1, "fro") / sqrt(n / 1.0);
    er[1] = arma::norm(temp2, "fro") / sqrt(pL / 1.0);
    er[2] = arma::norm(gamma - gammaold, "fro") / sqrt(pL / 1.0);
    
    if(max(er) <= tol)
      break;
    gammaold = gamma;
    xiold = xi;
    phiold = phi;
    theta1old = theta1;
    theta2old = theta2;
    zeta *= zetaincre;
  }
  
  iter++;
  if(iter == maxit)
    Rcpp::Rcout << "Not converge at lambda=" << lambda << "\n" << std::endl;
}

void qrinit(const arma::mat Y,
            const arma::mat W,
            arma::mat &Wt,
            arma::mat &gamma,
            arma::mat &xi,
            arma::mat &theta1,
            const double tau,
            double zeta,
            double zetaincre,
            int maxit,
            double tol){
  int n = Y.n_rows;
  int pL = W.n_cols;
  int iter = 0;
  arma::mat Winv_initial, IpL, gammaold, xiold;
  arma::vec er(2);
  Wt = arma::trans(W);
  IpL.eye(pL, pL);
  if(pL < n)
    Winv_initial = arma::inv_sympd(Wt * W) * Wt;
  else
    Winv_initial = arma::inv_sympd(Wt * W + log(pL) / (n * 1.0) * IpL) * Wt;
  
  gammaold = Winv_initial * Y;
  xiold.zeros(n, 1);
  arma::mat theta1old = theta1;
  arma::mat temp1;
  
  for(iter = 0; iter < maxit; iter++) {
    gamma = Winv_initial * (Y - xiold - theta1old/zeta);
    xi = soft1(zeta * Y - zeta * W * gamma - theta1old, tau) / zeta;
    temp1 = xi - Y + W * gamma; 
    theta1 = theta1old + zeta*(temp1);
    
    er[0] = arma::norm(temp1, "fro")/sqrt(n / 1.0);
    er[1] = arma::norm(gamma - gammaold, "fro") / sqrt(pL / 1.0);
    
    if(max(er) <= tol)
      break;
    
    gammaold = gamma;
    xiold = xi;
    theta1old = theta1;
    zeta *= zetaincre;
  }
  
  iter++;
  if(iter == maxit)
    Rcpp::Rcout << "Not converge with error" << max(er) << "\n" << std::endl;
}

//' Internal function: Quantile regression with adaptively group lasso with the input omega
//' @keywords internal
//' 
//' @param Y data matrix (n x 1)
//' @param W B-splines with covariates matrix (n x pL)
//' @param omega Weights for group lasso
//' @param lambda A sequence of tuning parameters 
//' @param tau A quantile of interest
//' @param qn A bound parameter for HDIC
//' @param zeta A step parameter
//' @param zetaincre An increment of each step
//' @param maxit The maximum number of iterations
//' @param tol A tolerance rate 
//' @return A list of selected parameters
// [[Rcpp::export]]
Rcpp::List awgl_omega(const arma::mat Y,
                      const arma::mat W,
                      const arma::mat omega,
                      const arma::vec lambda,
                      const double tau,
                      const int qn,
                      double zeta,
                      double zetaincre,
                      int maxit,
                      double tol) {
  /* Quantile regression with adaptively group lasso with the input omega
   * Returns
   *    - gamma: target estimate
   *    - xi, phi: auxiliary estimate in the ADMM algorithm
   *    - theta1, theta2: Lagrangian multipliers
   *    - BIC: BIC values of different lambdas
   */
  int pL = W.n_cols;
  int n_lambda = lambda.n_elem;
  int n = Y.n_rows;
  int p = omega.n_rows;
  int L = (int) pL / p;
  
  arma::mat gamma;
  arma::mat xi;
  arma::mat phi;
  arma::mat theta1;
  arma::mat theta2;
  arma::mat Winv, Wt, Winvt, gammaold, xiold, phiold, theta1old, theta2old, IpL, BIC_lambda;
  gamma.zeros(pL, n_lambda);
  xi.zeros(n, n_lambda);
  phi.zeros(pL, n_lambda);
  theta1.zeros(n , n_lambda);
  theta2.zeros(pL, n_lambda);
  IpL.eye(pL, pL);
  BIC_lambda.zeros(n_lambda, 2);
  
  /* compute the gamma, xi, theta1 w.r.t lambda = 0, and compute the correponding BIC*/
  gammaold = gamma.col(0);
  xiold = xi.col(0);
  phiold = phi.col(0);
  theta1old = theta1.col(0);
  theta2old = theta2.col(0);
  qrinit(Y, W, Wt,gammaold, xiold, theta1old, tau, zeta, zetaincre, maxit, tol);
  BIC_lambda.row(0) = BIC(xiold, gammaold, tau, L, qn);
  
  /*compute the gamma, xi, phi, theta1, theta2 w.r.t lambda, and compute the correponding BIC*/
  gamma.col(0) = gammaold;
  xi.col(0) = xiold;
  theta1.col(0) = theta1old;
  theta2.col(0) = theta2old;
  
  if(max(lambda) != 0) {
    if(pL < n)
      Winv = arma::inv_sympd(Wt * W + IpL);
    else
      Winv = arma::inv_sympd(Wt * W + (1 + (log(pL) / (n + 0.0))) * IpL);
    
    Winvt = Winv * Wt;
    for(int i = 1; i < n_lambda; i++) {
      gammaold = gamma.col(i - 1);
      xiold = xi.col(i - 1);
      phiold = phi.col(i - 1);
      theta1old = theta1.col(i - 1);
      theta2old = theta2.col(i - 1);
      qrcore(Y, W, Wt, Winv, Winvt, omega, gammaold, xiold, phiold, theta1old, theta2old, lambda[i], tau, zeta, zetaincre, maxit, tol);
      BIC_lambda.row(i) = BIC(xiold, phiold, tau, L, qn);
      gamma.col(i) = gammaold;
      xi.col(i) = xiold;
      phi.col(i) = phiold;
      theta1.col(i) = theta1old;
      theta2.col(i) = theta2old;
      
      /* If all elements of phi are zeros at some lambda, break out the loop.*/
      if(sum(abs(phi.col(i))) == 0) {
        for(int j = i + 1; j < n_lambda; j++)
          BIC_lambda.row(j) = BIC_lambda.row(i);
        Rcpp::Rcout << "All values of gamma are zeros when lambda >" << lambda[i] << "\n" << std::endl;
        break;
      }  
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("gamma") = gamma,
                            Rcpp::Named("xi") = xi, 
                            Rcpp::Named("phi") = phi,
                            Rcpp::Named("theta1") = theta1, 
                            Rcpp::Named("theta2") = theta2,
                            Rcpp::Named("BIC") = BIC_lambda);
  
}

//' Internal function: Quantile regression with adaptively group lasso without input Omega
//' 
//' @param Y data matrix (n x 1)
//' @param W B-splines with covariates matrix (n x pL)
//' @param lambda A sequence of tuning parameters 
//' @param tau A quantile of interest
//' @param L The number of groups
//' @param qn A bound parameter for HDIC
//' @param zeta A step parameter
//' @param zetaincre An increment of each step
//' @param maxit The maximum number of iterations
//' @param tol A tolerance rate 
//' @return A list of selected parameters
// [[Rcpp::export]]
Rcpp::List awgl(const arma::mat Y,
                const arma::mat W,
                const arma::vec lambda,
                const double tau,
                const int L,
                const int qn,
                double zeta,
                double zetaincre,
                int maxit,
                double tol) {
  /* Quantile regression with adaptively group lasso without input Omega
   * Returns
   *    - gamma: target estimate
   *    - xi, phi: auxiliary estimate in the ADMM algorithm
   *    - theta1, theta2: Lagrangian multipliers
   *    - BIC: BIC values of different lambdas
   *    - omega: estimate weights for group lasso
   */
  
  int pL = W.n_cols;
  int n_lambda = lambda.n_elem;
  int n = Y.n_rows;
  int p = (int) pL / L;
  float scad_weight = 3.7;
  
  arma::mat gamma, xi, phi, theta1, theta2, gammaold, xiold, phiold, theta1old, theta2old;
  arma::mat Winv, Wt, Winvt, omega_fake, IpL, BIC_lambda, omega;
  gamma.zeros(pL, n_lambda);
  xi.zeros(n, n_lambda);
  phi.zeros(pL, n_lambda);
  theta1.zeros(n , n_lambda);
  theta2.zeros(pL, n_lambda);
  IpL.eye(pL, pL);
  BIC_lambda.zeros(n_lambda, 2);
  omega_fake.ones(pL, 1);
  omega.zeros(p, 1);
  /* compute the gamma, xi, theta1 w.r.t lambda = 0, and compute the correponding BIC*/
  gammaold = gamma.col(0);
  xiold = xi.col(0);
  phiold = phi.col(0);
  theta1old = theta1.col(0);
  theta2old = theta2.col(0);
  qrinit(Y, W, Wt, gammaold, xiold, theta1old, tau, zeta, zetaincre, maxit, tol);
  BIC_lambda.row(0) = BIC(xiold, gammaold, tau, L, qn);
  
  /*compute the gamma, xi, phi, theta1, theta2 w.r.t lambda, and compute the correponding BIC*/
  gamma.col(0) = gammaold;
  xi.col(0) = xiold;
  theta1.col(0) = theta1old;
  theta2.col(0) = theta2old;
  
  if(max(lambda) != 0) {
    if(pL < n)
      Winv = arma::inv_sympd(Wt * W + IpL);
    else
      Winv = arma::inv_sympd(Wt * W + (1 + (log(pL) / (n + 0.0))) * IpL);
    
    Winvt = Winv * Wt;
    
    /* compute omega */
    for(int i = 1; i < n_lambda; i++) {
      gammaold = gamma.col(i - 1);
      xiold = xi.col(i - 1);
      phiold = phi.col(i - 1);
      theta1old = theta1.col(i - 1);
      theta2old = theta2.col(i - 1);
      qrcore(Y, W, Wt, Winv, Winvt, omega_fake, gammaold, xiold, phiold, theta1old, theta2old, lambda[i], tau, zeta, zetaincre, maxit, tol);
      BIC_lambda.row(i) = BIC(xiold, phiold, tau, L, qn);
      gamma.col(i) = gammaold;
      xi.col(i) = xiold;
      phi.col(i) = phiold;
      theta1.col(i) = theta1old;
      theta2.col(i) = theta2old;
      
      /* If all elements of phi are zeros at some lambda, break out the loop.*/
      if(sum(abs(phi.col(i))) == 0) {
        for(int j = i + 1; j < n_lambda; j++)
          BIC_lambda.row(j) = BIC_lambda.row(i);
        break;
      }
    }
    uword index1;
    // If .col(1), refer to BIC with log term
    (BIC_lambda.col(0)).min(index1);
    arma::mat weight_scad_deriv = scad_derivative(abs(gamma.col(index1)), lambda[index1], scad_weight);
    omega = omega_weight(weight_scad_deriv, p, L);
    
    /* main procedure */
    for(int i = 1; i < n_lambda; i++) {
      gammaold = gamma.col(i - 1);
      xiold = xi.col(i - 1);
      phiold = phi.col(i - 1);
      theta1old = theta1.col(i - 1);
      theta2old = theta2.col(i - 1);
      qrcore(Y, W, Wt, Winv, Winvt, omega, gammaold, xiold, phiold, theta1old, theta2old, lambda[i], tau, zeta, zetaincre, maxit, tol);
      BIC_lambda.row(i) = BIC(xiold, phiold, tau, L, qn);
      gamma.col(i) = gammaold;
      xi.col(i) = xiold;
      phi.col(i) = phiold;
      theta1.col(i) = theta1old;
      theta2.col(i) = theta2old;
      
      /* If all elements of phi are zeros at some lambda, break out the loop.*/
      if(sum(abs(phi.col(i))) == 0) {
        for(int j = i+1; j < n_lambda; j++)
          BIC_lambda.row(j) = BIC_lambda.row(i);
        
        Rcpp::Rcout << "All values of gamma are zeros when lambda >" << lambda[i] << "\n" << std::endl;
        break;
      }  
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("gamma") = gamma,
                            Rcpp::Named("xi") = xi, 
                            Rcpp::Named("phi") = phi,
                            Rcpp::Named("theta1") = theta1, 
                            Rcpp::Named("theta2") = theta2,
                            Rcpp::Named("BIC") = BIC_lambda,
                            Rcpp::Named("omega") = omega);
  
}