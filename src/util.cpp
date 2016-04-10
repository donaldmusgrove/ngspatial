#include "util.h"
#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace std;
using namespace arma;


// Draw a sample from a Wishart distribution
mat rwishart(int df, const mat& S){
  int m = S.n_rows;
  mat Z(m,m);
  
  for(int i = 0; i < m; i++){
    Z(i,i) = sqrt(R::rchisq(df-i));
  }
  
  for(int j = 0; j < m; j++){  
    for(int i = j+1; i < m; i++){    
      Z(i,j) = R::rnorm(0,1);
    }
  }
  
  mat C = trimatl(Z).t() * chol(S);
  
  return C.t()*C;
}

// Draw a sample from an Inverse-Wishart distribution
mat riwishart(int df, const mat & S){
  return rwishart(df,S.i()).i();
}

// Draw a sample from the posterior of the within-areal unit covariance matrix
mat SigmaSampler(const vec& a, const mat& Psi){
  int J     = Psi.n_rows;
  int q     = Psi.n_cols;
  double nu = 2;
  
  mat PP    = trans(Psi) * Psi;
  mat Adiag = diagmat(a);
  
  mat SigmaHat = PP + 2 * nu * Adiag;
  
  return riwishart(nu+q+J-1, SigmaHat);
}

// Draw a sample from the posterior of the diagonal of the within-areal unit covariance matrix
double aSampler(double Sinv, int J){
  double A  = pow(10,-10);
  double nu = 2;

  double aShp = (nu+J)/2;
  double aRte = nu * Sinv + A;
  
  return(1/R::rgamma(aShp, 1/aRte));
}

//Carry out batchmeans on a single vector
vec bmcpp(const vec& X){

  int n  = X.n_elem;
  int b  = floor(pow(n,0.5));
  int a  = floor(n/b);
  double xl, xu, xmean, muhat, varhat, e1;

  vec y(a);
  vec v(a);
  vec est(2);

  for (int i=1; i<(a+1); i++){
    xl      =  (i - 1) * b ;
    xu      =  i * b - 1;
    xmean   =  mean(X.rows(xl,xu));     
    y(i-1)  =  xmean;
  }

  muhat = mean(y);

  v.fill(muhat);
  v  = y - v;
  e1 = sum(trans(v) * v);

  varhat = b * e1 / (a - 1);

  est(0) = muhat;
  est(1) = pow( varhat/n, 0.5 );

  return(est) ;
}


//Carry out batchmeans on the columns of a matrix
mat bmmatcpp(const mat & M){

  int num = M.n_cols;
  mat bmvals(num,2);
  vec bmout(2);

  for (int i=0; i<num; i++){
    bmout         = bmcpp(M.col(i));
    bmvals.row(i) = trans(bmout);
  }

  return(bmvals);
}


// Join two arrays
cube arraybind(const cube & M1, const cube & M2){
 
  cube M = join_slices(M1,M2);

  return(M);
}


RCPP_MODULE(util){
  Rcpp::function("rwishart", &rwishart);
  Rcpp::function("riwishart", &riwishart);
  Rcpp::function("SigmaSampler", &SigmaSampler);
  Rcpp::function("aSampler", &aSampler);
  Rcpp::function("bmcpp", &bmcpp);
  Rcpp::function("bmmatcpp", &bmmatcpp);
  Rcpp::function("arraybind", &arraybind);
}
