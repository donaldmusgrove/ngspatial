#include "msglmmutil.h"
#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace std;
using namespace arma;


// Create a block diagonal matrix
mat blockDiag(const List& M) {

  int J = M.size();
  int dimen1 = 0 ;
  int dimen2 = 0 ;
  int idx1 = 0;
  int idx2 = 0;
  ivec dimvec1(J);
  ivec dimvec2(J);
  
  for(int i=0; i < J; i++) {
    mat Mmat   = M(i);
    dimvec1(i) = Mmat.n_rows ; 
    dimvec2(i) = Mmat.n_cols ; 
    dimen1    += dimvec1(i) ;
    dimen2    += dimvec2(i) ;
  }

  mat Mdiag(dimen1,dimen2,fill::zeros);
   
  for(int i=0; i<J; i++) {
    mat Mmat = M(i);
    Mdiag.submat(idx1, idx2, idx1 + dimvec1(i)-1, idx2 + dimvec2(i)-1 ) = Mmat;
    idx1 = idx1 + dimvec1(i);
    idx2 = idx2 + dimvec2(i);
  }

  return(Mdiag);
}


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
mat riwishart(int df, const mat& S){
  return rwishart(df,S.i()).i();
}

// Draw a sample from the posterior of the within-areal unit covariance matrix
mat SigmaSampler(const vec& a, const mat& Psi, double nu){
  int J     = Psi.n_rows;
  int q     = Psi.n_cols;
  
  mat PP    = trans(Psi) * Psi;
  mat Adiag = diagmat(a);
  
  mat SigmaHat = PP + 2 * nu * Adiag;
  
  return riwishart(nu+q+J-1, SigmaHat);
}

// Vectorized (Psi) version of SigmaSampler
mat SigmaSamplerVec(const vec& a, const vec& Psi, double nu){
  int J     = a.size();
  int q     = Psi.size();
      q     = q / J;
  
  mat Pmat  = mat(Psi);
  mat PP    = reshape(Pmat,q,J);
      PP    = trans(PP) * PP;
  mat Adiag = diagmat(a);
  
  mat SigmaHat = PP + 2 * nu * Adiag;
  
  return riwishart(nu+q+J-1, SigmaHat);
}

// Draw a sample from the posterior of the diagonal of the within-areal unit covariance matrix
double aSampler(double Sinv, int J, double nu, double taus){
  //double taus  = pow(10,-10);
  //double nu = 2;

  double aShp = (nu+J)/2;
  double aRte = nu * Sinv + taus;
  
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
mat bmmatcpp(const mat& M){

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
cube arraybind(const cube& M1, const cube& M2){
 
  cube M = join_slices(M1,M2);

  return(M);
}

// Check if an integer is in a vector
bool isin(int num, const vec& vals){
  int n = vals.n_elem;
  
  bool isequal = false;

  for (int i = 0; i < n; i++){   
  
    if(abs(vals(i) - num) < 0.001){
      isequal = true;
      break;
    }
  
  }
  return(isequal);
}

// Create a block diagonal matrix
mat blockDiagT(const List& M) {

  int J = M.size();
  int dimen1 = 0 ;
  int dimen2 = 0 ;
  int idx1 = 0;
  int idx2 = 0;
  ivec dimvec1(J);
  ivec dimvec2(J);
  
  for(int i=0; i < J; i++) {
    mat Mmat   = M(i);
    dimvec1(i) = Mmat.n_rows ; 
    dimvec2(i) = Mmat.n_cols ; 
    dimen1    += dimvec1(i) ;
    dimen2    += dimvec2(i) ;
  }

  mat Mdiag(dimen1,dimen2,fill::zeros);
   
  for(int i=0; i<J; i++) {
    mat Mmat = M(i);
    Mdiag.submat(idx1, idx2, idx1 + dimvec1(i)-1, idx2 + dimvec2(i)-1 ) = Mmat;
    idx1 = idx1 + dimvec1(i);
    idx2 = idx2 + dimvec2(i);
  }

  return(Mdiag.t());
}

// Transform covariance matrix to correlation matrix
mat cov2corcpp(const mat& V){
  vec Isv = V.diag();
      Isv = 1/sqrt(Isv);
  mat Is = diagmat(Isv);
  mat R  = Is*V*Is; 
  return(R);
}

// Apply cov2corcpp to each slice of an array and return an array
cube acov2corcpp(const cube& V){
  int  n = V.n_slices;
  cube R = V;
  
  for (int i = 0; i < n; i++){
    R.slice(i) = cov2corcpp(V.slice(i));
  }  

  return(R);
}

// Subset a vector between indices a and b
vec vectorSubset(const vec& x, int a, int b){
  vec subset = x.subvec(a, b);
  return(subset);
}


RCPP_MODULE(msglmmutil){
  Rcpp::function("rwishart", &rwishart);
  Rcpp::function("riwishart", &riwishart);
  Rcpp::function("SigmaSampler", &SigmaSampler);
  Rcpp::function("aSampler", &aSampler);
  Rcpp::function("bmcpp", &bmcpp);
  Rcpp::function("bmmatcpp", &bmmatcpp);
  Rcpp::function("arraybind", &arraybind);
  Rcpp::function("isin", &isin);
  Rcpp::function("blockDiagT", &blockDiagT);
  Rcpp::function("cov2corcpp", &cov2corcpp);
  Rcpp::function("acov2corcpp", &acov2corcpp);
  Rcpp::function("vectorSubset", &vectorSubset);
  Rcpp::function("SigmaSamplerVec", &SigmaSamplerVec);
  Rcpp::function("blockDiag", &blockDiag);
}
