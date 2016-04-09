#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Draw a sample from a Wishart distribution
mat rwishart(int df, mat S){
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
mat riwishart(int df, mat S){
  return rwishart(df,S.i()).i();
}


// Draw a sample from the posterior of the within-areal unit covariance matrix
mat SigmaSampler(vec a, mat Gamma){
  int J     = Gamma.n_rows;
  int q     = Gamma.n_cols;
  double nu = 2;
  
  mat PP    = trans(Gamma) * Gamma;
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
