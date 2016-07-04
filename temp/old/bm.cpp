#include <RcppArmadillo.h>
#include <utils.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
// Create a block diagonal matrix
mat blockDiagT(List M) {

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



// [[Rcpp::export]]
// Batchmeans on a single vector
vec bmcpp(vec X){

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


// [[Rcpp::export]]
// Batchmeans on a matrix
mat bmmatcpp(mat M){

  int num = M.n_cols;
  mat bmvals(num,2);
  vec bmout(2);

  for (int i=0; i<num; i++){
    bmout         = bmcpp(M.col(i));
    bmvals.row(i) = trans(bmout);
  }

  return(bmvals);
}


// [[Rcpp::export]]
// Join two arrays/cubes
cube arraybind(cube M1, cube M2){
 
  cube M = join_slices(M1,M2);

  return(M);
}



// [[Rcpp::export]]
// Transform covariance matrix to correlation matrix
mat cov2corcpp(mat V){
  vec Isv = V.diag();
      Isv = 1/sqrt(Isv);
  mat Is = diagmat(Isv);
  mat R  = Is*V*Is; 
  return(R);
}


// [[Rcpp::export]]
// Apply cov2corcpp to each slice of an array and return an array
cube acov2corcpp(cube V){
  int  n = V.n_slices;
  cube R = V;
  
  for (int i = 0; i < n; i++){
    R.slice(i) = cov2corcpp(V.slice(i));
  }  

  return(R);
}






