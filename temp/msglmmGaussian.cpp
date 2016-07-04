#include <RcppArmadillo.h>
#include <utils.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;



// Density of the multivariate SSGLMM for Gaussian data
// [[Rcpp::export]]
double dGaussMulti(mat Y, List X, List M, List beta, vec Psi, vec tau2e){

  int J = Y.n_cols;
  int n = Y.n_rows;
  int q = Psi.size()/J;
  vec ll(J);
  
  for(int j = 0; j < J; j++){
    mat Xj      = X(j);
    mat Mj      = M(j);
    vec Yj      = Y.col(j);
    vec betaj   = beta(j);
    vec psij    = vectorSubset(Psi, q*(j+1)-q, q*(j+1)-1);
    vec eta     = Yj - Xj*betaj - Mj*psij;
        eta     = dot(eta,eta);
    ll(j)       = sum(n * 0.5 * log(tau2e(j)) - 0.5 * tau2e(j) * eta);
  }
  
  return(sum(ll));
}


// Create a block diagonal matrix
// [[Rcpp::export]]
mat blockDiag(List M) {

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


// Sample the fixed effects for a single outcome
vec betaSgauss(double taub, vec Y, mat X, double tau2e){

  int J  = X.n_cols;
  mat Id = eye<mat>(J,J);

  mat bCov  = X.t() * X + taub / tau2e*Id;
      bCov  = inv_sympd( bCov ) ;
  vec bMean = bCov * X.t()*Y;
  
  mat bChol = chol(bCov, "lower");
  vec bNorm = Rcpp::rnorm(J,0,1);
  
  vec beta = bMean + sqrt(1/tau2e)*bChol*bNorm;
  
  return(beta);
}


// Sample the spatial effects for a single outcome
// [[Rcpp::export]]
mat psiSgauss(vec Y, mat Mdiag, vec tau2e, mat Sigma){

  int J       = Sigma.n_cols;
  int q       = Mdiag.n_cols / J;
  int n       = Mdiag.n_rows / J;
  mat Iq      = eye<mat>(q,q);
  mat In      = eye<mat>(n,n);
    
    
  mat S       = diagmat(1/tau2e);
      S       = kron(S, In);
  
  mat psiCov  = Mdiag.t() * S * Mdiag;
      psiCov  = psiCov + kron(inv(Sigma),Iq);
      psiCov  = inv(psiCov);
      
  mat psiMean = psiCov * Mdiag.t() * S * Y;
  
  mat psiChol = chol(psiCov, "lower");
  vec psiNorm = Rcpp::rnorm(q*J,0,1);
  vec Psiv    = psiMean + psiChol*psiNorm;
  
  return(Psiv);
}


// Sample the error variance for a single outcome
// [[Rcpp::export]]
double tau2eSgauss(vec Y, mat X, mat M, vec beta, vec psi, double as, double bs){

  int n       = X.n_rows;
  vec z       = Y - X*beta - M*psi;
  double zdot = dot(z,z);  
    
  double A = as + 0.5 * n;
  double B = 0.5 * zdot + 1/bs;
  
  return(R::rgamma(A,1/B));
}


// [[Rcpp::export]]
SEXP msglmmGaussian(int iterations, mat Y, List X, List M, 
                    double taub, double nu, double taus, double tau2ea, double tau2eb,
                    List beta0, vec Psi0, vec tau2e0, mat Sigma0,
                    int verbose, int runningIter, vec fivepct, int maxit){

  // Setup constants
  int J = Y.n_cols;
  int q = Psi0.size() / J;
  
  // Transform inputs
  vec Yvec  = vectorise(Y);
  mat Mdiag = blockDiag(M);
  mat Xdiag = blockDiag(X);
 
  
  //Initialize the chains
  mat  tau2eChain(J,iterations+1,fill::zeros);
  mat  PsiChain(q*J,iterations+1,fill::zeros);
  cube SigmaChain(J,J,iterations);
  mat  aChain(J,iterations+1,fill::zeros);

    
  //Fill the chains with initial values
  field<vec> bChain(iterations+1, J);
  List Betatemp(J);
  for(int j=0; j<J; j++){
    vec btemp   = beta0(j);
    bChain(0,j) = btemp;   
  }

  tau2eChain.col(0) = tau2e0;
  PsiChain.col(0)   = Psi0;
  SigmaChain        = join_slices(Sigma0,SigmaChain);
  aChain.col(0)     = ones<vec>(J);
    aChain.col(0)   = aChain.col(0)*4 / J;
 
  mat SigInv(J,J); 

   
  // Initialize DIC calculator
  vec v(iterations+1);
  v(0) = -2 * dGaussMulti(Y, X, M, beta0, PsiChain.col(0), tau2e0);
   
   
  Rcpp::checkUserInterrupt();
   
  //MCMC iterations
  for(int k=1; k < iterations+1; k++){
    Rcpp::checkUserInterrupt();
  
    if(verbose && isin(runningIter,fivepct)){
      Rprintf("Progress => %i", runningIter);
      Rprintf(" out of %i max iterations.\n", maxit);
      #ifdef Win32
      R_FlushConsole();
      #endif
    }
          
    mat Sig = SigmaChain.slice(k-1);
    SigInv  = inv(Sig);
    
    vec Psi         = psiSgauss(Yvec, Mdiag, tau2eChain.col(k-1), SigmaChain.slice(k-1));
    PsiChain.col(k) = Psi; 
    
    for(int j=0; j<J; j++){
      vec Psij        = vectorSubset(Psi, q*(j+1)-q, q*(j+1)-1);
      
      bChain(k,j)     = betaSgauss(taub, Y.col(j), X(j), tau2eChain(j,k-1));
      Betatemp(j)     = bChain(k,j); 
      tau2eChain(j,k) = tau2eSgauss(Y.col(j), X(j), M(j), bChain(k,j), Psij, tau2ea, tau2eb);
      aChain(j,k)     = aSampler(SigInv(j,j),J,nu,taus);
    }
       
    SigmaChain.slice(k) = SigmaSamplerVec(aChain.col(k), Psi, nu);
    v(k)                = -2 * dGaussMulti(Y, X, M, Betatemp, Psi, tau2eChain.col(k));
    
    runningIter = runningIter+1;
  } 
  
  
  return Rcpp::List::create(Rcpp::Named("bChain")      = bChain,
                            Rcpp::Named("PsiChain")    = PsiChain,
                            Rcpp::Named("SigmaChain")  = SigmaChain,
                            Rcpp::Named("aChain")      = aChain,
                            Rcpp::Named("tau2eChain")  = tau2eChain,
                            Rcpp::Named("v")           = v,
                            Rcpp::Named("runningIter") = runningIter);
}