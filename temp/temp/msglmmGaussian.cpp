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
mat SigmaSampler(vec a, mat Psi){
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



// Density of the multivariate SSGLMM for Gaussian data
double dGaussMulti(mat Y, List X, List M, mat beta, mat Psi, vec S){

  int J = Y.n_cols;
  int n = Y.n_rows;
  
  vec eta(n);
  vec ll(J);
  
  for(int j = 0; j < J; j++){
    double s2 = S(j);
    mat Xj    = X(j);
    mat Mj    = M(j);
    vec Yj    = Y.col(j);
    eta       = Yj-Xj*beta.col(j)-Mj*Psi.col(j);     
    ll(j) = - sum(0.5 * n * log(s2) + 0.5 * trans(eta)*eta / s2);
  }
  
  return(sum(ll));
}


// Sample the fixed effects for a single outcome
vec betaSgauss(vec y, mat X, double taub, double s2){

  int p = X.n_cols;

  mat Ip(p,p, fill::eye);
  
  mat Xty = trans(X) * y;
  mat XtX = trans(X) * X;

  mat Covb = XtX + Ip * s2/taub;
    Covb   = inv(Covb);
  
  vec bnorm = Rcpp::rnorm(p,0,1);
  
  vec beta = Covb * Xty + trans(chol(Covb))*bnorm;
    
  return(beta);
}


// Sample the spatial effects
vec psiSgauss(mat Psi, int num, double sigmas, vec Y, mat X, mat M, 
              mat SigmaI){

  int q = Psi.n_rows;

  vec l1(1); vec l2(1); vec l3(1);
  double uP, ll;
  uP = R::runif(0.0, 1.0);

  vec psi     = Psi.col(num);
  vec psiProp = Rcpp::rnorm(q, 0.0, sigmas);
  vec psiP    = psiProp + psi;

  mat PsiP      = Psi;
  PsiP.col(num) = psiP;
  
  vec Psiv  = vectorise(Psi);
  vec PsiPv = vectorise(PsiP);
  
  vec eta     = X*beta + M*psi;
  vec lambda  = exp(eta);
  
  vec etaP    = X*beta + M*psiP;
  vec lambdaP = exp(etaP);
  
  l1 = trans(Y) * ( log(lambdaP) - log(lambda)   );
  l2 = sum(lambda - lambdaP);
  l3 = 0.5 * ( trans(Psiv) * SigmaI * Psiv - trans(PsiPv) * SigmaI * PsiPv );
  ll = sum(l1+l2+l3);
  
  if(log(uP) < ll) psi = psiP;
    
  return(psi);
}


// [[Rcpp::export]]
SEXP msglmmGaussian(int iterations, mat Y, List X, List M, vec sigmab,
                    vec sigmas, double taub, mat beta0, mat Psi0, mat Sigma0,
                    int verbose){

  // Setup constants
  int J = Y.n_cols;
  int q = Psi0.n_rows;
  
      
  //Initialize the chains
  List PsiChain(iterations);
  List bChain(iterations);
  List SigmaChain(iterations);
  mat  aChain(J,iterations,fill::zeros);
  
  PsiChain(0)     = Psi0;
  bChain(0)       = beta0; 
  SigmaChain(0)   = Sigma0; 
  aChain.col(0)   = ones<vec>(J);
    aChain.col(0) = aChain.col(0)*4 / J;
 
  //Initialize covariance matrix
  mat SigInv(J,J);  
  mat SigI(J*q, J*q);
  mat Idq(q,q, fill::eye);

  // Initialize DIC calculator
  mat Beta     = bChain(0);
  mat Psi      = PsiChain(0);
  vec v(iterations);
  v(0) = -2 * dPoissMulti(Y, X, M, Beta, Psi);
  
  if(verbose=1){
    Rprintf("Sampling: MCMC can consume a lot of time.");
    #ifdef Win32
      R_FlushConsole();
    #endif
  }

  Rcpp::checkUserInterrupt();  
 
  //MCMC iterations
  for(int k=1; k < iterations; k++){
  
    mat Sig = SigmaChain(k-1);
    SigInv  = inv(Sig);
    SigI    = kron(SigInv, Idq);

    mat Beta = bChain(k-1);
    mat Psi  = PsiChain(k-1);
    
    Rcpp::checkUserInterrupt();
    
    for(int j=0; j<J; j++){
      Beta.col(j)  = betaSpois(Beta.col(j), sigmab(j), taub, Y.col(j), X(j), M(j), Psi.col(j)); 
      Psi.col(j)   = psiSpois(Psi, j, sigmas(j), Y.col(j), X(j), M(j), Beta.col(j), SigI);
      aChain(j,k)  = aSampler(SigInv(j,j),J);
    }
    
    SigmaChain(k)  = SigmaSampler(aChain.col(k), Psi);
    v(k)           = -2 * dPoissMulti(Y, X, M, Beta, Psi);

    bChain(k)     = Beta;
    PsiChain(k)   = Psi;  
  } 
  
  return Rcpp::List::create(Rcpp::Named("bChain")     = bChain,
                            Rcpp::Named("PsiChain")   = PsiChain,
                            Rcpp::Named("SigmaChain") = SigmaChain,
                            Rcpp::Named("aChain")     = aChain,
                            Rcpp::Named("v")          = v);
}

