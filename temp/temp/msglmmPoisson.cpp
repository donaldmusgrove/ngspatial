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



// Density (log) of the multivariate SSGLMM for Poisson data
double dPoissMulti(mat Y, List X, List M, mat offset, mat beta, mat Psi){

  int J = Y.n_cols;
  int n = Y.n_rows;
  
  vec lambda(n);
  vec eta(n);
  vec ll(J);
  
  for(int j = 0; j < J; j++){
    mat Xj = X(j);
    mat Mj = M(j);
    eta    = offset.col(j) + Xj*beta.col(j) + Mj*Psi.col(j);
    lambda = exp(eta);
    vec Yj = Y.col(j);
    ll(j) = sum( Yj % log(lambda) - lambda);
  }
  
  return(sum(ll));
}


// Sample the fixed effects for a single outcome
vec betaSpois(vec beta, double sigmab, double taub, vec Y, mat X, mat M, vec offset, vec psi){

  int J  = beta.n_elem;
  vec l1(1); vec l2(1); vec l3(1);
  double uP, ll;
  
  uP           = R::runif(0.0, 1.0);
  vec betaProp = Rcpp::rnorm(J, 0.0, sigmab);
  vec betaP    = betaProp + beta;

  vec eta    = offset + X*beta + M*psi;
  vec lambda = exp(eta);
  
  vec etaP    = offset + X*betaP + M*psi;
  vec lambdaP = exp(etaP);
  
  l1 = trans(Y) * ( log(lambdaP) - log(lambda) );
  l2 = sum(lambda - lambdaP);
  l3 = 0.5 * taub * ( dot(beta, beta) - dot(betaP, betaP) );
  ll = sum(l1+l2+l3);
  
  if(log(uP) < ll) beta = betaP;
  
  return(beta);
}


// Sample the spatial effects for a single outcome
vec psiSpois(mat Psi, int num, double sigmas, vec Y, mat X, mat M, vec offset, 
           vec beta, mat SigmaI){

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
  
  vec eta     = offset + X*beta + M*psi;
  vec lambda  = exp(eta);
  
  vec etaP    = offset + X*beta + M*psiP;
  vec lambdaP = exp(etaP);
  
  l1 = trans(Y) * ( log(lambdaP) - log(lambda)   );
  l2 = sum(lambda - lambdaP);
  l3 = 0.5 * ( trans(Psiv) * SigmaI * Psiv - trans(PsiPv) * SigmaI * PsiPv );
  ll = sum(l1+l2+l3);
  
  if(log(uP) < ll) psi = psiP;
    
  return(psi);
}


// [[Rcpp::export]]
SEXP msglmmPoisson(int iterations, mat Y, List X, List M, mat offset, vec sigmab,
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
  v(0) = -2 * dPoissMulti(Y, X, M, offset, Beta, Psi);
  
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
      Beta.col(j)  = betaSpois(Beta.col(j), sigmab(j), taub, Y.col(j), X(j), M(j), offset.col(j), Psi.col(j)); 
      Psi.col(j)   = psiSpois(Psi, j, sigmas(j), Y.col(j), X(j), M(j), offset.col(j), Beta.col(j), SigI);
      aChain(j,k)  = aSampler(SigInv(j,j),J);
    }
    
    SigmaChain(k)  = SigmaSampler(aChain.col(k), Psi);
    v(k)           = -2 * dPoissMulti(Y, X, M, offset, Beta, Psi);

    bChain(k)     = Beta;
    PsiChain(k) = Psi;  
  } 
  
  return Rcpp::List::create(Rcpp::Named("bChain")     = bChain,
                            Rcpp::Named("PsiChain")   = PsiChain,
                            Rcpp::Named("SigmaChain") = SigmaChain,
                            Rcpp::Named("aChain")     = aChain,
                            Rcpp::Named("v")          = v);
}


