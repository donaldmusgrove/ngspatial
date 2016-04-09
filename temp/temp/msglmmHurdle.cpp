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
// [[Rcpp::export]]
mat SigmaSampler(vec a, vec psi1, vec psi2){
  int q     = psi1.n_elem;;
  double nu = 2;
  
  mat PP    = join_rows(psi1,psi2);
      PP    = trans(PP) * PP;
  mat Adiag = diagmat(a);
  
  mat SigmaHat = PP + 2 * nu * Adiag;
  
  return riwishart(nu+q+2-1, SigmaHat);
}

// Draw a sample from the posterior of the diagonal of the within-areal unit covariance matrix
// [[Rcpp::export]]
double aSampler(double Sinv, int J){
  double A  = pow(10,-10);
  double nu = 2;

  double aShp = (nu+J)/2;
  double aRte = nu * Sinv + A;
  
  return(1/R::rgamma(aShp, 1/aRte));
}


// Hurdle model density 
// [[Rcpp::export]]
double dHurdle(vec y, vec z, List X, List M, vec offset, vec beta1, vec beta2,
               vec psi1, vec psi2){
   
  mat X1 = X(0); mat M1 = M(0);
  mat X2 = X(1); mat M2 = M(1);
  
  vec eta1 = X1*beta1 + M1*psi1;
  vec p    = 1 + exp(-eta1);
      p    = 1/p;

  vec eta2   = offset + X2*beta2 + M2*psi2;  
  vec lambda = exp(eta2);

  double l1 = sum(trans(1-z) * log(1 - p));
  double l2 = sum(trans(z) * ( log(p) - log(1-exp(-lambda)) + y % log(lambda) - lambda ));
  
  return(sum(l1+l2));
}


// ZIP model density 
// [[Rcpp::export]]
double dZIP(vec y, vec z, List X, List M, vec offset, vec beta1, vec beta2,
            vec psi1, vec psi2){
   
  mat X1 = X(0); mat M1 = M(0);
  mat X2 = X(1); mat M2 = M(1);
  
  vec eta1 = X1*beta1 + M1*psi1;
  vec p    = 1 + exp(-eta1);
      p    = 1/p;

  vec eta2   = offset + X2*beta2 + M2*psi2;  
  vec lambda = exp(eta2);

  //Reparameterize p for the ZIP model, i.e., p* = p(1-exp(lambda))
  p = p%(1-exp(lambda));
  
  double l1 = sum(trans(1-z) * log(1 - p));
  double l2 = sum(trans(z) * ( log(p) - log(1-exp(-lambda)) + y % log(lambda) - lambda ));
  
  return(sum(l1+l2));
}



// [[Rcpp::export]]
vec beta1Shurd(vec beta1, double taub, double tuneb1, vec z, mat X, mat M, vec psi1){

  int jj  = beta1.n_elem;
  vec l1(1); vec l2(1); vec l3(1); 
  double ll;
  
  double uP     = R::runif(0.0, 1.0);
  vec beta1Prop = Rcpp::rnorm(jj, 0.0, tuneb1);
  vec beta1P    = beta1Prop + beta1;

  vec eta  = X*beta1  + M*psi1;
  vec p    = 1 + exp(-eta);
      p    = 1/p;
  
  vec etaP = X*beta1P + M*psi1;
  vec pP   = 1 + exp(-etaP);
      pP   = 1/pP;
  
  l1 = trans(z)   * ( log(pP) - log(p) );
  l2 = trans(1-z) * ( log(1-pP) - log(1-p) );
  l3 = 0.5 * ( dot(beta1, beta1) - dot(beta1P, beta1P)  ) / taub;
  ll = sum(l1+l2+l3);
  
  if(log(uP) < ll) beta1 = beta1P;
  
  return(beta1);
}


// [[Rcpp::export]]
vec psi1Shurd(vec psi1, double tunep1, vec z, mat X, mat M, vec beta1, vec psi2,
              mat SigmaI){

  int q  = psi1.n_elem;
  vec l1(1); vec l2(1); vec l3(1); 
  double ll;
  
  double uP    = R::runif(0.0, 1.0);
  vec psi1Prop = Rcpp::rnorm(q, 0.0, tunep1);
  vec psi1P    = psi1Prop + psi1;
  
  mat PsiM  = join_cols(psi1,psi2);
  mat PsiPM = join_cols(psi1P,psi2);
  vec Psi   = vectorise(PsiM);
  vec PsiP  = vectorise(PsiPM);
  
  vec eta  = X*beta1 + M*psi1;
  vec p    = 1 + exp(-eta);
      p    = 1/p;
  
  vec etaP = X*beta1 + M*psi1P;
  vec pP   = 1 + exp(-etaP);
      pP   = 1/pP;
  
  l1 = trans(z)   * ( log(pP) - log(p) );
  l2 = trans(1-z) * ( log(1-pP) - log(1-p) );
  l3 = 0.5 * ( trans(Psi) * SigmaI * Psi - trans(PsiP) * SigmaI * PsiP);
  ll = sum(l1+l2+l3);
  
  if(log(uP) < ll) psi1 = psi1P;
  
  return(psi1);
}


// [[Rcpp::export]]
vec beta2Shurd(vec beta2, double taub, double tuneb2, vec y, vec z, mat X, mat M, vec offset, vec psi2){

  int jj  = beta2.n_elem;
  vec l1(1); vec l2(1); vec l3(1);  vec l4(1);
  double ll;
  
  double uP    = R::runif(0.0, 1.0);
  vec beta2Prop = Rcpp::rnorm(jj, 0.0, tuneb2);
  vec beta2P    = beta2Prop + beta2;

  vec eta    = offset + X*beta2  + M*psi2;
  vec lambda = exp(eta);
  
  vec etaP    = offset + X*beta2P  + M*psi2;
  vec lambdaP = exp(etaP);
  
  l1 = trans(z)   * ( log(1-exp(-lambda)) - log(1-exp(-lambdaP)) );
  l2 = trans(y%z) * ( log(lambdaP) - log(lambda)   );
  l3 = trans(z)   * (lambda - lambdaP);
  l4 = 0.5 * ( dot(beta2, beta2) - dot(beta2P, beta2P) ) / taub;
  ll = sum(l1+l2+l3+l4);
  
  if(log(uP) < ll) beta2 = beta2P;
  
  return(beta2);
}


// [[Rcpp::export]]
vec psi2Shurd(vec psi2, double tunep2, vec y, vec z, mat X, mat M, vec offset, vec beta2, 
          vec psi1, mat SigmaI){

  int q  = psi2.n_elem;
  vec l1(1); vec l2(1); vec l3(1); vec l4(1);
  double ll;
  
  double uP    = R::runif(0.0, 1.0);
  vec psi2Prop = Rcpp::rnorm(q, 0.0, tunep2);
  vec psi2P    = psi2Prop + psi2;

  mat PsiM  = join_cols(psi1,psi2);
  mat PsiPM = join_cols(psi1,psi2P);
  vec Psi   = vectorise(PsiM);
  vec PsiP  = vectorise(PsiPM); 
  
  vec eta    = offset + X*beta2  + M*psi2;
  vec lambda = exp(eta);
  
  vec etaP    = offset + X*beta2 + M*psi2P;
  vec lambdaP = exp(etaP);
  
  l1 = trans(z)   * ( log(1-exp(-lambda)) - log(1-exp(-lambdaP)) );
  l2 = trans(y%z) * ( log(lambdaP) - log(lambda)   );
  l3 = trans(z)   * (lambda - lambdaP);
  l4 = 0.5 * ( trans(Psi) * SigmaI * Psi - trans(PsiP) * SigmaI * PsiP);
  ll = sum(l1+l2+l3+l4);
  
  if(log(uP) < ll) psi2 = psi2P;
  
  return(psi2);
}


// [[Rcpp::export]]
List msglmmHurdle(int iterations, vec y, vec z, List X, List M, vec offset,
                  double tuneb1, double tuneb2, double tunep1, double tunep2,
			            double taub, vec beta1, vec beta2, vec psi1, vec psi2, mat Sigma,
                  int verbose){

  int p = beta1.n_elem;
  int q = psi1.n_elem;  
  
  // Initialize the chains
  mat psi1chain(q,iterations,fill::zeros);
  mat psi2chain(q,iterations,fill::zeros);
  mat beta1chain(p,iterations,fill::zeros);
  mat beta2chain(p,iterations,fill::zeros);
  List SigmaChain(iterations);
  mat  aChain(2,iterations,fill::zeros);
  
  psi1chain.col(0)  = psi1;
  psi2chain.col(0)  = psi2;
  beta1chain.col(0) = beta1;
  beta2chain.col(0) = beta2;
  SigmaChain(0)     = Sigma; 
  aChain.col(0)     = ones<vec>(2);
    aChain.col(0)   = aChain.col(0)*4 / 2;   
  
  //Initialize covariance matrix
  mat SigInv(2,2);  
  mat SigI(2*q, 2*q);
  mat Idq(q,q, fill::eye);
  
  // Initialize DIC calculator for Hurdle and ZIP models
  vec vHurdle(iterations);
  vec vZIP(iterations);
  vHurdle(0) = -2 * dHurdle(y, z, X, M, offset, beta1, beta2, psi1, psi2);
  vZIP(0)    = -2 * dZIP(y, z, X, M, offset, beta1, beta2, psi1, psi2);
  
  if(verbose=1){
    Rprintf("Sampling: MCMC can consume a lot of time.\n");
    #ifdef Win32
      R_FlushConsole();
    #endif
  }
  
  Rcpp::checkUserInterrupt();  
  
  //MCMC iterations
  for(int k=1; k < iterations; k++ ){
    mat Sig = SigmaChain(k-1);
    SigInv  = inv(Sig);
    SigI    = kron(SigInv, Idq);
      
    beta1chain.col(k) = beta1Shurd(beta1chain.col(k-1), taub, tuneb1, z, X(0), M(0), psi1chain.col(k-1));
    beta2chain.col(k) = beta2Shurd(beta2chain.col(k-1), taub, tuneb2, y, z, X(1), M(1), offset, psi2chain.col(k-1));
    psi1chain.col(k)  = psi1Shurd(psi1chain.col(k-1), tunep1, z, X(0), M(0), beta1chain.col(k), psi2chain.col(k-1), SigI);
    psi2chain.col(k)  = psi2Shurd(psi2chain.col(k-1), tunep2, y, z, X(1), M(1), offset, beta2chain.col(k),  psi1chain.col(k), SigI);
      
    aChain(0,k)   = aSampler(SigInv(0,0),2);
    aChain(1,k)   = aSampler(SigInv(1,1),2);
        
    SigmaChain(k) = SigmaSampler(aChain.col(k), psi1chain.col(k), psi2chain.col(k));
   
    vHurdle(k) = -2 * dHurdle(y, z, X, M, offset, beta1chain.col(k), beta2chain.col(k), psi1chain.col(k), psi2chain.col(k));    
    vZIP(k)    = -2 * dZIP(y, z, X, M, offset, beta1chain.col(k), beta2chain.col(k), psi1chain.col(k), psi2chain.col(k));
  }
    
  return Rcpp::List::create(Rcpp::Named("beta1chain") = beta1chain,
                            Rcpp::Named("beta2chain") = beta1chain,
                            Rcpp::Named("psi1chain")  = psi1chain,
                            Rcpp::Named("psi2chain")  = psi2chain,
                            Rcpp::Named("vHurdle")    = vHurdle,
                            Rcpp::Named("vZIP")       = vZIP);
}

