#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
// Density of the hurdle model
double dhurdle(vec Y, vec Z, mat X, mat M, vec offset, vec alpha, vec beta,
               vec phi1, vec phi2){

  int n = Y.n_elem;
   
  vec eta1(n);
  vec eta2(n);
  vec p(n);
  vec lambda(n);
  vec l1(1);
  vec l2(1);
  
  eta1 = X*alpha + M*phi1;
  p    = 1 + exp(-eta1);
  p    = 1/p;

  eta2   = offset + X*beta + M*phi2;  
  lambda = exp(eta2);

  l1 = trans(1-Z) * log(1 - p);
  l2 = trans(Z) * ( log(p) - log(1-exp(-lambda)) + Y % log(lambda) - lambda );
  
  return(sum(l1+l2));
}


// [[Rcpp::export]]
vec alphaS(vec alpha, double tunea, vec Z, mat X, mat M, vec phi1){

  int J  = alpha.n_elem;
  vec alphaP(J);
  vec alphaProp(J);
  vec eta(J);
  vec etaP(J);
  vec p(J);
  vec pP(J);
  vec l1(1); vec l2(1); vec l3(1); 
  double uP, ll;
  
  uP        = R::runif(0.0, 1.0);
  alphaProp = Rcpp::rnorm(J, 0.0, tunea);
  
  alphaP    = alphaProp + alpha;

  eta  = X*alpha  + M*phi1;
  p    = 1 + exp(-eta);
  p    = 1/p;
  
  etaP = X*alphaP + M*phi1;
  pP   = 1 + exp(-etaP);
  pP   = 1/pP;
  
  l1 = trans(Z)   * ( log(pP) - log(p) );
  l2 = trans(1-Z) * ( log(1-pP) - log(1-p) );
  l3 = 0.5 * ( dot(alpha, alpha) - dot(alphaP, alphaP)   ) / 10000;
  ll = sum(l1+l2+l3);
  
  if(log(uP) < ll) alpha = alphaP;
  
  return(alpha);
}


// [[Rcpp::export]]
vec phi1S(vec phi1, double tunep1, vec Z, mat X, mat M, vec alpha, double s1){

  int J  = phi1.n_elem;
  vec phi1P(J);
  vec phi1Prop(J);
  vec eta(J);
  vec etaP(J);
  vec p(J);
  vec pP(J);
  vec l1(1); vec l2(1); vec l3(1); 
  double uP, ll;
  
  uP       = R::runif(0.0, 1.0);
  phi1Prop = Rcpp::rnorm(J, 0.0, tunep1);
  phi1P    = phi1Prop + phi1;
   
  eta  = X*alpha  + M*phi1;
  p    = 1 + exp(-eta);
  p    = 1/p;
  
  etaP = X*alpha + M*phi1P;
  pP   = 1 + exp(-etaP);
  pP   = 1/pP;
  
  l1 = trans(Z)   * ( log(pP) - log(p) );
  l2 = trans(1-Z) * ( log(1-pP) - log(1-p) );
  l3 = 0.5 * ( dot(phi1, phi1) -  dot(phi1P, phi1P)) / s1;
  ll = sum(l1+l2+l3);
  
  if(log(uP) < ll) phi1 = phi1P;
  
  return(phi1);
}


// [[Rcpp::export]]
vec betaS(vec beta, double tuneb, vec Y, vec Z, mat X, mat M, vec logE, 
          vec phi2){

  int J  = beta.n_elem;
  vec betaP(J);
  vec betaProp(J);
  vec eta(J);
  vec etaP(J);
  vec lambda(J);
  vec lambdaP(J);
  vec l1(1); vec l2(1); vec l3(1);  vec l4(1);
  double uP, ll;
  
  uP       = R::runif(0.0, 1.0);
  betaProp = Rcpp::rnorm(J, 0.0, tuneb);
  betaP    = betaProp + beta;

  eta    = logE + X*beta  + M*phi2;
  lambda = exp(eta);
  
  etaP    = logE + X*betaP  + M*phi2;
  lambdaP = exp(etaP);
  
  l1 = trans(Z)   * ( log(1-exp(-lambda)) - log(1-exp(-lambdaP)) );
  l2 = trans(Y%Z) * ( log(lambdaP) - log(lambda)   );
  l3 = trans(Z)   * (lambda - lambdaP);
  l4 = 0.5 * ( dot(beta, beta) - dot(betaP, betaP)   ) / 10000;
  ll = sum(l1+l2+l3+l4);
  
  if(log(uP) < ll) beta = betaP;
  
  return(beta);
}


// [[Rcpp::export]]
vec phi2S(vec phi2, double tunep2, vec Y, vec Z, mat X, mat M, vec logE, vec beta, 
          double s2){

  int J  = phi2.n_elem;
  vec phi2P(J);
  vec phi2Prop(J);
  vec eta(J);
  vec etaP(J);
  vec lambda(J);
  vec lambdaP(J);
  vec l1(1); vec l2(1); vec l3(1); vec l4(1);
  double uP, ll;
  
  uP       = R::runif(0.0, 1.0);
  phi2Prop = Rcpp::rnorm(J, 0.0, tunep2);
  phi2P    = phi2Prop + phi2;
   
  eta    = logE + X*beta  + M*phi2;
  lambda = exp(eta);
  
  etaP    = logE + X*beta + M*phi2P;
  lambdaP = exp(etaP);
  
  l1 = trans(Z)   * ( log(1-exp(-lambda)) - log(1-exp(-lambdaP)) );
  l2 = trans(Y%Z) * ( log(lambdaP) - log(lambda)   );
  l3 = trans(Z)   * (lambda - lambdaP);
  l4 = 0.5 * ( dot(phi2, phi2) -  dot(phi2P, phi2P)) / s2;
  ll = sum(l1+l2+l3+l4);
  
  if(log(uP) < ll) phi2 = phi2P;
  
  return(phi2);
}


// [[Rcpp::export]]
double sSampler(vec phi){

  int q = phi.n_elem;
  double sShp  = q/2 + 0.001;
  double sRate = 0.5 * dot(phi,phi) + 0.001;
  
  return(1/R::rgamma(sShp, 1/sRate));
}



// [[Rcpp::export]]
List MCMCzip(int iterations, vec Y, vec Z, mat X, mat M, vec logE,
             double tunea, double tuneb, double tunep1, double tunep2){

  // Setup constants
  int q = M.n_cols;
  int p = X.n_cols;
  
  // Initialize the chains
  mat aChain(p, iterations);
    aChain.col(0) = randn<vec>(p);
  mat bChain(p, iterations);
    bChain.col(0) = randn<vec>(p);  
  mat p1Chain(q, iterations);
    p1Chain.col(0) = randn<vec>(q);
  mat p2Chain(q, iterations);
    p2Chain.col(0) = randn<vec>(q);
  mat sChain(2, iterations);
    sChain.col(0) = randu<vec>(2);	

  // Initialize DIC calculator
  vec v(iterations);
  v(0) = -2 * dhurdle(Y, Z, X, M, logE, aChain.col(0), bChain.col(0), p1Chain.col(0), p2Chain.col(0));

  //MCMC iterations
  for(int k=1; k < iterations; k++ ){
    aChain.col(k)  = alphaS(aChain.col(k-1), tunea, Z, X, M, p1Chain.col(k-1));
    bChain.col(k)  = betaS(bChain.col(k-1), tuneb, Y, Z, X, M, logE, p2Chain.col(k-1));
	p1Chain.col(k) = phi1S(p1Chain.col(k-1), tunep1, Z, X, M, aChain.col(k), sChain(0,k-1));
	p2Chain.col(k) = phi2S(p2Chain.col(k-1), tunep2, Y, Z, X, M, logE, bChain.col(k), sChain(1,k-1));
	sChain(0,k)    = sSampler(p1Chain.col(k));
	sChain(1,k)    = sSampler(p2Chain.col(k));
	v(k)           = -2 * dhurdle(Y, Z, X, M, logE, aChain.col(k), bChain.col(k), p1Chain.col(k), p2Chain.col(k));
  }
   
  return List::create(aChain, bChain, p1Chain, p2Chain, sChain, v);
}
