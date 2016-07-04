#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
// Density of the hurdle model
double dhurdle(vec Y, vec Z, mat X1, mat X2, mat M1, mat M2, vec offset, 
               vec alpha, vec beta, vec phi1, vec phi2){

  int n = Y.n_elem;
   
  vec eta1(n);
  vec eta2(n);
  vec p(n);
  vec lambda(n);
  vec l1(1);
  vec l2(1);
  
  eta1 = X1*alpha + M1*phi1;
  p    = 1 + exp(-eta1);
  p    = 1/p;

  eta2   = offset + X2*beta + M2*phi2;  
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
vec phi1S(vec phi1, double tunep1, vec Z, mat X, mat M, vec alpha, vec phi2,
          mat SigmaI){

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
  
  mat PhiM  = join_cols(phi1,phi2);
  mat PhiPM = join_cols(phi1P,phi2);
  vec Phi   = vectorise(PhiM);
  vec PhiP  = vectorise(PhiPM); 
  
  eta  = X*alpha + M*phi1;
  p    = 1 + exp(-eta);
  p    = 1/p;
  
  etaP = X*alpha + M*phi1P;
  pP   = 1 + exp(-etaP);
  pP   = 1/pP;
  
  l1 = trans(Z)   * ( log(pP) - log(p) );
  l2 = trans(1-Z) * ( log(1-pP) - log(1-p) );
  l3 = 0.5 * ( trans(Phi) * SigmaI * Phi - trans(PhiP) * SigmaI * PhiP);
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
          vec phi1, mat SigmaI){

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

  mat PhiM  = join_cols(phi1,phi2);
  mat PhiPM = join_cols(phi1,phi2P);
  vec Phi   = vectorise(PhiM);
  vec PhiP  = vectorise(PhiPM); 
  
  eta    = logE + X*beta  + M*phi2;
  lambda = exp(eta);
  
  etaP    = logE + X*beta + M*phi2P;
  lambdaP = exp(etaP);
  
  l1 = trans(Z)   * ( log(1-exp(-lambda)) - log(1-exp(-lambdaP)) );
  l2 = trans(Y%Z) * ( log(lambdaP) - log(lambda)   );
  l3 = trans(Z)   * (lambda - lambdaP);
  l4 = 0.5 * ( trans(Phi) * SigmaI * Phi - trans(PhiP) * SigmaI * PhiP);
  ll = sum(l1+l2+l3+l4);
  
  if(log(uP) < ll) phi2 = phi2P;
  
  return(phi2);
}

// [[Rcpp::export]]
double qlogisCPP(double rho){
  return( log(rho/(1-rho)) );
}

// [[Rcpp::export]]
double plogisCPP(double x){
  return( 1/(1+exp(-x)) );
}


// [[Rcpp::export]]
vec sSampler(vec S, vec tunes, vec phi1, vec phi2){

  //Initialize variables
  double s1  = S(0);
  double s2  = S(1);
  double rho = S(2);

  int q      = phi1.n_elem;  
  double e1, e2, e3, e4, ee, uP, s1T, s1TP, s1P;
  double s2T, s2TP, s2P;
  double rhoT, rhoTP, rhoP;
  vec s1prop(1);
  vec s2prop(1);
  vec rhoprop(1);
  
  //Update s1 
  s1T    = log(s1);
  s1prop = Rcpp::rnorm(1, 0.0, tunes(0));
  s1TP   = s1T + s1prop(0);
  s1P    = exp(s1TP);
  uP     = R::runif(0.0, 1.0);
  
  e1 = q/2 * (log(s1) - log(s1P));
  e2 = 0.5 / (1-rho*rho) * sum( phi1%phi1 / s1  - 2 * rho * phi1 % phi2 / sqrt(s1*s2) );
  e3 = 0.5 / (1-rho*rho) * sum( phi1%phi1 / s1P - 2 * rho * phi1 % phi2 / sqrt(s1P*s2) );  
  e4 = log(s1P) - log(s1);
  ee = e1 + e2 - e3 + e4;
  
  if(log(uP) < ee) s1 = s1P;

  //Update s2  
  s2T    = log(s2);
  s2prop = Rcpp::rnorm(1, 0.0, tunes(1));
  s2TP   = s2T + s2prop(0);
  s2P    = exp(s2TP);
  uP     = R::runif(0.0, 1.0);
  
  e1 = q/2 * (log(s2) - log(s2P));
  e2 = 0.5 / (1-rho*rho) * sum( phi2%phi2 / s2  - 2 * rho * phi1 % phi2 / sqrt(s1*s2) );
  e3 = 0.5 / (1-rho*rho) * sum( phi2%phi2 / s2P - 2 * rho * phi1 % phi2 / sqrt(s1*s2P) );  
  e4 = log(s2P) - log(s2);
  ee = e1 + e2 - e3 + e4;
  
  if(log(uP) < ee) s2 = s2P;
  
  //Update rho  
  rhoT    = qlogisCPP((rho+1)/2);
  rhoprop = Rcpp::rnorm(1, 0.0, tunes(2));
  rhoTP   = rhoT + rhoprop(0);
  rhoP    = 2 * plogisCPP(rhoTP) - 1;
  uP      = R::runif(0.0, 1.0);
  
  e1 = q/2 * ( log(1-rho*rho) - log(1-rhoP*rhoP) );
  e2 = 0.5 / (1-rho*rho)   * sum( phi1%phi1/s1 + phi2%phi2/s2  - 2 * rho * phi1 % phi2 / sqrt(s1*s2) );
  e3 = 0.5 / (1-rhoP*rhoP) * sum( phi1%phi1/s1 + phi2%phi2/s2  - 2 * rhoP * phi1 % phi2 / sqrt(s1*s2) );  
  e4 = log(rhoP + 1) + log(1-rhoP) - log(rho+1) - log(1-rho);
  ee = e1 + e2 - e3 + e4;
  
  if(log(uP) < ee) rho = rhoP;
  
  S(0) = s1;
  S(1) = s2;
  S(2) = rho;
  
  return(S);
}



// [[Rcpp::export]]
List MCMCzip(int iterations, vec Y, vec Z, mat X1, mat X2, mat M1, mat M2,
             vec logE, double tunea, double tuneb, double tunep1, double tunep2, vec tunes){

  int p1 = X1.n_cols;
  int p2 = X2.n_cols;
  int q1 = M1.n_cols;  
  int q2 = M2.n_cols; 
  
  // Initialize the chains
  mat aChain(p1, iterations);
    aChain.col(0) = randn<vec>(p1);
  mat bChain(p2, iterations);
    bChain.col(0) = randn<vec>(p2);  
  mat p1Chain(q1, iterations);
    p1Chain.col(0) = randn<vec>(q1);
  mat p2Chain(q2, iterations);
    p2Chain.col(0) = randn<vec>(q2);
  mat sChain(3, iterations);
    sChain.col(0) = randu<vec>(3);	

  // Initialize covariance matrix
  mat S(2,2, fill::zeros);
  mat R(2,2, fill::eye);
  mat SigInv(2,2);
  mat SigI(2*q1, 2*q1);
  mat Idq(q1,q1, fill::eye);
  
  // Initialize DIC calculator
  vec v(iterations);
  v(0) = -2 * dhurdle(Y, Z, X1, X2, M1, M2, logE, aChain.col(0), bChain.col(0), p1Chain.col(0), p2Chain.col(0));

  //MCMC iterations
  for(int k=1; k < iterations; k++ ){  
    S(0,0) = sqrt(sChain(0,k-1));
    S(1,1) = sqrt(sChain(1,k-1));	
    R(0,1) = sChain(2,k-1);
    R(1,0) = sChain(2,k-1);
    SigInv = S * R * S;
    SigInv = inv(SigInv);
    SigI   = kron(SigInv, Idq);
   
    aChain.col(k)  = alphaS(aChain.col(k-1), tunea, Z, X1, M1, p1Chain.col(k-1));
    bChain.col(k)  = betaS(bChain.col(k-1), tuneb, Y, Z, X2, M2, logE, p2Chain.col(k-1));
    p1Chain.col(k) = phi1S(p1Chain.col(k-1), tunep1, Z, X1, M1, aChain.col(k), p2Chain.col(k-1), SigI);
    p2Chain.col(k) = phi2S(p2Chain.col(k-1), tunep2, Y, Z, X2, M2, logE, bChain.col(k), p1Chain.col(k), SigI);
    sChain.col(k)  = sSampler(sChain.col(k-1), tunes, p1Chain.col(k), p2Chain.col(k));
    v(k)           = -2 * dhurdle(Y, Z, X1, X2, M1, M2, logE, aChain.col(k), bChain.col(k), p1Chain.col(k), p2Chain.col(k));   
  }
   
  return List::create(aChain, bChain, p1Chain, p2Chain, sChain, v);
}
