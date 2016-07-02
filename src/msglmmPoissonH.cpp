#include "msglmmPoissonH.h"
#include "util.h"
#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace std;
using namespace arma;


// Density of the multivariate SSGLMM for Poisson data
double dPoissMultiH(const mat& Y, const List& X, const List& M, 
                    const mat& offset, const List& beta, 
                    const mat& Psi, const mat& Delta){

  int J = Y.n_cols;
  vec ll(J);
  
  for(int j = 0; j < J; j++){
    mat Xj      = X(j);
    mat Mj      = M(j);
    vec betaj   = beta(j);
    vec psij    = Psi.col(j);
    vec deltaj  = Delta.col(j);
    vec offsetj = offset.col(j);
    vec eta     = offsetj + Xj*betaj + Mj*psij + Mj*deltaj;
    vec lambda  = exp(eta);
    vec Yj      = Y.col(j);
    ll(j)       = sum( Yj % log(lambda) - lambda);
  }
  
  return(sum(ll));
}


// Sample the fixed effects for a single outcome
vec betaSpoisH(vec beta, double sigmab, double taub, const vec& Y, 
               const mat& X, const mat& M, const vec& offset, 
               const vec& psi, const vec& delta){

  int J  = beta.n_elem;
  vec l1(1); vec l2(1); vec l3(1);
  double uP, ll;
  
  uP           = R::runif(0.0, 1.0);
  vec betaProp = Rcpp::rnorm(J, 0.0, sigmab);
  vec betaP    = betaProp + beta;

  vec eta    = offset + X*beta + M*psi + M*delta;
  vec lambda = exp(eta);
  
  vec etaP    = offset + X*betaP + M*psi + M*delta;
  vec lambdaP = exp(etaP);
  
  l1 = trans(Y) * ( log(lambdaP) - log(lambda) );
  l2 = sum(lambda - lambdaP);
  l3 = 0.5 * taub * ( dot(beta, beta) - dot(betaP, betaP) );
  ll = sum(l1+l2+l3);
  
  if(log(uP) < ll) beta = betaP;
  
  return(beta);
}


// Sample the spatial effects for a single outcome
vec psiSpoisH(const mat& Psi, int num, double sigmap, const vec& Y, 
              const mat& X, const mat& M, const vec& offset, 
              const vec& beta, const vec& delta, const mat& SigmaI){

  int q = Psi.n_rows;

  vec l1(1); vec l2(1); vec l3(1);
  double uP, ll;
  uP = R::runif(0.0, 1.0);

  vec psi     = Psi.col(num);
  vec psiProp = Rcpp::rnorm(q, 0.0, sigmap);
  vec psiP    = psiProp + psi;

  mat PsiP      = Psi;
  PsiP.col(num) = psiP;
  
  vec Psiv  = vectorise(Psi);
  vec PsiPv = vectorise(PsiP);
  
  vec eta     = offset + X*beta + M*psi + M*delta;
  vec lambda  = exp(eta);
  
  vec etaP    = offset + X*beta + M*psiP + M*delta;
  vec lambdaP = exp(etaP);
  
  l1 = trans(Y) * ( log(lambdaP) - log(lambda)   );
  l2 = sum(lambda - lambdaP);
  l3 = 0.5 * ( trans(Psiv) * SigmaI * Psiv - trans(PsiPv) * SigmaI * PsiPv );
  ll = sum(l1+l2+l3);
  
  if(log(uP) < ll) psi = psiP;
    
  return(psi);
}


// Sample the heterogeneous rfx for a single outcome
vec deltaSpois(vec delta, double sigmah, double tauh, const vec& Y, 
               const mat& X, const mat& M, const vec& offset, 
               const vec& beta, const vec& psi){

  int J  = delta.n_elem;
  vec l1(1); vec l2(1); vec l3(1);
  double uP, ll;
  
  uP            = R::runif(0.0, 1.0);
  vec deltaProp = Rcpp::rnorm(J, 0.0, sigmah);
  vec deltaP    = deltaProp + delta;

  vec eta    = offset + X*beta + M*psi + M*delta;
  vec lambda = exp(eta);
  
  vec etaP    = offset + X*beta + M*psi + M*deltaP;
  vec lambdaP = exp(etaP);
  
  l1 = trans(Y) * ( log(lambdaP) - log(lambda) );
  l2 = sum(lambda - lambdaP);
  l3 = 0.5 * tauh * ( dot(delta, delta) - dot(deltaP, deltaP) );
  ll = sum(l1+l2+l3);
  
  if(log(uP) < ll) delta = deltaP;
  
  return(delta);
}


List msglmmPoissonHCPP(int iterations, const mat& Y, const List& X, 
                       const List& M, const mat& offset, 
                       const vec& sigmab, const vec& sigmap, const vec& sigmah, 
                       double taub, double nu, double taus, double tauh,
                       const List& beta0, const mat& Psi0, const mat& Sigma0, 
                       const mat& Delta0,
                       int verbose, int runningIter, const vec& fivepct, int maxit){

  // Setup constants
  int J = Y.n_cols;
  int q = Psi0.n_rows;
  
  //Initialize the chains
  cube PsiChain(q,J,iterations);
  cube SigmaChain(J,J,iterations);
  cube DeltaChain(q,J,iterations);
  mat  aChain(J,iterations+1,fill::zeros);
    
  //Fill the chains with initial values
  field<vec> bChain(iterations+1, J);
  List Betatemp(J);
  for(int j=0; j<J; j++){
    vec btemp   = beta0(j);
    bChain(0,j) = btemp;   
  }

  PsiChain        = join_slices(Psi0,PsiChain);
  SigmaChain      = join_slices(Sigma0,SigmaChain);
  DeltaChain      = join_slices(Delta0,DeltaChain);
  aChain.col(0)   = ones<vec>(J);
    aChain.col(0) = aChain.col(0)*4 / J;
 
  //Initialize covariance matrix
  mat SigInv(J,J);  
  mat SigI(J*q, J*q);
  mat Idq(q,q, fill::eye);

  
  // Initialize DIC calculator
  vec v(iterations+1);
  v(0) = -2 * dPoissMultiH(Y, X, M, offset, beta0, Psi0, Delta0);
  
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
    SigI    = kron(SigInv, Idq);

    mat Psi     = PsiChain.slice(k-1);
    mat Delta   = DeltaChain.slice(k-1);    
    
    for(int j=0; j<J; j++){
      vec betaj    = bChain(k-1,j);
          betaj    = betaSpoisH(betaj, sigmab(j), taub, Y.col(j),  
                               X(j), M(j), offset.col(j), Psi.col(j), Delta.col(j));
      bChain(k,j)  = betaj; 
      Betatemp(j)  = betaj;
      Psi.col(j)   = psiSpoisH(Psi, j, sigmap(j), Y.col(j), X(j), M(j), offset.col(j), bChain(k,j), Delta.col(j), SigI);
      Delta.col(j) = deltaSpois(Delta.col(j), sigmah(j), tauh, Y.col(j), X(j), M(j), offset.col(j), bChain(k,j), Psi.col(j));
      aChain(j,k)  = aSampler(SigInv(j,j),J,nu,taus);
    }
    
    SigmaChain.slice(k) = SigmaSampler(aChain.col(k), Psi, nu);
    v(k)           = -2 * dPoissMultiH(Y, X, M, offset, Betatemp, Psi, Delta);

    PsiChain.slice(k) = Psi;
    
    runningIter = runningIter+1;
  } 
  
  return Rcpp::List::create(Rcpp::Named("bChain")      = bChain,
                            Rcpp::Named("PsiChain")    = PsiChain,
                            Rcpp::Named("DeltaChain")  = DeltaChain,
                            Rcpp::Named("SigmaChain")  = SigmaChain,
                            Rcpp::Named("aChain")      = aChain,
                            Rcpp::Named("v")           = v,
                            Rcpp::Named("runningIter") = runningIter);
}

RCPP_MODULE(msglmmPoissonH){
  Rcpp::function("dPoissMultiH", &dPoissMultiH);
  Rcpp::function("msglmmPoissonHCPP", &msglmmPoissonHCPP);
}
