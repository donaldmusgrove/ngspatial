#include <RcppArmadillo.h>
#include <utils.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// Density of the hurdle model 
// [[Rcpp::export]]
double dHurdle(mat Y, List X, List M, mat offset, List beta, mat Psi){
   
  // Declare input
  mat X1 = X(0); mat M1 = M(0);
  mat X2 = X(1); mat M2 = M(1);
  vec beta0 = beta(0); vec beta1 = beta(1); 
  vec Psi0 = Psi.col(0); vec Psi1 = Psi.col(1);
  
  vec z = Y.col(0);
  vec y = Y.col(1);
  
  // Binary piece
  vec eta1 = offset.col(0) + X1*beta0 + M1*Psi0;
  vec p    = 1 + exp(-eta1);
      p    = 1/p;

  // Count piece
  vec eta2   = offset.col(1) + X2*beta1 + M2*Psi1;  
  vec lambda = exp(eta2);

  double l1 = sum(trans(1-z) * log(1 - p));
  double l2 = sum(trans(z) * ( log(p) - log(1-exp(-lambda)) + y % log(lambda) - lambda ));
  
  return(sum(l1+l2));
}


// ZIP model density 
// [[Rcpp::export]]
double dZIP(mat Y, List X, List M, mat offset, List beta, mat Psi){

  mat X1 = X(0); mat M1 = M(0);
  mat X2 = X(1); mat M2 = M(1);
  vec beta0 = beta(0); vec beta1 = beta(1); 
  vec Psi0 = Psi.col(0); vec Psi1 = Psi.col(1);
  
  vec z = Y.col(0);
  vec y = Y.col(1);
  
  // Binary piece
  vec eta1 = offset.col(0) + X1*beta0 + M1*Psi0;
  vec p    = 1 + exp(-eta1);
      p    = 1/p;

  // Count piece
  vec eta2   = offset.col(1) + X2*beta1 + M2*Psi1;  
  vec lambda = exp(eta2);

  //Reparameterize p for the ZIP model, i.e., p* = p(1-exp(-lambda))
  p = p%(1-exp(-lambda));
  
  double l1 = sum(trans(1-z) * log(1 - p));
  double l2 = sum(trans(z) * ( log(p) - log(1-exp(-lambda)) + y % log(lambda) - lambda ));
  
  return(sum(l1+l2));
}


vec beta1Shurd(vec beta1, double taub, double tuneb1, vec z, mat X, 
               mat M, vec offset, vec psi1){

  int jj  = beta1.n_elem;
  vec l1(1); vec l2(1); vec l3(1); 
  double ll;
  
  double uP     = R::runif(0.0, 1.0);
  vec beta1Prop = Rcpp::rnorm(jj, 0.0, tuneb1);
  vec beta1P    = beta1Prop + beta1;

  vec eta  = offset + X*beta1  + M*psi1;
  vec p    = 1 + exp(-eta);
      p    = 1/p;
  
  vec etaP = offset + X*beta1P + M*psi1;
  vec pP   = 1 + exp(-etaP);
      pP   = 1/pP;
  
  l1 = trans(z)   * ( log(pP) - log(p) );
  l2 = trans(1-z) * ( log(1-pP) - log(1-p) );
  l3 = 0.5 * ( dot(beta1, beta1) - dot(beta1P, beta1P)  ) / taub;
  ll = sum(l1+l2+l3);
  
  if(log(uP) < ll) beta1 = beta1P;
  
  return(beta1);
}

vec psi1Shurd(mat Psi, double tunep1, vec z, mat X, mat M, vec offset, 
              vec beta1, mat SigmaI){

  int q  = Psi.n_rows;
  
  vec l1(1); vec l2(1); vec l3(1);
  double uP, ll;
  uP = R::runif(0.0, 1.0);
  
  vec psi1     = Psi.col(0);
  vec psi2     = Psi.col(1);
   
  vec psi1Prop = Rcpp::rnorm(q, 0.0, tunep1);
  vec psi1P    = psi1Prop + psi1;
  
  mat PsiM  = join_cols(psi1,psi2);
  mat PsiPM = join_cols(psi1P,psi2);
  vec Psiv  = vectorise(PsiM);
  vec PsivP = vectorise(PsiPM);
    
  vec eta  = offset + X*beta1 + M*psi1;
  vec p    = 1 + exp(-eta);
      p    = 1/p;
  
  
  vec etaP = offset + X*beta1 + M*psi1P;
  vec pP   = 1 + exp(-etaP);
      pP   = 1/pP;
  
  
  l1 = trans(z)   * ( log(pP) - log(p) );
  l2 = trans(1-z) * ( log(1-pP) - log(1-p) );
  l3 = 0.5 * ( trans(Psiv) * SigmaI * Psiv - trans(PsivP) * SigmaI * PsivP);
  ll = sum(l1+l2+l3);
  
  if(log(uP) < ll) psi1 = psi1P;
  
  return(psi1);
}


vec beta2Shurd(vec beta2, double taub, double tuneb2, vec z, vec y, mat X, 
               mat M, vec offset, vec psi2){

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


vec psi2Shurd(mat Psi, double tunep2, vec z, vec y, mat X, mat M, vec offset, vec beta2, 
          mat SigmaI){

  int q = Psi.n_rows;

  vec l1(1); vec l2(1); vec l3(1); vec l4(1);
  double uP, ll;
  uP = R::runif(0.0, 1.0);
  
  vec psi1     = Psi.col(0);
  vec psi2     = Psi.col(1);
  
  vec psi2Prop = Rcpp::rnorm(q, 0.0, tunep2);
  vec psi2P    = psi2Prop + psi2;

  mat PsiM  = join_cols(psi1,psi2);
  mat PsiPM = join_cols(psi1,psi2P);
  vec Psiv   = vectorise(PsiM);
  vec PsivP  = vectorise(PsiPM); 
  
  vec eta    = offset + X*beta2  + M*psi2;
  vec lambda = exp(eta);
  
  vec etaP    = offset + X*beta2 + M*psi2P;
  vec lambdaP = exp(etaP);
  
  l1 = trans(z)   * ( log(1-exp(-lambda)) - log(1-exp(-lambdaP)) );
  l2 = trans(y%z) * ( log(lambdaP) - log(lambda)   );
  l3 = trans(z)   * (lambda - lambdaP);
  l4 = 0.5 * ( trans(Psiv) * SigmaI * Psiv - trans(PsivP) * SigmaI * PsivP);
  ll = sum(l1+l2+l3+l4);
  
  if(log(uP) < ll) psi2 = psi2P;
  
  return(psi2);
}



// [[Rcpp::export]]
List msglmmHurdle(int iterations, mat Y, List X, List M, mat offset,
                  vec sigmab, vec sigmap, double taub, double nu, double taus,
                  List beta0, mat Psi0, mat Sigma0,
                  int verbose, int runningIter, vec fivepct, int maxit){

                  
  // Setup constants
  int J = Y.n_cols;
  int q = Psi0.n_rows;
  
  //Initialize the chains
  cube PsiChain(q,J,iterations);
  cube SigmaChain(J,J,iterations);
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
  aChain.col(0)   = ones<vec>(J);
    aChain.col(0) = aChain.col(0)*4 / J;
    
  
  //Initialize covariance matrix
  mat SigInv(J,J);  
  mat SigI(J*q, J*q);
  mat Idq(q,q, fill::eye);
  
 
  // Initialize DIC calculator for Hurdle and ZIP models
  vec vHurdle(iterations+1);
  vec vZIP(iterations+1);
  vHurdle(0) = -2 * dHurdle(Y, X, M, offset, beta0, Psi0);
  vZIP(0)    = -2 * dZIP(Y, X, M, offset, beta0, Psi0);
   
  
  Rcpp::checkUserInterrupt();
      
  //MCMC iterations
  for(int k=1; k < iterations+1; k++ ){
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
           
    vec beta1 = bChain(k-1, 0);
        beta1 = beta1Shurd(beta1, taub, sigmab(0), Y.col(0), X(0), M(0), offset.col(1), Psi.col(0));
    bChain(k,0) = beta1;
    Betatemp(0) = beta1;
    
    
    vec beta2 = bChain(k-1,1);
        beta2 = beta2Shurd(beta2, taub, sigmab(1), Y.col(0), Y.col(1), X(1), M(1), offset.col(1), Psi.col(1));
    bChain(k,1) = beta2;
    Betatemp(1) = beta2;
       
    Psi.col(0) = psi1Shurd(Psi, sigmap(0), Y.col(0), X(0), M(0), offset.col(1), beta1, SigI);  
    Psi.col(1) = psi2Shurd(Psi, sigmap(1), Y.col(0), Y.col(1), X(1), M(1), offset.col(1), beta1, SigI);
    PsiChain.slice(k) = Psi;    
      
    aChain(0,k)   = aSampler(SigInv(0,0),J,nu,taus);
    aChain(1,k)   = aSampler(SigInv(1,1),J,nu,taus);
        
    SigmaChain.slice(k) = SigmaSampler(aChain.col(k), Psi, nu);
      
    vHurdle(k) = -2 * dHurdle(Y, X, M, offset, Betatemp, Psi);    
    vZIP(k)    = -2 * dZIP(Y, X, M, offset, Betatemp, Psi);
    
    runningIter = runningIter+1;
  }
  
  return Rcpp::List::create(Rcpp::Named("bChain")      = bChain,
                            Rcpp::Named("PsiChain")    = PsiChain,
                            Rcpp::Named("SigmaChain")  = SigmaChain,
                            Rcpp::Named("aChain")      = aChain,
                            Rcpp::Named("vHurdle")     = vHurdle,
                            Rcpp::Named("vZIP")        = vZIP,
                            Rcpp::Named("runningIter") = runningIter);
  
}

