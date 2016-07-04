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
mat SigmaSampler(vec a, mat Psi, double nu){
  int J     = Psi.n_rows;
  int q     = Psi.n_cols;
  
  mat PP    = trans(Psi) * Psi;
  mat Adiag = diagmat(a);
  
  mat SigmaHat = PP + 2 * nu * Adiag;
  
  return riwishart(nu+q+J-1, SigmaHat);
}

// Draw a sample from the posterior of the diagonal of the within-areal unit covariance matrix
double aSampler(double Sinv, int J, double nu, double taus){
  //double taus  = pow(10,-10);
  //double nu = 2;

  double aShp = (nu+J)/2;
  double aRte = nu * Sinv + taus;
  
  return(1/R::rgamma(aShp, 1/aRte));
}


// Density of the multivariate SSGLMM for Binomial data
// [[Rcpp::export]]
double dBinomMulti(mat Y, mat Ntrials, List X, List M, mat offset, List beta, 
                   mat Psi){

  int J  = Y.n_cols;
  vec ll(J);
  
  for(int j = 0; j < J; j++){
    mat Xj       = X(j);
    mat Mj       = M(j);
    vec betaj    = beta(j);
    vec psij     = Psi.col(j);
    vec offsetj  = offset.col(j);
    vec eta      = offsetj + Xj*betaj + Mj*psij;
    vec lambda   = 1+exp(-eta);
        lambda   = 1/lambda;
    ll(j)        = sum(trans(Y.col(j)) * log(lambda) + trans(Ntrials.col(j) - Y.col(j)) * log(1-lambda));
  }
  
  return(sum(ll));
}

// Sample the fixed effects for a single outcome
vec betaSbinom(vec beta, double sigmab, double taub, vec Y, vec Ntrials, mat X, mat M, 
               vec offset, vec psi){

  int J  = beta.n_elem;
  vec l1(1); vec l2(1); vec l3(1);
  double uP, ll;
  
  uP           = R::runif(0.0, 1.0);
  vec betaProp = Rcpp::rnorm(J, 0.0, sigmab);
  vec betaP    = betaProp + beta;

  vec eta    = offset + X*beta + M*psi;
  vec lambda = 1 + exp(-eta);
      lambda = 1/lambda;
  
  vec etaP    = offset + X*betaP  + M*psi;
  vec lambdaP = 1 + exp(-etaP);
      lambdaP = 1/lambdaP;
  
  l1 = trans(Y) * ( log(lambdaP) - log(lambda)   );
  l2 = trans(Ntrials-Y) * ( log(1-lambdaP) - log(1-lambda) );
  l3 = 0.5 * taub * ( dot(beta, beta) - dot(betaP, betaP)   );
  ll = sum(l1+l2+l3);
  
  if(log(uP) < ll) beta = betaP;
  
  return(beta);
}


// Sample the spatial effects for a single outcome
vec psiSbinom(mat Psi, int num, double sigmap, vec Y, vec Ntrials, mat X, mat M, 
             vec offset, vec beta, mat SigmaI){

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
  
  vec eta    = offset + X*beta + M*psi;
  vec lambda = 1 + exp(-eta);
  lambda     = 1/lambda;
  
  vec etaP    = offset + X*beta + M*psiP;
  vec lambdaP = 1 + exp(-etaP);
  lambdaP     = 1/lambdaP;
   
  l1 = trans(Y) * ( log(lambdaP) - log(lambda)   );
  l2 = trans(Ntrials-Y) * ( log(1-lambdaP) - log(1-lambda) );
  l3 = 0.5 * ( trans(Psiv) * SigmaI * Psiv - trans(PsiPv) * SigmaI * PsiPv );
  ll = sum(l1+l2+l3);
  
  if(log(uP) < ll) psi = psiP;
    
  return(psi);
}


bool isin(int num, vec vals){
  int n = vals.n_elem;
  
  bool isequal = false;

  for (int i = 0; i < n; i++){   
  
    if(abs(vals(i) - num) < 0.001){
      isequal = true;
      break;
    }
  
  }
  return(isequal);
}

// [[Rcpp::export]]
List msglmmBinomial(int iterations, mat Y, mat Ntrials, List X, List M, mat offset, 
                    vec sigmab, vec sigmap, double taub, double nu, double taus,
                    List beta0, mat Psi0, mat Sigma0, 
                    int verbose, int runningIter, vec fivepct, int maxit){

  // Setup constants
  int J  = Y.n_cols;
  int q  = Psi0.n_rows;
  
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

  
  // Initialize DIC calculator
  vec v(iterations+1);
  v(0) = -2 * dBinomMulti(Y, Ntrials, X, M, offset, beta0, Psi0);

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
    
    for(int j=0; j<J; j++){
      vec betaj   = bChain(k-1,j);
          betaj   = betaSbinom(betaj, sigmab(j), taub, Y.col(j), Ntrials.col(j), 
                               X(j), M(j), offset.col(j), Psi.col(j));
      bChain(k,j) = betaj; 
      Betatemp(j) = betaj;  
      Psi.col(j)  = psiSbinom(Psi, j, sigmap(j), Y.col(j), Ntrials.col(j), X(j), M(j), offset.col(j), bChain(k,j), SigI);
      aChain(j,k) = aSampler(SigInv(j,j),J,nu,taus);
    }
        
    SigmaChain.slice(k) = SigmaSampler(aChain.col(k), Psi, nu);
    v(k)                = -2 * dBinomMulti(Y, Ntrials, X, M, offset, Betatemp, Psi);

    PsiChain.slice(k) = Psi;
    
    runningIter = runningIter+1;
  } 
    
  return Rcpp::List::create(Rcpp::Named("bChain")      = bChain,
                            Rcpp::Named("PsiChain")    = PsiChain,
                            Rcpp::Named("SigmaChain")  = SigmaChain,
                            Rcpp::Named("aChain")      = aChain,
                            Rcpp::Named("v")           = v,
                            Rcpp::Named("runningIter") = runningIter);
}
