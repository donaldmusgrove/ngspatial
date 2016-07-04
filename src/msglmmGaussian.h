#ifndef _ngspatial_MSGLMMGAUSSIAN_H
#define _ngspatial_MSGLMMGAUSSIAN_H

#include <RcppArmadillo.h>
          
         
double dGaussMulti(const arma::mat& Y, const Rcpp::List& X, const Rcpp::List& M, 
                   const Rcpp::List& beta, const arma::vec& Psi, 
                   const arma::vec& tau2e);                                 
Rcpp::List msglmmGaussianCPP(int iterations, const arma::mat& Y, 
                             const Rcpp::List& X, const Rcpp::List& M, 
                             double taub, double nu, double taus, double tau2ea, 
                             double tau2eb, const Rcpp::List& beta0, 
                             const arma::vec& Psi0, const arma::vec& tau2e0, 
                             const arma::mat& Sigma0, int verbose, int runningIter,
                             const arma::vec& fivepct, int maxit);

#endif
