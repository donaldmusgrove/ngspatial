#ifndef _ngspatial_MSGLMMHURDLE_H
#define _ngspatial_MSGLMMHURDLE_H

#include <RcppArmadillo.h>
             
double dHurdle(const arma::mat& Y, const Rcpp::List& X, const Rcpp::List& M, 
               const arma::mat& offset, const Rcpp::List& beta, 
               const arma::mat& Psi);
double dZIP(const arma::mat& Y, const Rcpp::List& X, const Rcpp::List& M, 
            const arma::mat& offset, const Rcpp::List& beta, const arma::mat& Psi);               
Rcpp::List msglmmHurdleCPP(int iterations, const arma::mat& Y, const Rcpp::List& X, 
                     const Rcpp::List& M, const arma::mat& offset,
                     const arma::vec& sigmab, const arma::vec& sigmap, 
                     double taub, double nu, double taus,
                     const Rcpp::List& beta0, const arma::mat& Psi0, const arma::mat& Sigma0,
                     int verbose, int runningIter, const arma::vec& fivepct, 
                     int maxit);

#endif
