#ifndef _ngspatial_MSGLMMPOISSONH_H
#define _ngspatial_MSGLMMPOISSONH_H

#include <RcppArmadillo.h>
          
         
double dPoissMultiH(const arma::mat& Y, const Rcpp::List& X, 
                   const Rcpp::List& M, const arma::mat& offset, 
                   const Rcpp::List& beta, const arma::mat& Psi,
                   const arma::mat& Delta);                                     
Rcpp::List msglmmPoissonHCPP(int iterations, const arma::mat& Y, 
                            const Rcpp::List& X, const Rcpp::List& M, 
                            const arma::mat& offset, const arma::vec& sigmab, 
                            const arma::vec& sigmap, const arma::vec& sigmah,
                            double taub, double nu, double taus, double tauh,
                            const Rcpp::List& beta0, const arma::mat& Psi0, 
                            const arma::mat& Sigma0, const arma::mat& Delta0,
                            int verbose, int runningIter,
                            const arma::vec& fivepct, int maxit);

#endif
