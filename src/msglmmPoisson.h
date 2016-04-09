#ifndef _ngspatial_MSGLMMPOISSON_H
#define _ngspatial_MSGLMMPOISSON_H

#include <RcppArmadillo.h>

double dPoissMulti(const arma::mat& Y, const Rcpp::List& X, const Rcpp::List& M, const arma::mat& offset, const arma::mat& beta, const arma::mat& Psi);
arma::vec betaSpois(const arma::vec& beta, double sigmab, double taub, const arma::vec& Y, const arma::mat& X, const arma::mat& M, const arma::vec& offset, const arma::vec& psi);
arma::vec psiSpois(const arma::mat& Psi, int num, double sigmas, const arma::vec& Y, const arma::mat& X, const arma::mat& M, const arma::vec& offset, const arma::vec& beta, const arma::mat& SigmaI);
Rcpp::List msglmmPoissonCPP(int iterations, const arma::mat& Y, const Rcpp::List& X, const Rcpp::List& M, const arma::mat& offset, const arma::vec& sigmab, const arma::vec& sigmas, double taub, const arma::mat& beta0, const arma::mat& Psi0, const arma::mat& Sigma0, int verbose);

#endif
