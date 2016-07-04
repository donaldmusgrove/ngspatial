#ifndef _ngspatial_MSGLMMUTIL_H
#define _ngspatial_MSGLMMUTIL_H

#include <RcppArmadillo.h>

arma::mat rwishart(int df, const arma::mat& S);
arma::mat riwishart(int df, const arma::mat& S);
arma::mat SigmaSampler(const arma::vec& a, const arma::mat& Psi, double nu);
double aSampler(double Sinv, int J, double nu, double taus);
arma::vec bmcpp(const arma::vec& X);
arma::mat bmmatcpp(const arma::mat& M);
arma::cube arraybind(const arma::cube& M1, const arma::cube& M2);
bool isin(int num, const arma::vec& vals);
arma::mat blockDiagT(const Rcpp::List& M);
arma::mat cov2corcpp(const arma::mat& V);
arma::cube acov2corcpp(const arma::cube& V);
arma::vec vectorSubset(const arma::vec& x, int a, int b);
arma::mat SigmaSamplerVec(const arma::vec& a, const arma::vec& Psi, double nu);
arma::mat blockDiag(const Rcpp::List& M);

#endif
