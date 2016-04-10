#ifndef _ngspatial_UTIL_H
#define _ngspatial_UTIL_H

#include <RcppArmadillo.h>

arma::mat rwishart(int df, const arma::mat& S);
arma::mat riwishart(int df, const arma::mat & S);
arma::mat SigmaSampler(const arma::vec& a, const arma::mat& Psi);
double aSampler(double Sinv, int J);
arma::vec bmcpp(const arma::vec& X);
arma::mat bmmatcpp(const arma::mat& M);
arma::cube arraybind(const arma::cube& M1, const arma::cube& M2);

#endif
