#ifndef _ngspatial_UTIL_H
#define _ngspatial_UTIL_H

#include <RcppArmadillo.h>

arma::mat rwishart(int df, const arma::mat& S);
arma::mat riwishart(int df, const arma::mat & S);
arma::mat SigmaSampler(const arma::vec& a, const arma::mat& Psi);
double aSampler(double Sinv, int J);

#endif
