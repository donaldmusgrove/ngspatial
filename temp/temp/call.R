library(RcppArmadillo)

setwd("D://Dropbox//github//ngspatial//temp")
setwd("C://Dropbox//github//ngspatial//temp")


Rcpp::sourceCpp('temp/msglmmBinomial.cpp')
Rcpp::sourceCpp('msglmmPoisson.cpp')
Rcpp::sourceCpp('msglmmGaussian.cpp')
Rcpp::sourceCpp('msglmmHurdle.cpp')


setwd("C:/Dropbox/github")
devtools::document("ngspatial")
detach("package:ngspatial", unload=TRUE)
install.packages("ngspatial_1.1.tar.gz", repos = NULL, type="source")
library(ngspatial)



n <- 500
J <- 3
q <- 25

Y       <- matrix(rbinom(n*J, 6, 0.8),ncol=J)
Ntrials <- matrix(6,n,J)
X       <- cbind(rep(1,n),seq(0,6,,n))
  X     <- lapply(1:J, function(x) X) 
M       <- matrix(rnorm(n*q),ncol=q)
  M     <- lapply(1:J, function(x) M) 
offset  <- matrix(1,n,J)
beta    <- matrix(1:(2*J)/(2*J),ncol=J)
Phi     <- matrix(1:(J*q)/(J*q),ncol=J)


hh <- sparse.msglmm(family="binomial",1e3, Y, X, M, offset, rep(0.001,J), rep(0.1,J), 1/1000,
              beta, Phi, diag(J), verbose=1)
plot(sapply(hh$bChain, function(x) x[1,3]))



hh <- msglmmBinomial(1e4, Y, Ntrials, X, M, offset, 
                    rep(0.1,J), rep(0.1,J), 1/1000,
                    beta, Phi, diag(J), 1)


h1 <- msglmmPoisson(1e4, Y, X, M, offset, 
                    rep(0.1,J), rep(0.1,J), 1/1000,
                    beta, Phi, diag(J), 1)

plot(sapply(hh$bChain, function(x) x[[1]]))
plot(sapply(h1$bChain, function(x) x[[1]]))


h2 <- msglmmHurdle(1e4, Y[,1], rep(1,n), X[1:2], M[1:2], offset[,1], 0.1, 0.1, 
                   0.1, 0.1, 1e4, beta[,1], beta[,2], Phi[,1], Phi[,2], 
                   diag(2), 1)



