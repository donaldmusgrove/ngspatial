library(batchmeans)
library(ngspatial)
library(RcppArmadillo)

setwd("C:/Users/donnie/Dropbox/github")
setwd("C:/Dropbox/github")


################################################################################
### Build and test the ngspatial package
################################################################################
setwd("C:/Dropbox/github/")
devtools::document("ngspatial")
devtools::build("ngspatial")
#devtools::build("ngspatial", binary=TRUE)


detach("package:ngspatial", unload=TRUE)
install.packages("ngspatial_1.1.tar.gz", repos = NULL, type="source")
#install.packages("ngspatial_1.1.tar.gz", repos = NULL, type="source")

library(ngspatial)



fbinom <- sparse.msglmm(family="binomial", Y, X, M, offset, Ntrials, 
                    beta.start, Psi.start, Sigma.start, Delta0=NULL,
                    tol, minit=1000, maxit=1e4, 
                    tune, hyper, verbose=1)


fpois <- sparse.msglmm(family="poisson", Y=Ycount, X=X, M=M, offset=offset, Ntrials=NULL, 
                       beta0=beta.start, Psi0=Psi.start, Sigma0=Sigma.start,  Delta0=NULL,
                       tol=tol, minit=1e3, maxit=2e4, 
                       tune=tune, hyper=hyper, verbose=1)



tuneH  <- c(tune, sigma.h=list(c(0.00001,0.00001,0.00001)))
hyperH <- c(hyper, tau.h = 1e-5)

fHpois <- sparse.msglmm(family="poisson", Y=Ycount, X=X, M=M, offset=offset, Ntrials=NULL, 
                       beta0=beta.start, Psi0=Psi.start, Sigma0=Sigma.start,  Delta0=matrix(rnorm(q*J),q,J),
                       tol=tol, minit=100, maxit=1e3, 
                       tune=tuneH , hyper=hyperH, verbose=1)



tuneH <- tune; tuneH$sigma.b <- c(0.001,0.001); tuneH$sigma.p <- c(0.1,0.1)


fHurdle <- sparse.msglmm(family="hurdle", Y=Ycount[,1], X=X[1:2], M=M[1:2], offset=offset[,1:2], Ntrials=NULL, 
                        beta0=beta.start[1:2], Psi0=Psi.start[1:2], Sigma0=Sigma.start[1:2,1:2],  Delta0=NULL,
                        tol=tol, minit=100, maxit=1e4, 
                        tune=tuneH, hyper=hyper, verbose=1)

fHurdle$beta.mcse
fHurdle$beta.accept
fHurdle$Rho.est
fHurdle$coef
fHurdle$dicZ



hyperG <- c(hyper, tau2e.a = list(1/2), tau2e.b = list(1e6))
tau2e.start <- rep(1,J)



ff <- sparse.msglmm.fit.gaussian(Ygauss, X, M, 
  beta.start, Psi.start, Sigma.start, tau2e.start, 
  tol, minit=100, maxit=1e4, 
  hyperG, verbose)



################################################################################
### Test the R function that calls the cpp functions
################################################################################
source("_sparseTemp.R")

### Make fake data
n   <- 400
m   <- sqrt(n)
J   <- 3
q   <- 50
rho <- 0.8
A   <- adjacency.matrix(sqrt(n))
X   <- cbind(x = rep(0:(m - 1) / (m - 1), times = m), rep(0:(m - 1) / (m - 1), each = m))
  colnames(X) = c("x1", "x2")
P   <- diag(n) - X%*%solve(t(X)%*%X)%*%t(X)
M0  <- P%*%A%*%P
M   <- eigen(M0)$vectors[,1:q]


Qsim   <- t(M) %*% (diag(rowSums(A), n) - A) %*% M
Lsim   <- t(chol(solve(Qsim)))
Mssim  <- M%*%Lsim
Rhosim <- (1-rho)*diag(J) + rho*matrix(1,J,J)
Ssim   <- diag(sqrt(c(4:6)))
Sigsim <- Ssim%*%Rhosim%*%Ssim
P      <- sapply(1:q, function(xx) t(chol(Sigsim)) %*%rnorm(J))

etasim  <- sapply(1:J, function(j) X%*%c(1,0.4) + Mssim%*%P[j,])
probsim <- plogis(etasim)
lambsim <- exp(etasim)

Y       <- apply(probsim, 2, function(x)  rbinom(n, 1,x))
Ycount  <- apply(lambsim, 2, function(x)  rpois(n,x))
Ygauss  <- etasim + matrix(rnorm(n*J,0,1),ncol=J)

X       <- lapply(1:J, function(xx) X)
M       <- lapply(1:J, function(xx) Mssim)
Ntrials <- matrix(1,n,J)
offset  <- matrix(0,n,J)
beta.start  <- list(c(1,0.4),c(1,0.4),c(1,0.4))
Psi.start   <- t(P)
Sigma.start <- Sigsim
tol         <- 0.01
minit       <- maxit <- 1e4
verbose     <- 1

tune <- list(sigma.b = rep(0.05,J),
             sigma.p = rep(0.05,J))
hyper <- list(tau.b = 1/1e6,
              nu    = 2,
              tau.s = 10^-5)


Times <- numeric(5)

### Run multivariate binomial model
start <- Sys.time()
ff <- sparse.msglmm.fit.binomial(Y, X, Ntrials, offset, M, 
                                 beta.start, Psi.start, Sigma.start, 
                                 tol, minit=minit, maxit=maxit, 
                                 tune, hyper, verbose)
Times[1] <- Sys.time()-start


### Run multivariate poisson model
start <- Sys.time()
ff <- sparse.msglmm.fit.poisson(Ycount, X, offset, M, 
                                beta.start, Psi.start, Sigma.start,  
                                tol, minit=minit, maxit=maxit, 
                                tune, hyper, verbose)
Times[2] <- Sys.time()-start


### Run multivariate Poisson model
Delta.start <- matrix(0,q,J)
tuneH  <- c(tune, sigma.h=list(c(0.01,0.01,0.01)))
hyperH <- c(hyper, tau.h = 1e6)
start <- Sys.time()
ff <- sparse.msglmm.fit.poissonH(Ycount, X, offset, M, 
                                beta.start, Psi.start, Sigma.start, Delta.start, 
                                tol, minit=minit, maxit=maxit, 
                                tuneH, hyperH, verbose)
Times[3] <- Sys.time()-start


### Run multivariate hurdle model
Ybin <- ifelse(Ycount[,1]==0,0,1)
YH <- cbind(Ybin, Ycount[,1])
XH <- X[1:2]
offsetH <- offset[,1:2]
MH <- M[1:2]
start <- Sys.time()
ff <- sparse.msglmm.fit.hurdle(YH, XH, offsetH, MH, 
                                beta.start[1:2], Psi.start[,1:2], Sigma.start[1:2,1:2],  
                                tol, minit=minit, maxit=maxit, 
                                tune, hyper, verbose)
Times[4] <- Sys.time()-start


ff$Sigma.est
ff$Rho.est



### Run multivariate Gaussian model
hyperG <- c(hyper, tau2e.a = list(1/2), tau2e.b = list(1e6))
tau2e.start <- rep(1,J)

start <- Sys.time()
ff <- sparse.msglmm.fit.gaussian(Ygauss, X, M, 
                                 beta.start, Psi.start, Sigma.start, tau2e.start, 
                                 tol, minit=100, maxit=1e4, 
                                 hyperG, verbose)
Times[5] <- Sys.time()-start

names(ff)
ff$coef
ff$Psi.mcse
ff$Sigma.mcse
ff$tau2e.est
ff$Rho.est
ff$Sigma.est
ff$dic

hist(ff$beta.sample[[3]][,1])



beta0 <- lm(Y[,j] ~ X[[j]] -1 + M[[j]])$coef[1:2]
psi0 <- lm(Y[,j] ~ X[[j]] -1 + M[[j]])$coef[-c(1:2)]

ttt <- replicate(1000,tau2eSgauss(Y[,j], X[[j]], M[[j]], beta0, psi0, 0.5, 1e6))

e1 <- sum((Y[,j] - X[[j]]%*%beta0 - M[[j]]%*%psi0)^2)
hist(sqrt(1/rgamma(1e4,0.5*n+0.5, 0.5*e1+1/1e6)))





###
ff$coefficients
ff$Rho.est
diag(ff$Sigma.est)

ff$beta.accept
ff$Psi.accept
unlist(ff$beta.mcse)
max(unlist(ff$Psi.mcse))

PP <- rnorm(q*J,1:q,1)
PP <- rep(PP,5)

array(PP,c(q,J,5))


################################################################################
### Test the individual cpp functions
################################################################################
setwd("C:/Dropbox/github/ngspatial/temp")

Rcpp::sourceCpp('msglmmBinomial.cpp')
Rcpp::sourceCpp('msglmmPoisson.cpp')
Rcpp::sourceCpp('msglmmPoissonHetero.cpp')
Rcpp::sourceCpp('msglmmHurdle.cpp')
Rcpp::sourceCpp('msglmmGaussian.cpp')
Rcpp::sourceCpp('bm.cpp')


  running.iter = 1
  five.pct     = round(0.05*maxit, 0)
  five.pct     = round(seq(five.pct, maxit, five.pct), 0)


### Run multivariate Gaussian model
tuneG  <- tune
hyperG <- c(hyper, taue.a = list(1), taue.b = list(1))
taue.start <- rep(1,J)

start <- Sys.time()

ff <- msglmmGaussian(iterations  = 100, 
                     Y           = Ycount, 
                     X           = X,  
                     M           = M, 
                     sigmab      = tuneG$sigma.b, 
                     sigmap      = tuneG$sigma.p, 
                     taub        = hyperG$tau.b, 
                     nu          = hyperG$nu, 
                     taus        = hyperG$tau.s,  
                     tauea       = hyperG$taue.a,
                     taueb       = hyperG$taue.b,
                     beta0       = beta.start, 
                     Psi0        = Psi.start, 
                     taue0       = taue.start,
                     Sigma0      = Sigma.start, 
                     verbose     = verbose, 
                     runningIter = running.iter, 
                     fivepct     = five.pct, 
                     maxit       = maxit)

Sys.time()-start
rm(ff)
gc()

bmmat(t(ff$tau2eChain))[,2]


### Run multivariate poisson model with hetero rfx
tuneH  <- c(tune, sigma.h=list(c(0.05,0.05,0.05)))
hyperH <- c(hyper, tau.h = list(1e-6))
ff <- msglmmPoissonHetero(iterations  = minit, 
                          Y           = Ycount, 
                          X           = X,  
                          M           = M, 
                          offset      = offset, 
                          sigmab      = tuneH$sigma.b, 
                          sigmap      = tuneH$sigma.p, 
                          sigmah      = tuneH$sigma.h, 
                          taub        = hyperH$tau.b, 
                          nu          = hyperH$nu, 
                          taus        = hyperH$tau.s,  
                          tauh        = hyperH$tau.h,
                          beta0       = beta.start, 
                          Psi0        = Psi.start, 
                          Sigma0      = Sigma.start, 
                          Delta0      = matrix(rnorm(q*J)*0.1,ncol=J),
                          verbose     = verbose, 
                          runningIter = running.iter, 
                          fivepct     = five.pct, 
                          maxit       = maxit)


### Run hurdle/zip model
YH <- cbind(Y[,1], Ycount[,1])
XH <- X[1:2]
MH <- M[1:2]
offsetH <- offset[,1:2]
hh <- msglmmHurdle(iterations  = iterations, 
                    Y           = YH, 
                    X           = XH, 
                    M           = MH, 
                    offset      = offset[,1:2], 
                    sigmab      = tune$sigma.b[1:2],
                    sigmap      = tune$sigma.p[1:2],
                    taub        = hyper$tau.b,
                    nu          = hyper$nu,
                    taus        = hyper$tau.s,
                    beta0       = beta.start[1:2],
                    Psi0        = Psi.start[,1:2],
                    Sigma0      = Sigma.start[1:2,1:2],
                    verbose     = verbose,
                    runningIter = running.iter,
                    fivepct     = five.pct,
                    maxit       = maxit)

SigInv  = solve(Sigma.start[1:2,1:2]);
SigI    = kronecker(SigInv, diag(q));

psi1Shurd(Psi.start[,1:2], tune$sigma.p[1], YY[,1], X[[1]], M[[1]], offset[,1], 
          beta.start[[1]], SigI)













### Reformat beta.chain and store estimates
beta.chain   <- lapply(1:J, function(j) sapply(1:running.iter, function(i) beta[,j][[i]]))
betanames    <- lapply(X, colnames)
coefficients <- beta.mcse <- vector("list", length=J)
for(j in 1:J){
  if(!is.null(nrow(beta.chain[[j]]))){
    temp              <- bmmat(t(beta.chain[[j]]))
    coefficients[[j]] <- temp[,1]
    beta.mcse[[j]]    <- temp[,2]
  } else{
    temp              <- bm(beta.chain[[j]]) 
    coefficients[[j]] <- temp$est
    beta.mcse[[j]]    <- temp$se
  }
  names(coefficients[[j]]) <- betanames[[j]]
  names(beta.mcse[[j]]) <- betanames[[j]]
}


### Store estimates of Psi
Psi.est = Psi.se = matrix(, q,J)
for(j in 1:J){
  temp         = bmmat(t(Psi[,j,]))
  Psi.est[,j] = temp[,1]
  Psi.se[,j]  = temp[,2]
}

### Store estimates of Sigma
Sigma.est <- apply(Sigma,c(1,2),function(x) unlist(bm(x)))[1,,]
Sigma.se  <- apply(Sigma,c(1,2),function(x) unlist(bm(x)))[2,,]


### Store estimates of Correlation matrix (cov2cor(Sigma))
Rho     <- acov2corcpp(Sigma)
Rho.est <- apply(Rho,c(1,2),function(x) unlist(bm(x)))[1,,]
Rho.se  <- apply(Rho,c(1,2),function(x) unlist(bm(x)))[2,,]
diag(Rho.se) <- NA

### Compute acceptance rates
beta.accept <- Psi.accept <- numeric(J)

  for(j in 1:J){
    if(!is.null(nrow(beta.chain[[j]]))){
      beta.accept[j] = sum(diff(beta.chain[[j]][1,]) != 0) / running.iter
    } else{
      beta.accept[j] = sum(diff(beta.chain[[j]]) != 0) / running.iter
    }

    Psi.accept[j] = sum(diff(Psi[1,j,]) != 0) / running.iter
  }


### Compute linear predictors and fitted values
linear.predictors = matrix(0,n,J)
fitted.values     = matrix(0,n,J)

j = 1 
i = 1

for(j in 1:J){
  for(i in 1:running.iter){
    eta = offset[,j] + X[[j]]%*%beta[[i]][[j]] + M[[j]]%*%Psi[,j,i]
    linear.predictors[,j] = linear.predictors[,j] + eta/running.iter
    prob = pnorm(eta)
    fitted.values[,j] = fitted.values[,j] + prob/running.iter
  }
}

Y-fitted.values





### Check convergence of each set of variables
sapply(beta.chain, function(x) bmmat(t(x)))
betaSE     <- sapply(beta.chain, function(x) bmmat(t(x))[,2])
any(betaSE>tol)

PsiSE   <- apply(Psi,c(1,2),function(x) unlist(bm(x)))[2,,]
any(PsiSE>tol)

SigmaSE <- apply(Sigma,c(1,2),function(x) unlist(bm(x)))[2,,]
any(SigmaSE>tol)

#####



h1 <- msglmmPoisson(1e4, Y, X, M, offset, 
                    rep(0.1,J), rep(0.1,J), 1/1000,
                    beta, Phi, diag(J), 1)

plot(sapply(hh$bChain, function(x) x[[1]]))
plot(sapply(h1$bChain, function(x) x[[1]]))


h2 <- msglmmHurdle(1e4, Y[,1], rep(1,n), X[1:2], M[1:2], offset[,1], 0.1, 0.1, 
                   0.1, 0.1, 1e4, beta[,1], beta[,2], Phi[,1], Phi[,2], 
                   diag(2), 1)






