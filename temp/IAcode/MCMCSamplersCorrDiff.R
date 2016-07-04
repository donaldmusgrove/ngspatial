
Rcpp::sourceCpp('code/bm.cpp')
Rcpp::sourceCpp('code/MCMCSamplersCorrDiff.cpp')

HPD <- function(obj, prob = 0.95){
  obj   <- as.matrix(obj)
  vals  <- apply(obj, 2, sort)
  nsamp <- nrow(vals)
  npar  <- ncol(vals)
  gap   <- max(1, min(nsamp - 1, round(nsamp * prob)))
  init  <- 1:(nsamp - gap)
  inds  <- apply(vals[init + gap, ,drop=FALSE] - vals[init, ,drop=FALSE],
                2, which.min)
  ans   <- cbind(vals[cbind(inds, 1:npar)],
                 vals[cbind(inds + gap, 1:npar)])
  dimnames(ans) <- list(colnames(obj), c("lower", "upper"))
  return(ans)
}

################################################################################
zipMCMC <- function(iterations, Y, Z, X1, X2, M1, M2, logE, Qs1, Qs2, tunea=0.1, tuneb=0.1, 
                    tunep1=0.01, tunep2=0.01, tunes=c(0.1,0.1, 0.1)){

  ### Setup a few constants
  L1  <- t(solve(chol(Qs1)))
  Ms1 <- M1%*%L1
  L2  <- t(solve(chol(Qs2)))
  Ms2 <- M2%*%L2
  
  ### Draw from the posterior
  PP <- MCMCzip(iterations, Y, Z, X1, X2, Ms1, Ms2, logE, tunea, tuneb, tunep1, tunep2,
                tunes)

  ### Extract the chains
  alpha.chain <- PP[[1]]
  beta.chain  <- PP[[2]]
  phi1.chain  <- PP[[3]]
  phi2.chain  <- PP[[4]]
  s.chain     <- PP[[5]]
  v           <- PP[[6]]

  ### Calculate acceptance rates
  phi1.acc   <- sum(diff(phi1.chain[1,]) != 0)/iterations
  phi2.acc   <- sum(diff(phi2.chain[1,]) != 0)/iterations
  alpha.acc  <- sum(diff(alpha.chain[1,]) != 0)/iterations
  beta.acc   <- sum(diff(beta.chain[1,]) != 0)/iterations
  s1.acc     <- sum(diff(s.chain[1,]) != 0)/iterations
  s2.acc     <- sum(diff(s.chain[2,]) != 0)/iterations
  rho.acc    <- sum(diff(s.chain[3,]) != 0)/iterations
  accept     <- c(phi1  = phi1.acc,
                  phi2  = phi2.acc, 
                  alpha = alpha.acc, 
                  beta  = beta.acc,
                  s1    = s1.acc,
                  s2    = s2.acc,
                  rho   = rho.acc)


  ### Extract elements of Sigma 
  s1.chain    <- s.chain[1,]
  s2.chain    <- s.chain[2,]
  rho.chain   <- s.chain[3,]


  ### Calculate batch mean estimates
  alphahat <- bmmatcpp(t(alpha.chain))
  betahat  <- bmmatcpp(t(beta.chain))
  phi1hat  <- bmmatcpp(t(phi1.chain))
  phi2hat  <- bmmatcpp(t(phi2.chain))
  s1hat    <- bmcpp(s1.chain)
  s2hat    <- bmcpp(s2.chain)
  rhohat   <- bmcpp(rho.chain)
  
  medians  <- list(s1    = median(s1.chain),
                   s2    = median(s2.chain),
                   rho   = median(rho.chain))


  est      <- list(alpha = alphahat[,1],
                   beta  = betahat[,1],
                   phi1  = phi1hat[,1],
                   phi2  = phi2hat[,1],
                   s1    = s1hat[1],
                   s2    = s2hat[1],
                   rho   = rhohat[1])

  MCSE     <- list(alpha = alphahat[,2],
                   beta  = betahat[,2],
                   phi1  = phi1hat[,2],
                   phi2  = phi2hat[,2],
                   s1    = s1hat[2],
                   s2    = s2hat[2],
                   rho   = rhohat[2])

  ### Calculate HPDs
  alphaHPD <- apply(alpha.chain, 1, HPD)
  betaHPD  <- apply(beta.chain, 1, HPD)
  s1HPD    <- HPD(s1.chain)
  s2HPD    <- HPD(s2.chain)
  rhoHPD   <- HPD(rho.chain)
  
  HPDs     <- list(alpha = alphaHPD,
                   beta  = betaHPD,
                   s1    = s1HPD,
                   s2    = s2HPD,
                   rho   = rhoHPD)

  ### Calculate DIC
  D.bar <- mean(v)
  pD    <- D.bar + 2 * dhurdle(Y, Z, X1, X2, Ms1, Ms2, logE, alphahat[,1], betahat[,1], phi1hat[,1], phi2hat[,1])
  dic   <- D.bar + pD

  ### Calculate fitted values
  eta1 <- X1%*%alphahat[,1] + Ms1%*%phi1hat[,1]
  eta2 <- X2%*%betahat[,1]  + Ms2%*%phi2hat[,1]

  pHurdle <- plogis(eta1)
  lambda  <- exp(logE + eta2)
  pZip    <- pHurdle * (1 - exp(-lambda))
  
  fitted  <- list(eta1    = eta1,
                  eta2    = eta2,
                  pHurdle = pHurdle,
                  pZip    = pZip,
                  lambda  = lambda)
  
  ### Return a list
  list(est     = est,
       medians = medians,
       MCSE    = MCSE,
       HPDs    = HPDs,
       dic     = dic,
       pD      = pD,
       accept  = accept,
       fitted  = fitted)
}



