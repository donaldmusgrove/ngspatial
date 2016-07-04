
sparse.msglmm.fit.binomial = function(Y, X, Ntrials, offset, M, 
                                      beta.start, Psi.start, Sigma.start, 
                                      tol, minit, maxit, 
                                      tune, hyper, verbose)
{

  iterations = minit
  n          = nrow(Y)
  J          = ncol(Y)
  q          = ncol(M[[1]])
  pJ         = sapply(X,ncol) 
    
  beta       = matrix(beta.start,1,J)
  Psi        = array(Psi.start, c(q,J,1))
  Sigma      = array(Sigma.start, c(J,J,1))
  v          = c()   
    
  running.iter = 1
  five.pct     = round(0.05*maxit, 0)
  five.pct     = round(seq(five.pct, maxit, five.pct), 0)

  if(verbose){
    cat("Now sampling. Caution, MCMC can be time consuming. \n")
  }
  
  
  repeat{
    #### MCMC sampling
    chains = msglmmBinomial(iterations  = iterations, 
                            Y           = Y, 
                            Ntrials     = Ntrials, 
                            X           = X, 
                            M           = M, 
                            offset      = offset, 
                            sigmab      = tune$sigma.b,
                            sigmap      = tune$sigma.p,
                            taub        = hyper$tau.b,
                            nu          = hyper$nu,
                            taus        = hyper$tau.s,
                            beta0       = beta[running.iter,],
                            Psi0        = Psi[,,running.iter],
                            Sigma0      = Sigma[,,running.iter],
                            verbose     = verbose,
                            runningIter = running.iter,
                            fivepct     = five.pct,
                            maxit       = maxit)
                            
    running.iter = chains$runningIter

    ### Combine latest chains with previous chains
    beta  = rbind(beta, chains$bChain[-1,])
    Psi   = arraybind(Psi, chains$PsiChain[,,-1])
    Sigma = arraybind(Sigma, chains$SigmaChain[,,-1])
    v     = c(v,chains$v)
        
    if (running.iter >= maxit){
      break
    }
    
    done = TRUE
            
    ### Check if spatial effects converged
    PsiSE     = apply(Psi,c(1,2),function(x) unlist(bm(x)))[2,,]
    if (any(PsiSE>tol)){
      done = FALSE
    }
    
    ### Check if covariance matrix converged
    if (done){
      SigmaSE = apply(Sigma,c(1,2),function(x) unlist(bm(x)))[2,,]
      if (any(SigmaSE>tol)){
        done = FALSE
      }
    }
        
    ### Check if regression coefficients converged
    if (done){
      beta.chain = lapply(1:J, function(j) matrix(unlist(beta[,j]), ncol=pJ[j], byrow=TRUE))
      betaSE     = sapply(beta.chain, function(x) bmmat(x)[,2])
      if (any(betaSE>tol)){
        done = FALSE
      }
    }  
        
    ### If all converged, break the repeat
    if (done){
      break
    }      
        
    ### If at least one did not converge and maxit has not been reached, 
    ### continue sampling until maxit is reached
    else{
      iterAdd    = running.iter + minit
      iterations = ifelse(iterAdd>maxit, maxit-running.iter+1, iterations)
    }
  }
  
  if(verbose){
    cat("Sampling complete. Now compiling results. \n")
  }
  
    
  ### Extract regression coefficients
  if(verbose){
    cat("Now extracting regression coefficients. \n")
  }
  
  beta.chain = lapply(1:J, function(j) matrix(unlist(beta[,j]), ncol=pJ[j], byrow=TRUE))
    
  betanames    = lapply(X, colnames)
  coefficients = beta.mcse = vector("list", length=J)
  for(j in 1:J){
  if(!is.null(ncol(beta.chain[[j]]))){
    temp              = bmmat(beta.chain[[j]])
    coefficients[[j]] = temp[,1]
    beta.mcse[[j]]    = temp[,2]
  } else{
    temp              = bm(beta.chain[[j]]) 
    coefficients[[j]] = temp$est
    beta.mcse[[j]]    = temp$se
  }
    names(coefficients[[j]]) = betanames[[j]]
    names(beta.mcse[[j]])    = betanames[[j]]
  }
    
    
  ### Extract spatial effects
  if(verbose){
    cat("=> Extracting spatial random effects. \n")
  }
  Psi.est = Psi.se = matrix(, q,J)
  for(j in 1:J){
    temp        = bmmat(t(Psi[,j,]))
    Psi.est[,j] = temp[,1]
    Psi.se[,j]  = temp[,2]
  }
    
  ### Extract covariance matrix estimate
  if(verbose){
    cat("=> Extracting covariance matrix. \n")
  }
  Sigma.est = apply(Sigma,c(1,2),function(x) unlist(bm(x)))[1,,]
  Sigma.se  = apply(Sigma,c(1,2),function(x) unlist(bm(x)))[2,,]
  
  ### Create correlation matrix and extract estimates
  if(verbose){
    cat("=> Computing and extracting correlation matrix. \n")
  }
  Rho     = acov2corcpp(Sigma)
  Rho.est = apply(Rho,c(1,2),function(x) unlist(bm(x)))[1,,]
  Rho.se  = apply(Rho,c(1,2),function(x) unlist(bm(x)))[2,,]
  diag(Rho.se) = NA


  ### Acceptance rates
  if(verbose){
    cat("=> Calculating parameter acceptance rates. \n")
  }
  beta.accept = Psi.accept = numeric(J)
  for(j in 1:J){
    if(!is.null(ncol(beta.chain[[j]]))){
      beta.accept[j] = sum(diff(beta.chain[[j]][,1]) != 0) / running.iter
    } else{
      beta.accept[j] = sum(diff(beta.chain[[j]]) != 0) / running.iter
    }
      
    Psi.accept[j] = sum(diff(Psi[1,j,]) != 0) / running.iter
  }

        
  ### Fitted values and other things
  if(verbose){
    cat("=> Calculating linear predictors and fitted values. \n")
  }
  start <- Sys.time()
  
  linear.predictors = matrix(0,n,J)
  fitted.values     = matrix(0,n,J)
    
  for(j in 1:J){
    eta  = offset[,j] + X[[j]]%*%t(beta.chain[[j]]) + M[[j]]%*%Psi[,j,]
    prob = plogis(eta)
    linear.predictors[,j] = rowMeans(eta)  
    fitted.values[,j]    = rowMeans(prob)
  }
  fftime <- Sys.time()-start  
  if(verbose){
    fftime <- round(fftime,3)
    cat(paste0("=> Calculation of linear predictors and fitted values took ",fftime," seconds. \n"))
  } 
    
  ### Calculate residuals
  residuals = Y - fitted.values
    
  ### Calculate DIC
  D.bar = mean(v)
  pD = D.bar + 2 * dBinomMulti(Y, Ntrials, X, M, offset, coefficients, Psi.est)
  dic = D.bar + pD
    
  object = list(coefficients      = coefficients, 
                fitted.values     = fitted.values,
                linear.predictors = linear.predictors, 
                residuals         = residuals,
                beta.sample       = beta.chain, 
                Psi.sample        = Psi, 
                Sigma.sample      = Sigma,
                beta.mcse         = beta.mcse, 
                Psi.mcse          = Psi.se, 
                Sigma.mcse        = Sigma.se,
                Psi.est           = Psi.est, 
                Sigma.est         = Sigma.est,
                Rho.est           = Rho.est,
                Rho.se            = Rho.se,
                iter              = running.iter-1, 
                dic               = dic,
                D.bar             = D.bar, 
                pD                = pD, 
                beta.accept       = beta.accept, 
                Psi.accept        = Psi.accept)
  class(object) = c("sparse.msglmm")
  object
}



### Multivariate non-hetero poisson
sparse.msglmm.fit.poisson = function(Y, X, offset, M, 
                                     beta.start, Psi.start, Sigma.start, 
                                     tol, minit, maxit, 
                                     tune, hyper, verbose)
{


  iterations = minit
  n          = nrow(Y)
  J          = ncol(Y)
  q          = ncol(M[[1]])
  pJ         = sapply(X,ncol) 
    
  beta       = matrix(beta.start,1,J)
  Psi        = array(Psi.start, c(q,J,1))
  Sigma      = array(Sigma.start, c(J,J,1))
  v          = c()   
  
  running.iter = 1
  five.pct     = round(0.05*maxit, 0)
  five.pct     = round(seq(five.pct, maxit, five.pct), 0)

  if(verbose){
    cat("Now sampling. Caution, MCMC can be time consuming. \n")
  }
  
  
  repeat{
    #### MCMC sampling
    chains = msglmmPoisson(iterations  = iterations, 
                           Y           = Y, 
                           X           = X, 
                           M           = M, 
                           offset      = offset, 
                           sigmab      = tune$sigma.b,
                           sigmap      = tune$sigma.p,
                           taub        = hyper$tau.b,
                           nu          = hyper$nu,
                           taus        = hyper$tau.s,
                           beta0       = beta[running.iter,],
                           Psi0        = Psi[,,running.iter],
                           Sigma0      = Sigma[,,running.iter],
                           verbose     = verbose,
                           runningIter = running.iter,
                           fivepct     = five.pct,
                           maxit       = maxit)
                            
    running.iter = chains$runningIter

    ### Combine latest chains with previous chains
    beta  = rbind(beta, chains$bChain[-1,])
    Psi   = arraybind(Psi, chains$PsiChain[,,-1])
    Sigma = arraybind(Sigma, chains$SigmaChain[,,-1])
    v     = c(v,chains$v)
        
    if (running.iter >= maxit){
      break
    }
    
    done = TRUE
            
    ### Check if spatial effects converged
    PsiSE     = apply(Psi,c(1,2),function(x) unlist(bm(x)))[2,,]
    if (any(PsiSE>tol)){
      done = FALSE
    }
    
    ### Check if covariance matrix converged
    if (done){
      SigmaSE = apply(Sigma,c(1,2),function(x) unlist(bm(x)))[2,,]
      if (any(SigmaSE>tol)){
        done = FALSE
      }
    }
        
    ### Check if regression coefficients converged
    if (done){
      beta.chain = lapply(1:J, function(j) matrix(unlist(beta[,j]), ncol=pJ[j], byrow=TRUE))
      betaSE     = sapply(beta.chain, function(x) bmmat(x)[,2])
      if (any(betaSE>tol)){
        done = FALSE
      }
    }  
        
    ### If all converged, break the repeat
    if (done){
      break
    }      
        
    ### If at least one did not converge and maxit has not been reached, 
    ### continue sampling until maxit is reached
    else{
      iterAdd    = running.iter + minit
      iterations = ifelse(iterAdd>maxit, maxit-running.iter+1, iterations)
    }
  }
  
  if(verbose){
    cat("Sampling complete. Now compiling results. \n")
  }
  
    
  ### Extract regression coefficients
  if(verbose){
    cat("Now extracting regression coefficients. \n")
  }
  
  beta.chain = lapply(1:J, function(j) matrix(unlist(beta[,j]), ncol=pJ[j], byrow=TRUE))
    
  betanames    = lapply(X, colnames)
  coefficients = beta.mcse = vector("list", length=J)
  for(j in 1:J){
  if(!is.null(ncol(beta.chain[[j]]))){
    temp              = bmmat(beta.chain[[j]])
    coefficients[[j]] = temp[,1]
    beta.mcse[[j]]    = temp[,2]
  } else{
    temp              = bm(beta.chain[[j]]) 
    coefficients[[j]] = temp$est
    beta.mcse[[j]]    = temp$se
  }
    names(coefficients[[j]]) = betanames[[j]]
    names(beta.mcse[[j]])    = betanames[[j]]
  }
    
    
  ### Extract spatial effects
  if(verbose){
    cat("=> Extracting spatial random effects. \n")
  }
  Psi.est = Psi.se = matrix(, q,J)
  for(j in 1:J){
    temp        = bmmat(t(Psi[,j,]))
    Psi.est[,j] = temp[,1]
    Psi.se[,j]  = temp[,2]
  }
    
  ### Extract covariance matrix estimate
  if(verbose){
    cat("=> Extracting covariance matrix. \n")
  }
  Sigma.est = apply(Sigma,c(1,2),function(x) unlist(bm(x)))[1,,]
  Sigma.se  = apply(Sigma,c(1,2),function(x) unlist(bm(x)))[2,,]
  
  ### Create correlation matrix and extract estimates
  if(verbose){
    cat("=> Computing and extracting correlation matrix. \n")
  }
  Rho     = acov2corcpp(Sigma)
  Rho.est = apply(Rho,c(1,2),function(x) unlist(bm(x)))[1,,]
  Rho.se  = apply(Rho,c(1,2),function(x) unlist(bm(x)))[2,,]
  diag(Rho.se) = NA


  ### Acceptance rates
  if(verbose){
    cat("=> Calculating parameter acceptance rates. \n")
  }
  beta.accept = Psi.accept = numeric(J)
  for(j in 1:J){
    if(!is.null(ncol(beta.chain[[j]]))){
      beta.accept[j] = sum(diff(beta.chain[[j]][,1]) != 0) / running.iter
    } else{
      beta.accept[j] = sum(diff(beta.chain[[j]]) != 0) / running.iter
    }
      
    Psi.accept[j] = sum(diff(Psi[1,j,]) != 0) / running.iter
  }

        
  ### Fitted values and other things
  if(verbose){
    cat("=> Calculating linear predictors and fitted values. \n")
  }
  start <- Sys.time()
  
  linear.predictors = matrix(0,n,J)
  fitted.values     = matrix(0,n,J)
    
  for(j in 1:J){
    eta                   = offset[,j] + X[[j]]%*%t(beta.chain[[j]]) + M[[j]]%*%Psi[,j,]
    lambda                = exp(eta)
    linear.predictors[,j] = rowMeans(eta)  
    fitted.values[,j]     = rowMeans(lambda)
  }
  fftime <- Sys.time()-start  
  if(verbose){
    fftime <- round(fftime,3)
    cat(paste0("=> Calculation of linear predictors and fitted values took ",fftime," seconds. \n"))
  } 
    
  ### Calculate residuals
  residuals = Y - fitted.values
    
  ### Calculate DIC
  D.bar = mean(v)
  pD    = D.bar + 2 * dPoissMulti(Y, X, M, offset, coefficients, Psi.est)
  dic   = D.bar + pD
    
  object = list(coefficients      = coefficients, 
                fitted.values     = fitted.values,
                linear.predictors = linear.predictors, 
                residuals         = residuals,
                beta.sample       = beta.chain, 
                Psi.sample        = Psi, 
                Sigma.sample      = Sigma,
                beta.mcse         = beta.mcse, 
                Psi.mcse          = Psi.se, 
                Sigma.mcse        = Sigma.se,
                Psi.est           = Psi.est, 
                Sigma.est         = Sigma.est,
                Rho.est           = Rho.est,
                Rho.se            = Rho.se,
                iter              = running.iter-1, 
                dic               = dic,
                D.bar             = D.bar, 
                pD                = pD, 
                beta.accept       = beta.accept, 
                Psi.accept        = Psi.accept)
  class(object) = c("sparse.msglmm")
  object
}





### Multivariate hetero
sparse.msglmm.fit.poissonH = function(Y, X, offset, M, 
                                     beta.start, Psi.start, Sigma.start, Delta.start,
                                     tol, minit, maxit, 
                                     tune, hyper, verbose)
{


  iterations = minit
  n          = nrow(Y)
  J          = ncol(Y)
  q          = ncol(M[[1]])
  pJ         = sapply(X,ncol) 
    
  beta       = matrix(beta.start,1,J)
  Psi        = array(Psi.start, c(q,J,1))
  Sigma      = array(Sigma.start, c(J,J,1))
  Delta      = array(Delta.start, c(q,J,1))
  v          = c()   
  
  running.iter = 1
  five.pct     = round(0.05*maxit, 0)
  five.pct     = round(seq(five.pct, maxit, five.pct), 0)

  if(verbose){
    cat("Now sampling. Caution, MCMC can be time consuming. \n")
  }
  
  
  repeat{
    #### MCMC sampling
    chains = msglmmPoissonHetero(iterations  = iterations, 
                                 Y           = Y, 
                                 X           = X, 
                                 M           = M, 
                                 offset      = offset, 
                                 sigmab      = tune$sigma.b,
                                 sigmap      = tune$sigma.p,
                                 sigmah      = tune$sigma.h,
                                 taub        = hyper$tau.b,
                                 nu          = hyper$nu,
                                 taus        = hyper$tau.s,
                                 tauh        = hyper$tau.h,
                                 beta0       = beta[running.iter,],
                                 Psi0        = Psi[,,running.iter],
                                 Sigma0      = Sigma[,,running.iter],
                                 Delta0      = Delta[,,running.iter],
                                 verbose     = verbose,
                                 runningIter = running.iter,
                                 fivepct     = five.pct,
                                 maxit       = maxit)
                            
    running.iter = chains$runningIter

    ### Combine latest chains with previous chains
    beta  = rbind(beta, chains$bChain[-1,])
    Psi   = arraybind(Psi, chains$PsiChain[,,-1])
    Sigma = arraybind(Sigma, chains$SigmaChain[,,-1])
    Delta = arraybind(Delta, chains$DeltaChain[,,-1])
    v     = c(v,chains$v)
        
    if (running.iter >= maxit){
      break
    }
    
    done = TRUE
            
    ### Check if spatial effects converged
    PsiSE     = apply(Psi,c(1,2),function(x) unlist(bm(x)))[2,,]
    if (any(PsiSE>tol)){
      done = FALSE
    }
    
    ### Check if covariance matrix converged
    if (done){
      SigmaSE = apply(Sigma,c(1,2),function(x) unlist(bm(x)))[2,,]
      if (any(SigmaSE>tol)){
        done = FALSE
      }
    }
      
    ### Check if hetero rfx converged
    if (done){
      DeltaSE     = apply(Delta,c(1,2),function(x) unlist(bm(x)))[2,,]
      if (any(DeltaSE>tol)){
        done = FALSE
      }
    }
    
    
    ### Check if regression coefficients converged
    if (done){
      beta.chain = lapply(1:J, function(j) matrix(unlist(beta[,j]), ncol=pJ[j], byrow=TRUE))
      betaSE     = sapply(beta.chain, function(x) bmmat(x)[,2])
      if (any(betaSE>tol)){
        done = FALSE
      }
    }  
        
    ### If all converged, break the repeat
    if (done){
      break
    }      
        
    ### If at least one did not converge and maxit has not been reached, 
    ### continue sampling until maxit is reached
    else{
      iterAdd    = running.iter + minit
      iterations = ifelse(iterAdd>maxit, maxit-running.iter+1, iterations)
    }
  }
  
  if(verbose){
    cat("Sampling complete. Now compiling results. \n")
  }
  
    
  ### Extract regression coefficients
  if(verbose){
    cat("Now extracting regression coefficients. \n")
  }
  
  beta.chain = lapply(1:J, function(j) matrix(unlist(beta[,j]), ncol=pJ[j], byrow=TRUE))
    
  betanames    = lapply(X, colnames)
  coefficients = beta.mcse = vector("list", length=J)
  for(j in 1:J){
  if(!is.null(ncol(beta.chain[[j]]))){
    temp              = bmmat(beta.chain[[j]])
    coefficients[[j]] = temp[,1]
    beta.mcse[[j]]    = temp[,2]
  } else{
    temp              = bm(beta.chain[[j]]) 
    coefficients[[j]] = temp$est
    beta.mcse[[j]]    = temp$se
  }
    names(coefficients[[j]]) = betanames[[j]]
    names(beta.mcse[[j]])    = betanames[[j]]
  }
    
    
  ### Extract spatial effects
  if(verbose){
    cat("=> Extracting spatial random effects. \n")
  }
  Psi.est = Psi.se = matrix(, q,J)
  for(j in 1:J){
    temp        = bmmat(t(Psi[,j,]))
    Psi.est[,j] = temp[,1]
    Psi.se[,j]  = temp[,2]
  }
  
  
  ### Extract hetero rfx
  if(verbose){
    cat("=> Extracting spatial random effects. \n")
  }
  Delta.est = Delta.se = matrix(,q,J)
  for(j in 1:J){
    temp          = bmmat(t(Delta[,j,]))
    Delta.est[,j] = temp[,1]
    Delta.se[,j]  = temp[,2]
  }
  
    
  ### Extract covariance matrix estimate
  if(verbose){
    cat("=> Extracting covariance matrix. \n")
  }
  Sigma.est = apply(Sigma,c(1,2),function(x) unlist(bm(x)))[1,,]
  Sigma.se  = apply(Sigma,c(1,2),function(x) unlist(bm(x)))[2,,]
  
  ### Create correlation matrix and extract estimates
  if(verbose){
    cat("=> Computing and extracting correlation matrix. \n")
  }
  Rho     = acov2corcpp(Sigma)
  Rho.est = apply(Rho,c(1,2),function(x) unlist(bm(x)))[1,,]
  Rho.se  = apply(Rho,c(1,2),function(x) unlist(bm(x)))[2,,]
  diag(Rho.se) = NA


  ### Acceptance rates
  if(verbose){
    cat("=> Calculating parameter acceptance rates. \n")
  }
  beta.accept = Psi.accept = Delta.accept = numeric(J)
  for(j in 1:J){
    if(!is.null(ncol(beta.chain[[j]]))){
      beta.accept[j] = sum(diff(beta.chain[[j]][,1]) != 0) / running.iter
    } else{
      beta.accept[j] = sum(diff(beta.chain[[j]]) != 0) / running.iter
    }
      
    Psi.accept[j]   = sum(diff(Psi[1,j,]) != 0) / running.iter
    Delta.accept[j] = sum(diff(Delta[1,j,]) != 0) / running.iter
  }

        
  ### Fitted values and other things
  if(verbose){
    cat("=> Calculating linear predictors and fitted values. \n")
  }
  start <- Sys.time()
  
  linear.predictors = matrix(0,n,J)
  fitted.values     = matrix(0,n,J)
    
  for(j in 1:J){
    eta                   = offset[,j] + X[[j]]%*%t(beta.chain[[j]]) + M[[j]]%*%Psi[,j,] + M[[j]]%*%Delta[,j,]
    lambda                = exp(eta)
    linear.predictors[,j] = rowMeans(eta)  
    fitted.values[,j]     = rowMeans(lambda)
  }
  fftime <- Sys.time()-start  
  if(verbose){
    fftime <- round(fftime,3)
    cat(paste0("=> Calculation of linear predictors and fitted values took ",fftime," seconds. \n"))
  } 
    
  ### Calculate residuals
  residuals = Y - fitted.values
    
  ### Calculate DIC
  D.bar = mean(v)
  pD    = D.bar + 2 * dPoissMultiH(Y, X, M, offset, coefficients, Psi.est, Delta.est)
  dic   = D.bar + pD
    
  object = list(coefficients      = coefficients, 
                fitted.values     = fitted.values,
                linear.predictors = linear.predictors, 
                residuals         = residuals,
                beta.sample       = beta.chain, 
                Psi.sample        = Psi, 
                Sigma.sample      = Sigma,
                Delta.sample      = Delta,
                beta.mcse         = beta.mcse, 
                Psi.mcse          = Psi.se, 
                Sigma.mcse        = Sigma.se,
                Psi.est           = Psi.est, 
                Sigma.est         = Sigma.est,
                Rho.est           = Rho.est,
                Rho.se            = Rho.se,
                iter              = running.iter-1, 
                dic               = dic,
                D.bar             = D.bar, 
                pD                = pD, 
                beta.accept       = beta.accept, 
                Psi.accept        = Psi.accept)
  class(object) = c("sparse.msglmm")
  object
}






### Hurdle/ZIP model
sparse.msglmm.fit.hurdle = function(Y, X, offset, M, 
                                    beta.start, Psi.start, Sigma.start, 
                                    tol, minit, maxit, 
                                    tune, hyper, verbose)
{


  iterations = minit
  n          = nrow(Y)
  J          = ncol(Y)
  q          = ncol(M[[1]])
  pJ         = sapply(X,ncol) 
    
  beta       = matrix(beta.start,1,J)
  Psi        = array(Psi.start, c(q,J,1))
  Sigma      = array(Sigma.start, c(J,J,1))
  vH         = c()   ### DIC calculator for hurdle model
  vZ         = c()   ### DIC calculator for ZIP model
    
  running.iter = 1
  five.pct     = round(0.05*maxit, 0)
  five.pct     = round(seq(five.pct, maxit, five.pct), 0)

  if(verbose){
    cat("Now sampling. Caution, MCMC can be time consuming. \n")
  }
  
  
  repeat{
    #### MCMC sampling
    chains = msglmmHurdle(iterations  = iterations, 
                          Y           = Y, 
                          X           = X, 
                          M           = M, 
                          offset      = offset, 
                          sigmab      = tune$sigma.b,
                          sigmap      = tune$sigma.p,
                          taub        = hyper$tau.b,
                          nu          = hyper$nu,
                          taus        = hyper$tau.s,
                          beta0       = beta[running.iter,],
                          Psi0        = Psi[,,running.iter],
                          Sigma0      = Sigma[,,running.iter],
                          verbose     = verbose,
                          runningIter = running.iter,
                          fivepct     = five.pct,
                          maxit       = maxit)
                            
    running.iter = chains$runningIter

    ### Combine latest chains with previous chains
    beta  = rbind(beta, chains$bChain[-1,])
    Psi   = arraybind(Psi, chains$PsiChain[,,-1])
    Sigma = arraybind(Sigma, chains$SigmaChain[,,-1])
    vH     = c(vH,chains$vHurdle)
    vZ     = c(vZ,chains$vZIP)
    
    if (running.iter >= maxit){
      break
    }
    
    done = TRUE
            
    ### Check if spatial effects converged
    PsiSE     = apply(Psi,c(1,2),function(x) unlist(bm(x)))[2,,]
    if (any(PsiSE>tol)){
      done = FALSE
    }
    
    ### Check if covariance matrix converged
    if (done){
      SigmaSE = apply(Sigma,c(1,2),function(x) unlist(bm(x)))[2,,]
      if (any(SigmaSE>tol)){
        done = FALSE
      }
    }
        
    ### Check if regression coefficients converged
    if (done){
      beta.chain = lapply(1:J, function(j) matrix(unlist(beta[,j]), ncol=pJ[j], byrow=TRUE))
      betaSE     = sapply(beta.chain, function(x) bmmat(x)[,2])
      if (any(betaSE>tol)){
        done = FALSE
      }
    }  
        
    ### If all converged, break the repeat
    if (done){
      break
    }      
        
    ### If at least one did not converge and maxit has not been reached, 
    ### continue sampling until maxit is reached
    else{
      iterAdd    = running.iter + minit
      iterations = ifelse(iterAdd>maxit, maxit-running.iter+1, iterations)
    }
  }
  
  if(verbose){
    cat("Sampling complete. Now compiling results. \n")
  }
  
    
  ### Extract regression coefficients
  if(verbose){
    cat("Now extracting regression coefficients. \n")
  }
  
  beta.chain = lapply(1:J, function(j) matrix(unlist(beta[,j]), ncol=pJ[j], byrow=TRUE))
    
  betanames    = lapply(X, colnames)
  coefficients = beta.mcse = vector("list", length=J)
  for(j in 1:J){
  if(!is.null(ncol(beta.chain[[j]]))){
    temp              = bmmat(beta.chain[[j]])
    coefficients[[j]] = temp[,1]
    beta.mcse[[j]]    = temp[,2]
  } else{
    temp              = bm(beta.chain[[j]]) 
    coefficients[[j]] = temp$est
    beta.mcse[[j]]    = temp$se
  }
    names(coefficients[[j]]) = betanames[[j]]
    names(beta.mcse[[j]])    = betanames[[j]]
  }
    
    
  ### Extract spatial effects
  if(verbose){
    cat("=> Extracting spatial random effects. \n")
  }
  Psi.est = Psi.se = matrix(, q,J)
  for(j in 1:J){
    temp        = bmmat(t(Psi[,j,]))
    Psi.est[,j] = temp[,1]
    Psi.se[,j]  = temp[,2]
  }
    
  ### Extract covariance matrix estimate
  if(verbose){
    cat("=> Extracting covariance matrix. \n")
  }
  Sigma.est = apply(Sigma,c(1,2),function(x) unlist(bm(x)))[1,,]
  Sigma.se  = apply(Sigma,c(1,2),function(x) unlist(bm(x)))[2,,]
  
  ### Create correlation matrix and extract estimates
  if(verbose){
    cat("=> Computing and extracting correlation matrix. \n")
  }
  Rho     = acov2corcpp(Sigma)
  Rho.est = apply(Rho,c(1,2),function(x) unlist(bm(x)))[1,,]
  Rho.se  = apply(Rho,c(1,2),function(x) unlist(bm(x)))[2,,]
  diag(Rho.se) = NA


  ### Acceptance rates
  if(verbose){
    cat("=> Calculating parameter acceptance rates. \n")
  }
  beta.accept = Psi.accept = numeric(J)
  for(j in 1:J){
    if(!is.null(ncol(beta.chain[[j]]))){
      beta.accept[j] = sum(diff(beta.chain[[j]][,1]) != 0) / running.iter
    } else{
      beta.accept[j] = sum(diff(beta.chain[[j]]) != 0) / running.iter
    }
      
    Psi.accept[j] = sum(diff(Psi[1,j,]) != 0) / running.iter
  }

        
  ### Fitted values and other things
  if(verbose){
    cat("=> Calculating linear predictors and fitted values. \n")
  }
  start <- Sys.time()
  
  linear.predictors = matrix(0,n,J)
  fitted.values     = matrix(0,n,J)
    
  for(j in 1:J){
    eta                   = offset[,j] + X[[j]]%*%t(beta.chain[[j]]) + M[[j]]%*%Psi[,j,]
    lambda                = exp(eta)
    linear.predictors[,j] = rowMeans(eta)  
    fitted.values[,j]     = rowMeans(lambda)
  }
  fftime <- Sys.time()-start  
  if(verbose){
    fftime <- round(fftime,3)
    cat(paste0("=> Calculation of linear predictors and fitted values took ",fftime," seconds. \n"))
  } 
    
  ### Calculate residuals
  residuals = Y - fitted.values
  
  ### Calculate DIC for hurdle model
  D.barH = mean(vH)
  pDH    = D.barH + 2 * dHurdle(Y, X, M, offset, coefficients, Psi.est)
  dicH   = D.barH + pDH
  
    
  ### Calculate DIC for ZIP model
  D.barZ = mean(vZ)
  pDZ    = D.barZ + 2 * dZIP(Y, X, M, offset, coefficients, Psi.est)
  dicZ   = D.barZ + pDZ
    
  object = list(coefficients      = coefficients, 
                fitted.values     = fitted.values,
                linear.predictors = linear.predictors, 
                residuals         = residuals,
                beta.sample       = beta.chain, 
                Psi.sample        = Psi, 
                Sigma.sample      = Sigma,
                beta.mcse         = beta.mcse, 
                Psi.mcse          = Psi.se, 
                Sigma.mcse        = Sigma.se,
                Psi.est           = Psi.est, 
                Sigma.est         = Sigma.est,
                Rho.est           = Rho.est,
                Rho.se            = Rho.se,
                iter              = running.iter-1, 
                dicH              = dicH,
                D.barH            = D.barH, 
                pDH               = pDH, 
                dicZ              = dicZ,
                D.barZ            = D.barZ, 
                pDZ               = pDZ, 
                beta.accept       = beta.accept, 
                Psi.accept        = Psi.accept)
  class(object) = c("sparse.msglmm")
  object
}


### Multivariate Gaussian
sparse.msglmm.fit.gaussian = function(Y, X, M, 
                                     beta.start, Psi.start, Sigma.start, tau2e.start, 
                                     tol, minit, maxit, 
                                     hyper, verbose)
{


  iterations = minit
  n          = nrow(Y)
  J          = ncol(Y)
  q          = ncol(M[[1]])
  pJ         = sapply(X,ncol) 
    
  beta       = matrix(beta.start,1,J)
  Psi        = matrix(c(Psi.start), q*J,1)
  tau2e      = matrix(tau2e.start, J, 1)
  Sigma      = array(Sigma.start, c(J,J,1))
  v          = c()   
  
  running.iter = 1
  five.pct     = round(0.05*maxit, 0)
  five.pct     = round(seq(five.pct, maxit, five.pct), 0)

  if(verbose){
    cat("Now sampling. Caution, MCMC can be time consuming. \n")
  }
  
  
  repeat{
    #### MCMC sampling
    chains = msglmmGaussian(iterations  = iterations, 
                            Y           = Y, 
                            X           = X, 
                            M           = M, 
                            taub        = hyper$tau.b,
                            nu          = hyper$nu,
                            taus        = hyper$tau.s,
                            tau2ea      = hyper$tau2e.a,
                            tau2eb      = hyper$tau2e.b,
                            beta0       = beta[running.iter,],
                            Psi0        = Psi[,running.iter],
                            tau2e0      = tau2e[,running.iter],
                            Sigma0      = Sigma[,,running.iter],
                            verbose     = verbose,
                            runningIter = running.iter,
                            fivepct     = five.pct,
                            maxit       = maxit)
                            
    running.iter = chains$runningIter

    ### Combine latest chains with previous chains
    beta  = rbind(beta, chains$bChain[-1,])
    Psi   = cbind(Psi, chains$PsiChain[,-1])
    tau2e = cbind(tau2e, chains$tau2eChain[,-1])
    Sigma = arraybind(Sigma, chains$SigmaChain[,,-1])
    v     = c(v,chains$v)
        
    if (running.iter >= maxit){
      break
    }
    
    done = TRUE
            
    ### Check if spatial effects converged
    PsiSE     = bmmat(t(Psi))[,2]
    if (any(PsiSE>tol)){
      done = FALSE
    }
    
    ### Check if covariance matrix converged
    if (done){
      SigmaSE = apply(Sigma,c(1,2),function(x) unlist(bm(x)))[2,,]
      if (any(SigmaSE>tol)){
        done = FALSE
      }
    }
    
    ### Check if error variance tau2e converged
    if (done){
      tau2eSE     = bmmat(t(tau2e))[,2]
      if (any(tau2eSE>tol)){
        done = FALSE
      }
    }
        
    ### Check if regression coefficients converged
    if (done){
      beta.chain = lapply(1:J, function(j) matrix(unlist(beta[,j]), ncol=pJ[j], byrow=TRUE))
      betaSE     = sapply(beta.chain, function(x) bmmat(x)[,2])
      if (any(betaSE>tol)){
        done = FALSE
      }
    }  
        
    ### If all converged, break the repeat
    if (done){
      break
    }      
        
    ### If at least one did not converge and maxit has not been reached, 
    ### continue sampling until maxit is reached
    else{
      iterAdd    = running.iter + minit
      iterations = ifelse(iterAdd>maxit, maxit-running.iter+1, iterations)
    }
  }
  
  if(verbose){
    cat("Sampling complete. Now compiling results. \n")
  }
  
    
  ### Extract regression coefficients
  if(verbose){
    cat("Now extracting regression coefficients. \n")
  }
  
  beta.chain = lapply(1:J, function(j) matrix(unlist(beta[,j]), ncol=pJ[j], byrow=TRUE))
    
  betanames    = lapply(X, colnames)
  coefficients = beta.mcse = vector("list", length=J)
  for(j in 1:J){
  if(!is.null(ncol(beta.chain[[j]]))){
    temp              = bmmat(beta.chain[[j]])
    coefficients[[j]] = temp[,1]
    beta.mcse[[j]]    = temp[,2]
  } else{
    temp              = bm(beta.chain[[j]]) 
    coefficients[[j]] = temp$est
    beta.mcse[[j]]    = temp$se
  }
    names(coefficients[[j]]) = betanames[[j]]
    names(beta.mcse[[j]])    = betanames[[j]]
  }
    
    
  ### Extract spatial effects
  if(verbose){
    cat("=> Extracting spatial random effects. \n")
  }
  
  ### Put spatial effects into array form
  iterP   = ncol(Psi)
  Psi     = array(Psi,c(q,J,iterP))
  Psi.est = Psi.se = matrix(, q,J)
  for(j in 1:J){
    temp        = bmmat(t(Psi[,j,]))
    Psi.est[,j] = temp[,1]
    Psi.se[,j]  = temp[,2]
  }
  
  
  ### Extract error variance
  if(verbose){
    cat("=> Extracting error variance. \n")
  }
  temp      = bmmat(t(tau2e))
  tau2e.est = temp[,1]
  tau2e.se  = temp[,2]

  
  
  ### Extract covariance matrix estimate
  if(verbose){
    cat("=> Extracting covariance matrix. \n")
  }
  Sigma.est = apply(Sigma,c(1,2),function(x) unlist(bm(x)))[1,,]
  Sigma.se  = apply(Sigma,c(1,2),function(x) unlist(bm(x)))[2,,]
  
  ### Create correlation matrix and extract estimates
  if(verbose){
    cat("=> Computing and extracting correlation matrix. \n")
  }
  Rho     = acov2corcpp(Sigma)
  Rho.est = apply(Rho,c(1,2),function(x) unlist(bm(x)))[1,,]
  Rho.se  = apply(Rho,c(1,2),function(x) unlist(bm(x)))[2,,]
  diag(Rho.se) = NA

        
  ### Fitted values and other things
  if(verbose){
    cat("=> Calculating linear predictors and fitted values. \n")
  }
  start <- Sys.time()
  
  linear.predictors = matrix(0,n,J)
  fitted.values     = matrix(0,n,J)
    
  for(j in 1:J){
    eta                   = X[[j]]%*%t(beta.chain[[j]]) + M[[j]]%*%Psi[,j,]
    lambda                = eta
    linear.predictors[,j] = rowMeans(eta)  
    fitted.values[,j]     = rowMeans(lambda)
  }
  fftime <- Sys.time()-start  
  if(verbose){
    fftime <- round(fftime,3)
    cat(paste0("=> Calculation of linear predictors and fitted values took ",fftime," seconds. \n"))
  } 
    
  ### Calculate residuals
  residuals = Y - fitted.values
    
  ### Calculate DIC
  D.bar = mean(v)
  pD    = D.bar + 2 * dGaussMulti(Y, X, M, coefficients, c(Psi.est), tau2e.est)
  dic   = D.bar + pD
    
  object = list(coefficients      = coefficients, 
                fitted.values     = fitted.values,
                linear.predictors = linear.predictors, 
                residuals         = residuals,
                beta.sample       = beta.chain, 
                Psi.sample        = Psi, 
                tau2e.sample      = tau2e, 
                Sigma.sample      = Sigma,
                beta.mcse         = beta.mcse, 
                Psi.mcse          = Psi.se, 
                tau2e.mcse        = tau2e.se, 
                Sigma.mcse        = Sigma.se,
                Psi.est           = Psi.est, 
                tau2e.est         = tau2e.est,      
                Sigma.est         = Sigma.est,
                Rho.est           = Rho.est,
                Rho.se            = Rho.se,
                iter              = running.iter-1, 
                dic               = dic,
                D.bar             = D.bar, 
                pD                = pD)
  class(object) = c("sparse.msglmm")
  object
}


