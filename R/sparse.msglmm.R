
sparse.msglmm.fit.binomial = function(Y, X, A, M, beta.start, V, offset, tol, minit, maxit, sigma.s, sigma.b, verbose)
{
    iterations = minit
    n = length(Y)
    Q = t(M) %*% (diag(rowSums(A), n) - A) %*% M
    p = ncol(X)
    beta = matrix(0, iterations + 1, p)
    beta[1, ] = beta.start
    tau.s = numeric(iterations + 1)
    tau.s[1] = 0.1
    q = ncol(M)
    gamma = matrix(0, iterations + 1, q)
    gamma[1, ] = rnorm(q, 0, 1)
    a.s = 0.5
    b.s = 2000

    
    if (is.null(offset))
        offset = rep(0, n)
    start = 1
    k = 1
    five.pct = round(0.05 * maxit, 0)
    repeat
    {
        for (j in (start + 1):(start + iterations))
        {
            k = k + 1
            if (verbose && k %% five.pct == 0)
                cat("Progress => ", round(k / maxit * 100, 0), "%\n", sep = "")


            #### MCMC sampling
            
            
            if (j == maxit)
                break
        }
        if (j == maxit)
        {
            beta = as.matrix(beta[1:maxit, ])
            gamma = as.matrix(gamma[1:maxit, ])
            tau.s = tau.s[1:maxit]
            break
        }
        done = TRUE
        for (j in 1:p)
        {
            temp = bm(beta[, j])
            if (temp$se > tol)
            {
                done = FALSE
                break
            }
        }
        if (done)
        {
            for (j in 1:q)
            {
                temp = bm(gamma[, j])
                if (temp$se > tol)
                {
                    done = FALSE
                    break
                }
            }
        }
        if (done)
        {
            temp = bm(tau.s)
            if (temp$se > tol)
                done = FALSE
        }
        if (done)
            break
        else
        {
            start = start + iterations
            temp = matrix(0, iterations, p)
            beta = rbind(beta, temp)
            temp = matrix(0, iterations, q)
            gamma = rbind(gamma, temp)
            tau.s = c(tau.s, rep(0, iterations))
        }
    }
    coefficients = numeric(p)
    beta.mcse = numeric(p)
    names(coefficients) = names(beta.mcse) = colnames(X)
    for (j in 1:p)
    {
        temp = bm(beta[, j])
        coefficients[j] = temp$est
        beta.mcse[j] = temp$se
    }
    gamma.est = numeric(q)
    gamma.mcse = numeric(q)
    for (j in 1:q)
    {
        temp = bm(gamma[, j])
        gamma.est[j] = temp$est
        gamma.mcse[j] = temp$se
    }
    temp = bm(tau.s)
    tau.s.est = temp$est
    tau.s.mcse = temp$se
    linear.predictors = numeric(n)
    fitted.values = numeric(n)
    iter = length(tau.s)
    v = numeric(iter)
    for (j in 1:iter)
    {
        eta = offset + X %*% beta[j, ] + M %*% gamma[j, ]
        linear.predictors = linear.predictors + eta / iter
        p = linkinv(eta)
        fitted.values = fitted.values + p / iter
        v[j] = -2 * sum(dbinom(Y, 1, p, log = TRUE))
    }
    D.bar = mean(v)
    pD = D.bar + 2 * sum(dbinom(Y, 1, fitted.values, log = TRUE))
    dic = D.bar + pD
    residuals = Y - fitted.values
    beta.accept = sum(diff(beta[, 1]) != 0) / iter
    gamma.accept = sum(diff(gamma[, 1]) != 0) / iter
    object = list(coefficients = coefficients, fitted.values = fitted.values,
                  linear.predictors = linear.predictors, residuals = residuals,
                  beta.sample = beta, gamma.sample = gamma, tau.s.sample = tau.s,
                  beta.mcse = beta.mcse, gamma.mcse = gamma.mcse, tau.s.mcse = tau.s.mcse,
                  gamma.est = gamma.est, tau.s.est = tau.s.est, iter = iter, dic = dic,
                  D.bar = D.bar, pD = pD, beta.accept = beta.accept, gamma.accept = gamma.accept)
    class(object) = c("sparse.sglmm")
    object
}




#' Fit a sparse multivariate SGLMM.
#'
#' @details This function does things
#' @param X the design matrix.
#' @param A the adjacency matrix for the underlying graph. The matrix need not be binary, but it must be numeric and symmetric.
#' @param theta the vector of parameter values: \eqn{\theta = (\beta^\prime, \eta)^\prime}{\theta = (\beta', \eta)'}.
#' @return A vector that is distributed exactly according to the centered autologistic model with the given design matrix and parameter values.
#' @export
sparse.msglmm <- function(family, iterations, Y, X, M, offset, sigmab, sigmas, taub, 
                          beta0, Psi0, Sigma0, verbose=1){

  
  if(family=="poisson"){
    msglmmPoisson$msglmmPoissonCPP(iterations, Y, X, M, offset, 
                                   sigmab, sigmas, taub, 
                                   beta0, Psi0, Sigma0, verbose)
  } else if(family=="binomial"){
    n       <- nrow(Y); J <- ncol(Y)
    Ntrials <- matrix(1,n,J)
    msglmmBinomial$msglmmBinomialCPP(iterations, Y, Ntrials, X, M, offset, 
                                     sigmab, sigmas, taub, 
                                     beta0, Psi0, Sigma0, verbose)
  } 
}


