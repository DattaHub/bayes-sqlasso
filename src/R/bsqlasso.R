library(statmod)
library(MASS)
library(mvtnorm)

bsqlasso <- function(Y, X, method.tau = "fixed", r = 1, delta = 0.5, nmc=100, burn=100, 
                             eff_zero=1e-5, verbose=TRUE) {
  n <- nrow(X)
  p <- ncol(X)
  t = 1
  Tau = 1
  Theta = rep(0, p)
  Lambda = rep(1, p)
  tau.fixed = max(c(sum(Y>=sqrt(2*log(p)))/p, 1/p))
  # tau.fixed = 0.1
  
  ThetaSave = matrix(0, nrow=nmc, ncol=p)
  LambdaSave = matrix(0, nrow=nmc, ncol=p)
  KappaSave = matrix(0, nrow=nmc, ncol=p)
  TauSave = rep(0, nmc)
  tsave = rep(0, nmc)
  
  for(iter in 1:(nmc+burn)) {
    if(isTRUE(verbose) && iter %% 200 == 0) 
      cat("Iteration ",iter, "\n")

    A = crossprod(X)*t + diag(1/Lambda^2)
    sigma_t = chol2inv(A)
    mu_t = sigma_t%*%crossprod(X,Y)*t
    Theta_new = t(chol(sigma_t))%*%rnorm(p)+mu_t
    ## Theta_new = rmvnorm(n = 1, mu_t, sigma_t)
    
    # kappa = t * lam_t
    kappa = sigma_t%*%crossprod(X,rep(1,n))*t
    
    stopifnot(!any(is.nan(Theta_new)))
    Theta = Theta_new
    
    ## Now update Lambda, comment out for global shrinkage only
    Lambda2_new = rinvgauss(p, mean = abs(Theta/Tau), shape = Theta^2)
    
    ## Following Bayesian Lasso 
    # Lambda2_new = 1/rinvgauss(n, mean = abs(Tau/Theta), shape = Tau^2)
    
    # XXX: A hack to avoid 0/0 in Lambda / Tau.
    Lambda2_new = ifelse(Lambda2_new < eff_zero, eff_zero, Lambda2_new)
    
    # DEBUG: Remove. 
    #cat(sprintf("t=%d, mean(Lambda2 > 0 (eff))=%g \n", t, 
    #            mean(Lambda2_new > eff_zero)))
    stopifnot(!any(is.nan(Lambda2_new)))
    Lambda2 = Lambda2_new
    Lambda = sqrt(Lambda2)
    
    #cat(sum(Lambda2),"\n")
    
    ## browser()
    # Now Update t
    # Y_pred = drop(X%*%Theta)
    # t = rinvgauss(1, mean = sqrt(1/crossprod((Y-Y_pred))), shape = 1)
    # cat(t,"\n")
    
    # Now Update tau
    
    # Tau = ifelse(method.tau =="fixed", tau.fixed,
    #              sqrt(rgamma(1,n + r, delta+sum(Lambda2)/2)))
    
    if (method.tau == "fixed"){
      Tau <- tau.fixed
    } else if (method.tau == "gamma") {
      Tau <- sqrt(rgamma(1,p + r, delta+sum(Lambda2)/2))
    } else if (method.tau == "Cauchy") {
      
      eta = Tau^2
      u = runif(1, 0, 1/(eta + 1))
      ub = ifelse((1 - u) / u < 1, (1-u)/u, 1)
      a = n - 1
      b = 0.5*crossprod(Lambda)
      ub2 = pgamma(ub, a, rate=b)
      u2 = runif(1, 0, ub2)
      eta = qgamma(u2, a, rate=b)
      #browser()
      Tau_new = sqrt(eta)
      Tau_new = ifelse(Tau_new < eff_zero, eff_zero, Tau_new)
      
      # DEBUG: Remove. 
      cat(Tau_new,"\n")
      stopifnot(!any(is.nan(Tau_new)))
      Tau = Tau_new
    } else { Tau <- 1}
    
    # cat(delta+crossprod(Lambda)/2,tau_2)
    # cat(Tau,"\n")
    #if(anyNA(Nu)){ browser() }
    
    if(iter > burn)
    {
      ThetaSave[iter-burn,] = Theta
      LambdaSave[iter-burn,] = Lambda
      KappaSave[iter-burn,] = kappa
      TauSave[iter-burn] = Tau
      tsave[iter-burn] = t
    }
  }
  ThetaHat = apply(ThetaSave, 2, mean)
  LambdaHat = apply(abs(LambdaSave), 2, mean)
  TauHat = mean(TauSave)
  return(list(ThetaSave=ThetaSave, 
              LambdaSave=LambdaSave,
              KappaSave = KappaSave,
              TauSave=TauSave, 
              ThetaHat=ThetaHat,
              LambdaHat=LambdaHat,
              TauHat=TauHat
  ))
}
