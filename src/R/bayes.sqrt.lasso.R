library(statmod)

bayes.sqrt.lasso <- function(Y, method.tau = "fixed", r = 1, delta = 0.5, nmc=100, burn=100, 
                     eff_zero=1e-5, verbose=TRUE) {
  n <- length(Y)
  t = 1
  Tau = 1
  Theta = Y #rep(1, n)
  Lambda = Y #rep(1, n)
  tau.fixed = max(c(sum(Y>=sqrt(2*log(n)))/n, 1/n))
  # tau.fixed = 0.1
  
  ThetaSave = matrix(0, nrow=nmc, ncol=n)
  LambdaSave = matrix(0, nrow=nmc, ncol=n)
  KappaSave = matrix(0, nrow=nmc, ncol=n)
  TauSave = rep(0, nmc)
  tsave = rep(0, nmc)
  
  for(iter in 1:(nmc+burn)) {
    if(isTRUE(verbose) && iter %% 200 == 0) 
      cat("Iteration ",iter, "\n")
    
    
    lam_t = Lambda^2. / (1 + t*Lambda^2)
    Theta_new = rnorm(n, t * Y*lam_t , sqrt(lam_t)) 
    
    kappa = t * lam_t
    
    stopifnot(!any(is.nan(Theta_new)))
    Theta = Theta_new
    
    # Now update Lambda, comment out for global shrinkage only
    Lambda2_new = rinvgauss(n, mean = abs(Theta/Tau), shape = Theta^2)
    
    ## Following Bayesian Lasso 
    ## Lambda2_new = 1/rinvgauss(n, mean = abs(Tau/Theta), shape = Tau^2)
    
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
    ## t = rinvgauss(1, mean = sqrt(1/sum((Y-Theta)^2)), shape = 1)
    #cat(t,"\n")
    
    # Now Update tau
    
    # Tau = ifelse(method.tau =="fixed", tau.fixed,
    #              sqrt(rgamma(1,n + r, delta+sum(Lambda2)/2)))
    
    if (method.tau == "fixed"){
      Tau <- tau.fixed
    } else if (method.tau == "gamma") {
      Tau <- sqrt(rgamma(1,n + r, delta+sum(Lambda2)/2))
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