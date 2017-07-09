library(statmod)

sql.normal.means <- function(Y, r = 1, delta = 2, nmc=100, burn=100, 
                     eff_zero=1e-5, verbose=TRUE) {
  n <- length(Y)
  t = 1
  Tau = 1
  Theta = Y #rep(1, n)
  Lambda = rep(1, n)
  r = 1
  delta  = 2
  
  ThetaSave = matrix(0, nrow=nmc, ncol=n)
  LambdaSave = matrix(0, nrow=nmc, ncol=n)
  TauSave = rep(0, nmc)
  tsave = rep(0, nmc)
  
  for(iter in 1:(nmc+burn)) {
    if(isTRUE(verbose) && iter %% 200 == 0) 
      cat("Iteration ",iter, "\n")
    
    
    lam_t = Lambda^2. / (1 + t*Lambda^2)
    Theta_new = rnorm(n, t * Y*lam_t , sqrt(lam_t)) 
    
    stopifnot(!any(is.nan(Theta_new)))
    Theta = Theta_new
    
    # Now update Lambda, comment out for global shrinkage only
    Lambda2_new = rinvgauss(1, mean = abs(Theta/Tau), shape = Theta^2)
    # XXX: A hack to avoid 0/0 in Lambda / Tau.
    Lambda2_new = ifelse(Lambda2_new < eff_zero, eff_zero, Lambda2_new)
    # 
    # Lambda2_new = rep(0,n)
    # for(i in 1:n)
    # {
    #   Lambda2_new[i] = rinvgauss(1, mean = abs(Theta[i]/Tau), shape = Theta[i]^2)
    #   # XXX: A hack to avoid 0/0 in Lambda / Tau.
    #   Lambda2_new[i] = ifelse(Lambda2_new[i] < eff_zero, eff_zero, Lambda2_new[i]) 
    # }
    
    # DEBUG: Remove. 
    #cat(sprintf("t=%d, mean(Lambda2 > 0 (eff))=%g \n", t, 
    #            mean(Lambda2_new > eff_zero)))
    stopifnot(!any(is.nan(Lambda2_new)))
    Lambda2 = Lambda2_new
    Lambda = sqrt(Lambda2)
    
    #cat(sum(Lambda2),"\n")
    
    ## browser()
    # Now Update t
    t = rinvgauss(1, mean = sqrt(1/sum((Y-Theta)^2)), shape = 1)
    #cat(t,"\n")
    
    # Now Update tau 
    tau_2 = rgamma(1,n + r, delta+sum(Lambda2)/2)
    # cat(delta+crossprod(Lambda)/2,tau_2)
    Tau = sqrt(tau_2)
    #if(anyNA(Nu)){ browser() }
    
    if(iter > burn)
    {
      ThetaSave[iter-burn,] = Theta
      LambdaSave[iter-burn,] = Lambda
      TauSave[iter-burn] = Tau
      tsave[iter-burn] = t
    }
  }
  
  ThetaHat = apply(ThetaSave, 2, mean)
  LambdaHat = apply(abs(LambdaSave), 2, mean)
  TauHat = mean(TauSave)
  return(list(ThetaSave=ThetaSave, 
              LambdaSave=LambdaSave, 
              TauSave=TauSave, 
              ThetaHat=ThetaHat,
              LambdaHat=LambdaHat,
              TauHat=TauHat
  ))
}