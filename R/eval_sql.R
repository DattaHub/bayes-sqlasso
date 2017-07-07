library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(statmod)

eval_sql <- function(Y, r = 1, delta = 2, nmc=100, burn=100, 
                        eff_zero=1e-5, verbose=TRUE) {
  n <- length(Y)
  t = 1
  Tau = 1
  Theta = rep(1, n)
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
    
    
    lam_t = 2*Lambda^2. / (2 + t*Lambda^2)
    Theta_new = rnorm(n, t * Y*lam_t / 2, 
                      sqrt(lam_t)) 
    
    stopifnot(!any(is.nan(Theta_new)))
    Theta = Theta_new
    
    # Now update Lambda, comment out for global shrinkage only
    Lambda2_new = rinvgauss(n, mean = abs(Theta/Tau), shape = Theta^2)
    
    # XXX: A hack to avoid 0/0 in Lambda / Tau.
    Lambda2_new = ifelse(Lambda2_new < eff_zero, eff_zero, Lambda2_new)
    
    # DEBUG: Remove. 
    #cat(sprintf("t=%d, mean(Lambda2 > 0 (eff))=%g \n", t, 
    #            mean(Lambda2_new > eff_zero)))
    stopifnot(!any(is.nan(Lambda2_new)))
    Lambda2 = Lambda2_new
    
    Lambda = sqrt(Lambda2)
    ## browser()
    # Now Update t
    t = rinvgauss(1, mean = 2/sum((Y-Theta)^2), shape = 1)
    
    # Now Update tau 
    tau_2 = rgamma(1,n + r, delta+sum(Lambda^2)/2)
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

# Ex 1
theta = c(rep(3,10),rep(0,90))
Y = rnorm(100,theta,0.1)
ans1 = eval_sql(Y)

plot(Y)
lines(ans1$ThetaHat)

# Ex 2
theta = c(rep(7,10),rep(3,10),rep(0,80))
Y = rnorm(100,theta,1)
ans2 = eval_sql(Y)

plot(Y)
lines(ans2$ThetaHat)


##############################
########## Fancy Plots #######
##############################


set.seed(198)
theta_1 = c(rep(7, 10), rep(0, 90))
Y_1 = rnorm(100, theta_1, 1)
Y_2 = rnorm(100, theta_1, 0.5)

theta_2 = c(rep(7, 10), rep(3,10), rep(0, 80))
Y_3 = rnorm(100, theta_2, 1)
Y_4 = rnorm(100, theta_2, 0.5)

prob_set = list(list(name='theta_1, sigma=1', obs=Y_1),
                list(name='theta_1, sigma=2', obs=Y_2),
                list(name='theta_2, sigma=1', obs=Y_3),
                list(name='theta_3, sigma=2', obs=Y_4))  

theta.data = lapply(prob_set,
                    function(set) 
                      list('settings'=set$name,
                           'samples'=eval_sql(set$obs, 
                                                 nmc=5000)))

theta.data = lapply(theta.data,
                    function(sample.data)
                      cbind(melt(sample.data$samples$ThetaSave), 
                            'settings'=sample.data$settings, 
                            'var'='theta'))

theta.data = Reduce(rbind, theta.data)

colnames(theta.data)[1:2] = c('sample', 'component')

obs.data = lapply(prob_set, 
                  function(set) 
                    data.frame(component=seq_along(set$obs), 
                               value=set$obs,
                               sample=NA, 
                               settings=set$name, 
                               var="obs")) 
obs.data = Reduce(rbind, obs.data)

all.data = rbind(obs.data, theta.data)

plot.data = all.data %>% 
  group_by(component, settings, var) %>%
  summarise(upper = quantile(value, prob=0.95), 
            lower = quantile(value, prob=0.05), 
            middle = median(value))

theta.plot = ggplot(plot.data,
                    aes(x=component, y=middle, group=component,
                        colour=var)) + theme_bw() +   
  geom_pointrange(aes(ymin=lower, ymax=upper), size=0.2, alpha=0.5) + 
  facet_grid(settings ~ ., scales="free_y") + #, labeller = label_both) +
  scale_colour_manual(values=c("#D55E00","#0072B2")) +
  xlab("") + ylab("")

print(theta.plot)
#setwd("C:/Users/Jyotishka/OneDrive/Documents/R/sqlasso")
ggsave("sparse-means-sql-1.pdf", plot=theta.plot)

