setwd("C:/Users/admin/Desktop/Horseshoe")

set.seed(123)
remove(list = ls())
X = matrix(rnorm(40*20 , 2 , 1) ,nrow = 40 , ncol =  20)
beta0 <- c(rep(5,6) , rep(0 , 14))
y = crossprod(t(X) , beta0) + rnorm(40 , 0 ,sqrt(5))
source("Bayesian_sqrtLasso.R")

res <- BSqrt_lasso(y , X , 0 , 0 , 0 , 2000 , 400)

beta_hat <- apply(res$Beta_samples[-(1:400),] , 2 , mean)
beta_hat
plot(res$Beta_samples[,1])
hist(Beta_samples[-(1:400),1] , "FD" , freq = F)


#############################


#Bayesian lasso codes
set.seed(123)
#remove(list = ls())
source("Bayes_lasso.R")

res <- Bayes_lasso(y , X , lambda.meth = "gamma" , 10 , 1.78 , 600 , 3000 , 0)

attach(res)
plot(res$Beta_hat)
plot(res$Beta_samples[,3])
plot(res$Lambda2_samples)
plot(res$Tau2_samples[,1])
