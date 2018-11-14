## Linear Regression
setwd("C:/Users/Jyotishka/OneDrive/Documents/R/sqlasso")

library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

## Ex 1

# (note that bayes.sqrt.lasso is much faster in this case)
set.seed(244)
X <- diag(100)
beta <- c(rep(0, 80), rep(8, 20))
Y <- beta + rnorm(100)

# 
source("bsqlasso.R")
ans1 = bsqlasso(Y, X, method.tau ="fixed", r = 1, delta = 1/2, burn = 1000, nmc = 5000)

# 
par(mfrow=c(1,1))
plot(beta)
lines(ans1$ThetaHat)
points(apply(ans1$ThetaSave,2,median),col="red",pch=15)
# 
