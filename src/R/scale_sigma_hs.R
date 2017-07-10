## MWE

set.seed(123)
beta = 20; sigmasq = 1
X = matrix(c(1,1),byrow=T)
y = X%*%beta + rnorm(2)
library(horseshoe)

res <- horseshoe(y, X, method.tau = "truncatedCauchy",
                       method.sigma = "Jeffreys",
                 burn = 1000, nmc = 5000, alpha = 0.05)
library(ggplot2)
samples.hs <- data.frame(values = t(unlist(res$BetaSamples)))
ggplot(data=samples.hs,aes(x=values))+
  geom_density(position="identity")

setwd("C:/Users/Jyotishka/OneDrive/Documents/R/sqlasso")
source("eval_sql.R")

Y = y
res2 <- eval_sql(Y,r = 1, delta = 2, burn = 1000, nmc = 5000)
samples.sql <- data.frame(values = as.vector(res2$ThetaSave))
ggplot(data=samples.sql,aes(x=values))+
  geom_density(position="identity")

plot(res2$ThetaSave[,1])

setwd("C:/Users/Jyotishka/OneDrive/Documents/R/sqlasso")
dev.copy2pdf(file="theta_bimodal.pdf")