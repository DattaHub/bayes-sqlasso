## MWE

set.seed(123)
beta = 20; sigmasq = 1
X = matrix(c(1,1),byrow=T)
y = X%*%beta + rnorm(2)
#library(horseshoe)

res <- horseshoe(y, X, method.tau = "truncatedCauchy",
                       method.sigma = "Jeffreys",
                 burn = 1000, nmc = 5000, alpha = 0.05)
#library(ggplot2)
samples.hs <- data.frame(values = t(unlist(res$BetaSamples)))
ggplot(data=samples.hs,aes(x=values))+
  geom_density(position="identity")

