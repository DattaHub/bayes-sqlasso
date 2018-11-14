library(GeneralizedHyperbolic)
set.seed(123)
n = 100
p = 100
prop = tail(head(seq(0 , 1 , by = .1) , -1),-1)
matbeta = matrix(0 , nrow = p , ncol =9 )
X = matrix(rnorm(n*p , 0 , sqrt(2)) , ncol = p) 
for(i in 1:9)  matbeta[,i] = c(rep(5 , prop[i]*p) , rep(0 , ceiling(p-prop[i]*p)) ) 

L <- as.list(1:9)
build <- function(X , beta) {
return(list(beta0 = beta , y = (crossprod(t(X) , beta) + rnorm(n , 0 , sqrt(5)))))
						}
L <- lapply(L  , function(x) build( X , matbeta[,x] ))



setwd("C:/Users/admin/Desktop/Horseshoe")

source("DL_sqrt.R")
#x <- L[[1]]$y 
res <- lapply(L , function(x) DL_square( y = x[[2]] , X =  X ,a =1/p , nmc =8000 , nburn = 1000 , thin = 2 ))
res2<- lapply(L , function(x) horseshoe(y = x[[2]], X, method.tau = "truncatedCauchy",
                       method.sigma = "Jeffreys", burn = 1000, nmc = 8000 , thin = 2))
mean_beta <- lapply(res , function(x) x$BetaHat )
median_beta <- lapply(res , function(x) apply(x$BetaSamples , 2 , median))
tau_dl <- lapply(res, function(x) 1/sqrt(x$Tau2Samples))
tau_hs <- lapply(res2, function(x) x$TauSamples)
#tau_dl <- t(1/simplify2array(tau_dl))
#tau_hs <- t(simplify2array(tau_hs))
taudl <- as.data.frame(tau_dl) ; taudl[,"Method"] <- rep("Sqrt-DL" , 4000) ; names(taudl)[1:9] <- (1:9)/10
tauhs <- as.data.frame(tau_hs) ; tauhs[,"Method"] <- rep("Horseshoe" , 4000) ; names(tauhs)[1:9] <- (1:9)/10
dd<- rbind(taudl , tauhs)
#library(reshape2)
dd<- melt(dd , id.vars = "Method") ; names(dd)[c(2,3)] <- c("Sparsity" , "Tau") 
library(ggplot2)
x11()
tau_box <- ggplot(dd[(dd[,"Method"]=="Sqrt-DL"),] , aes(by = Method , x=Sparsity , y = Tau )) + geom_boxplot() + facet_grid(Method~.) + theme_bw()
tau_box
str(median_beta , 1)

out <- list(Means = mean_beta , Medians = median_beta )

j = 1
f <- function(x ,j){ plot(x , xlab = "Index" , ylab = '' 
 ,main = bquote(paste( hat(beta)," with sparsity " )== .(prop[j] )))
evalq(j <- j+1 , envir = parent.frame(3) ) }

k = 1
g <- function(x , j) { x11()
					  par(mfrow = c(3,3) , mar = c(4.1 , 2.5 , 3.5 ,2 ) , 
					 oma = c(0,0,2,0) , cex.main = 1 , las = 1)
					 lapply(x,function(x) f(x ,j)) ; evalq(j <-1 , envir = parent.frame(2) )  
					 mtext(names(out)[k] , outer = TRUE , cex = 1) ; evalq(k <- k+ 1 , envir = parent.frame(3) )
}

lapply(out , function(x) g(x , j))


leftCI <- lapply(res , function(x) apply(x$BetaSamples , 2 ,function(z) quantile(z , 0.025)))
rightCI <- lapply(res , function(x) apply(x$BetaSamples , 2 ,function(z) quantile(z , 0.975)))

j = 1
f1 <- function(x ,j){ plot(x , xlab = "Index" , ylab = '' , ylim = c(-14 ,14)
 ,main = bquote(paste( 'CI for ', hat(beta)," with sparsity " )== .(prop[j] )))
points(rightCI[[j]] , col = "blue")
evalq(j <- j+1 , envir = parent.frame(3) ) }

k = 1
g1 <- function(x , j) { x11()
					  par(mfrow = c(2,2) , mar = c(4.1 , 2.5 , 3.5 ,2 ) , 
					 oma = c(0,0,2,0) , cex.main = 1 , las = 1)
					 lapply(x,function(x) f(x ,j)) ; evalq(j <-1 , envir = parent.frame(2) )  
					 mtext(names(out)[k] , outer = TRUE , cex = 1) ; evalq(k <- k+ 1 , envir = parent.frame(3) )
}
par(mfrow = c(2,2))
lapply(leftCI , function(x) f1(x , j))



###
### Some diagnostics
###

Tau2 <- lapply(res , function(x) x$Tau2Samples)

Tau2 <- as.data.frame(Tau2)

boxplot(Tau2)