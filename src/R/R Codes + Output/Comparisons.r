setwd("C:/Users/admin/Desktop/Horseshoe")

set.seed(123)
library(horseshoe)
library(monomvn)
source("Bayesian_sqrtLasso.R")
library(scalreg)
library(glmnet)
source("DL_sqrt.R")
start_values = 1

compare <-function(a,b){
n = a
p = b
prop = tail(head(seq(0 , 1 , by = .1) , -1) , -1)
matbeta = matrix(0 , nrow = p , ncol = 9)
X = matrix(rnorm(n*p , 0 , sqrt(2)) , ncol = p) 
for(i in 1:9)  matbeta[,i] = c(rep(5 , prop[i]*n) , rep(0 , ceiling(p-prop[i]*n)) ) 

L <- as.list(1:9)
build <- function(X , beta) {
return(list(beta0 = beta , y = (crossprod(t(X) , beta) + rnorm(n , 0 , sqrt(10)))))
						}
L <- lapply(L  , function(x) build( X , matbeta[,x] ))



estim <- function(y){
out <- list()
#####
##### Horseshoe
#####
out$HS <- horseshoe(y , X , method.tau = "halfCauchy" , method.sigma = 'Jeffreys' , 
					burn = 1000 , nmc = 5000)
#Beta_Median <- apply(out$HS$BetaSamples , 1 , median)
#out$HS[['Beta_Median' 

#out$HSf <- horseshoe(y , X , method.tau = "fixed" , method.sigma = 'Jeffreys' , 
#					burn = 1000 , nmc = 5000)



####
#### Bayesian Lasso
####
#out$blasso <- blasso(X , y , T = 5000 , icept = FALSE )
#betahat<- apply(out.blasso$beta , 2 , mean)
#beta_median <- apply(out.blasso , 2 , median)
#out$blasso <- list(BetaSamples = out.blasso$beta , BetaHat = betahat , Sigma2 = out.blasso$s2 , 
#					Lambda2 = out.blasso$lambda2 , Tau2inverse = out.blasso$tau2i)

####
#### Lasso
####
lasso <- cv.glmnet(X , y , family = 'gaussian' , alpha = 1 , intercept = FALSE)
out$lasso <- coef.glmnet(lasso , s = lasso$lambda.1se)
out$lasso<- tail(as.vector(coef.glmnet(lasso , s = lasso$lambda.1se) ) , -1)
if(length(out$lasso) != p){ browser() }

####
#### Square_root Lasso
####
#slasso <- scalreg(X , y , lam0 = "univ")
#out$sqrtlasso <- slasso$coefficients

####
#### Bayesian square root lasso
####
#setwd("C:/Users/Desktop/Horseshoe/")

out$bsqlasso <- BSqrt_lasso( y , X , 1 , 1 , start_values , nmc = 5000 , nburn =  1000 , thin = 1)


###
### Square root Dirichlet Laplace
###

out$DL_s <- DL_square( y =y , X =  X ,a =1/p , nmc =5000 , nburn = 1000 , thin = 1 )
return(out)
}

output <- lapply(L , function(x) estim(x[[2]]))

ResMean<- list()
ResMean$HS <- lapply(output , function(x) x$HS$BetaHat)
#res$HSf <- lapply(output , function(x) x$HSf$BetaHat)
ResMean$DL_s <- lapply(output , function(x) x$DL_s$BetaHat )
ResMean$bsqlasso <- lapply(output , function(x) x$bsqlasso$BetaHat) 
ResMean$Lasso <- lapply(output , function(x) x$lasso)
#res$SqrtLasso <- lapply(output , function(x) x$sqrtlasso)

ResMedian<- list()
ResMedian$HS <- lapply(output , function(x) x$HS$BetaMedian)
#res$HSf <- lapply(output , function(x) x$HSf$BetaMedian)
ResMedian$DL_s <- lapply(output , function(x) x$DL_s$BetaMedian )
ResMedian$bsqlasso <- lapply(output , function(x) x$bsqlasso$BetaMedian) 
ResMedian$Lasso <- lapply(output , function(x) x$lasso)
#res$SqrtLasso <- lapply(output , function(x) x$sqrtlasso)


j = 1
f <- function(x ,j){ plot(x , xlab = "Index" , ylab = '' 
 ,main = bquote(paste( hat(beta)," with sparsity " )== .(prop[j] )))
evalq(j <- j+1 , envir = parent.frame(3) ) }

k = 1
g <- function(x , j) { x11()
					  par(mfrow = c(3,3) , mar = c(4.1 , 2.5 , 3.5 ,2 ) , 
					 oma = c(0,0,2,0) , cex.main = 1 , las = 1)
					 lapply(x,function(x) f(x ,j)) ; evalq(j <-1 , envir = parent.frame(2) )  
					 mtext(names(ResMean)[k] , outer = TRUE , cex = 1) ; evalq(k <- k+ 1 , envir = parent.frame(3) )
}

lapply(ResMean , function(x) g(x , j))

(list(Betahat = ResMean , Median = ResMedian , Output = output , data = list(L ,X) ))


#########################
#########################
##Some unfinished plots##
#########################
#########################
# pdf("Histograms.pdf" , paper = 'a4' , height = 11 )
# #output = res$Output
# j = 1
# i = 1
# f <- function(x ,j){ hist(x ,"FD"  , freq = F , xlab = '' , ylab = ''
 # ,main = bquote(paste( beta," with sparsity " )== .(prop[j] )))
 # }


# ff <- function(x ,j) { apply(x , 1 ,function(x) f(x ,j)) ; evalq(j <- j+1 , envir = parent.frame(6) )}

 
# k = 1
# g <- function(x , j) { #x11()
					  # par(mfrow = c(5,2) , mar = c(4.1 , 2.5 , 3.5 ,2 ) , 
					 # oma = c(0,0,2,0) , cex.main = 1 , las = 1)
					 # lapply(x,function(x) ff(x,j) ) ; evalq(j <-1 , envir = parent.frame(2) )  
					 # mtext(names(res$Betahat)[k] , outer = TRUE , cex = 1) ; evalq(k <-k+1 , envir = parent.frame(3) )
# }
# samples <- list()


# samples$HS <- lapply(output , function(x) x$HS$BetaSamples)
# samples$HSf <- lapply(output , function(x) x$HSf$BetaSamples)
# samples$blasso <- lapply(output , function(x) t(x$blasso$betaSamples) )
# samples$bsqlasso <- lapply(output , function(x) t(x$bsqlasso$Beta_samples)) 

# lapply(samples , function(x) g(x , j))
# dev.off()
}

# For n < p 

Sim1 <- compare(100 , 200 )

str(Sim1 , 1)

j = 1
f <- function(x ,j){ plot(x , xlab = "Index" , ylab = '' , col = "blue"
 ,main = bquote(paste( hat(beta)," with sparsity " )== .(prop[j] )))
evalq(j <- j+1 , envir = parent.frame(3) ) }

k = 1
g <- function(x , j) { x11()
					  par(mfrow = c(3,3) , mar = c(4.1 , 2.5 , 3.5 ,2 ) , 
					 oma = c(0,0,2,0) , cex.main = 1 , las = 1)
					 lapply(x,function(x) f(x ,j)) ; evalq(j <-1 , envir = parent.frame(2) )  
					 mtext(names(Sim1$Median)[k] , outer = TRUE , cex = 1) ; evalq(k <- k+ 1 , envir = parent.frame(3) )
}

lapply(Sim1$Median , function(x) g(x, j ))

