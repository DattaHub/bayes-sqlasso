#############
############# Code for Bayesian Square root Lasso 
#############			Gibbs sampler			

library(statmod)
library(MASS)

#####
#####Build data then apply following code 
####




BSqrt_lasso <- function( y , X , r , delta , start_values , nmc , nburn, thin){

n = length(y)
p = ncol(X)
stopifnot(nrow(X) == n)
y <- matrix(y , ncol = 1)
Tau <- rep(1 , p)
lambda = 1
v_2 <- 1
eff_zero=1e-5

beta_save <- matrix(0 , nmc , p)
tau_save <- matrix(0 , nmc  , p)
v_save <- rep(0 , nmc )
lambda_save <- rep(0 , nmc )

for(i in 1:(nmc+nburn)) {
	if(i %% 2000 == 0) 
	 cat("Iteration " , i , "\n")
	A <- chol2inv( chol((1/v_2)*crossprod(X , X) + diag( 1/Tau ) ))
	mu_beta = crossprod(t(A) , crossprod(X,y)) / v_2
	Sigma_beta <- A

	beta_new <- t(chol(A)) %*%matrix(rnorm(p) , ncol=1) + mu_beta

	v_2_new <- 1 / rinvgauss( 1 , 1 / base:::norm(y-X%*%beta_new , type = "2") , shape = 1)

	#if( any( beta_new==0 ) ) browser()

	tau_new <- 1/ rinvgauss(p , mean = lambda / abs(beta_new) , shape = lambda^2)

	if( any (tau_new==0)||any( is.finite(tau_new) ==F)) browser()
	#tau_new = ifelse(tau_new < eff_zero , eff_zero , tau_new)

	lambda_new <- sqrt(rgamma(1 , shape =  p+r , rate = delta + sum(tau_new) /2 ))
if(i>nburn ){ j = i-nburn

if (xor(thin==1,j%%thin == 1)) { k = j%/%thin
	beta_save[k,] <- beta_new
	tau_save[k,] <- tau_new
	v_save[k] <- v_2_new
	lambda_save[k] <- lambda_new
}
}
	v_2 <- v_2_new
	Tau <- tau_new
	lambda <- lambda_new





}


beta_hat <- apply(beta_save , 2 , mean)
beta_median <- apply(beta_save , 2 , median)
leftCI <- apply(beta_save , 2 ,function(x) quantile(x , .025 ))
rightCI <- apply(beta_save , 2 ,function(x) quantile(x , .975 ))


return(output <- list(BetaHat = beta_hat, BetaMedian = beta_median ,BetaSamples = beta_save , Tau2Samples = tau_save , V2Samples = 
				v_save , Lambda2Samples = lambda_save , LeftCI = leftCI , RightCI = rightCI))





}


class_vector <- function(x){ t=x
					  km = head(kmeans(abs(x) ,2) , 2)
					  x = km$clust
					  #center = km$centers
					  #clust  = km$clust
					  signal <- which.max(km$centers)
		    			  noise <- which.min(km$centers)
					  t[which(x==noise)] <- 0
		                  t[which(x==signal)] <- 1
					  return(t)
}

