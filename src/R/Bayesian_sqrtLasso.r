#############
############# Code for Bayesian Square root Lasso 
#############			Gibbs sampler			

library(statmod)
library(MASS)

#####
#####Build data then apply following code 
####




BSqrt_lasso <- function( y , X , r , delta , start_values , nmc , nburn){

n = length(y)
p = ncol(X)
stopifnot(nrow(X) == n)
y <- matrix(y , ncol = 1)
Tau <- rep(1 , p)
lambda = 1
v_2 <- 1
eff_zero=1e-5

beta_save <- matrix(0 , nmc + nburn, p)
tau_save <- matrix(0 , nmc  + nburn, p)
v_save <- rep(0 , nmc + nburn)
lambda_save <- rep(0 , nmc + nburn)

for(i in 1:(nmc+nburn)) {
	if(i %% 200 == 0) 
	 cat("Iteration " , i , "\n")
	A <- chol2inv( chol((1/v_2)*crossprod(X , X) + diag( 1/Tau ) ))
	mu_beta = crossprod(t(A) , crossprod(X,y)) / v_2
	Sigma_beta <- A

	beta_new <- t(chol(A)) %*%matrix(rnorm(p) , ncol=1) + mu_beta

	v_2_new <- 1 / rinvgauss( 1 , 1 / norm(y-X%*%beta_new , type = "2") , shape = 1)

	#if( any( beta_new==0 ) ) browser()

	tau_new <- 1/ rinvgauss(p , mean = lambda / abs(beta_new) , shape = lambda^2)

	if( any (tau_new==0)||any( is.finite(tau_new) ==F)) browser()
	#tau_new = ifelse(tau_new < eff_zero , eff_zero , tau_new)

	lambda_new <- sqrt(rgamma(1 , shape =  p+r , rate = delta + sum(tau_new) /2 ))

	beta_save[i,] <- beta_new
	tau_save[i,] <- tau_new
	v_save[i] <- v_2_new
	lambda_save[i] <- lambda_new

	v_2 <- v_2_new
	Tau <- tau_new
	lambda <- lambda_new





}


beta_hat <- apply(beta_save[-(1:nburn),] , 2 , mean)


return(output <- list(Beta_hat = beta_hat ,Beta_samples = beta_save , Tau2_samples = tau_save , v2_samples = 
				v_save , Lambda2_samples = lambda_save))





}

