library(statmod)
library(MASS)
library(GIGrvg)
#####
#####Build data then apply following code 
####




DL_square <- function( y , X , a  , nmc , nburn , thin){

n = length(y)
p = ncol(X)
stopifnot(nrow(X) == n)
y <- matrix(y , ncol = 1)
Phi <- rep(1 , p)
Psi = rep(1,p)
Tau = 1
v_2 <- 1
eff_zero=1e-5

beta_save <- matrix(0 , nmc%/%thin , p)
Phi_save <- matrix(0 ,  nmc%/%thin  , p)
Psi_save <- matrix(0 ,  nmc%/%thin  , p)
v_save <- rep(0 ,  nmc%/%thin )
Tau_save <- rep(0 ,  nmc%/%thin )

for(i in 1:(nmc+nburn)) {
	if(i %% 200 == 0) 
	 cat("Iteration " , i , "\n")
	A <- chol2inv( chol((1/v_2)*crossprod(X , X) + diag( 1/(Tau *Phi*Psi) )))
	mu_beta = crossprod(t(A) , crossprod(X,y)) / v_2
	Sigma_beta <- A

	beta_new <- t(chol(A)) %*%matrix(rnorm(p) , ncol=1) + mu_beta
	#if(i > nburn)  browser()

	v_2_new <- 1 / rinvgauss( 1 , 1 / (1/Tau + base:::norm(y-X%*%beta_new , type = "2")) , shape = 1)

	#if( any( beta_new==0 ) ) browser()
	Psi_new <- 1/ rinvgauss(p , mean = Phi*Tau / abs(beta_new) , shape =1)

	tau_new <- 1/rgig(1, chi = 2*sum(abs(beta_new)/Phi) + 1/v_2_new, psi = 1, lambda = 1-p)

	 T <- matrix(beta_new , nrow =1)
	 T <- as.vector(apply(T , 2 ,function(x) rgig(1, chi = 2*abs(x), psi = 1, lambda = 1/p - 1)))
	 T <- T/sum(T)
	 Phi_new <- 1/T


	#if( any (tau_new==0)||any( is.finite(tau_new) ==F)) browser()
	#tau_new = ifelse(tau_new < eff_zero , eff_zero , tau_new)

#	lambda_new <- sqrt(rgamma(1 , shape =  p+r , rate = delta + sum(tau_new) /2 ))
if(i>nburn ){ j = i-nburn

if (j%%thin == 0) { k = j%/%thin 
	beta_save[k,] <- beta_new
	Tau_save[k] <- tau_new
	v_save[k] <- v_2_new
	Psi_save[k,] <- Psi_new
	Phi_save[k,] <- Phi_new
}

}
	v_2 <- v_2_new
	Tau <- tau_new
	Psi <- Psi_new
	Phi <- Phi_new





}


beta_hat <- apply(beta_save , 2 , mean)
beta_median <- apply(beta_save , 2 , median)
leftCI <- apply(beta_save , 2 ,function(x) quantile(x , .025 ))
rightCI <- apply(beta_save , 2 ,function(x) quantile(x , .975 ))

return(output <- list(BetaHat = beta_hat, BetaMedian = beta_median ,BetaSamples = beta_save , Tau2Samples = Tau_save , V2Samples = 
				v_save , PsiSamples = Psi_save , PhiSamples = Phi_save , LeftCI = leftCI , RightCI = rightCI))





}
