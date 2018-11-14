#######
####### Bayesian Lasso Gibbs sampler
#######

library(statmod)
library(MASS)

Bayes_lasso <- function(y , X , lambda.meth = "gamma"  , r , delta ,
nburn , nmcmc , niterEM) {
	 
	 n = length(y)
	 p = ncol(X)
	 Tau = rep(1,p)
	 if(n>p ){
	    beta_ls <- crossprod(t( chol2inv(chol(crossprod(X,X))) ) , crossprod(X , y))
	    sigma2 = sum( crossprod( ( y - crossprod(t(X) , beta_ls) ) ) )/(n-p)
	     lambda = p*sqrt(sigma2)/sum(abs(beta_ls))
		 } else{
	    sigma2 = 1
	     lambda = 1
		 }
	
	beta_save <- matrix(0 , nmcmc + nburn, p)
	tau_save <- matrix(0 , nmcmc  + nburn, p)
	sigma2_save <- rep(0 , nmcmc + nburn)
	lambda_save <- rep(0 , nmcmc + nburn)
	
	 y = y - mean(y)*rep(1,n)
	 for(i in 1:(nburn + nmcmc)){
		 if(i %% 200 ==0)
		  cat("Iteration" , i , "\n")
		 A = chol2inv( chol(crossprod(X,X) + diag(1/Tau) ))
		 Sigma = sigma2 * A
		 mu = crossprod( t(A) , crossprod(X,y) )
		 
		 beta_new = matrix(sqrt(sigma2)*crossprod(chol(A) , rnorm(p) ) + mu , ncol=1)
		 
		 T = rinvgauss( p, mean =  sqrt((lambda^2)* sigma2)/abs( beta_new )  , shape = lambda^2)
		 #if(length(T)!=p)) browser()
		 Tau_new = 1/T
		 Tau_new = ifelse(Tau_new > 1e2 , 1e2 , Tau_new)
		 if(length(Tau_new)!=p) browser()
		 a = y - crossprod(t(X) , beta_new)
		 a2 = t(a)%*%a
		 b = t(beta_new)%*% crossprod(diag(T), beta_new ) 
		 sigma2_new = 1/rgamma(1 , shape = (n+p-1)/2 , scale = .5*(a2 + b))
		 
		 
		beta_save[i,] <- beta_new
		tau_save[i,] <- Tau_new
		#lambda_save[i] <- lambda_new
		sigma2_save [i] <- sigma2_new
		
		Tau <- Tau_new
		#
		 sigma2 = sigma2_new
		 #lamda_new = 1
#		if (lambda.meth =="gamma"){
			
			lambda_new = sqrt(rgamma(1 , shape = p+r , rate = delta + .5*sum(Tau)))
#			}
		 lambda_save[i] <- lambda_new
		 lambda <- lambda_new
		 
		 
		 }
	 beta_hat <- apply(beta_save[-(1:nburn),] , 2 , mean)
	lambda_hat <- tail(lambda_save , 1)
return(output <- list(Beta_hat = beta_hat ,Beta_samples = beta_save , Tau2_samples = tau_save
				, Lambda2_samples = lambda_save , lambda_hat = lambda_hat , Sigma2 = sigma2_save))
	 }
	 
	 
	 		 # if(lambda.meth == "Ebayes"){
						 # gibbs_lambda <- function(sigma2 , Tau , lambda , y , X ,niter){
			 # for(i in 1:niter){
			# #if(i %% 200 ==0)
			  # # # cat("Iteration" , i , "\n")
			 # A = chol2inv( chol(crossprod(X,X) + 1/Tau ) )
			 # Sigma = sigma2 * A
			 # mu = crossprod( t(A) , crossprod(X,y) )
			 
			 # beta_new = sqrt(sigma2)*crossprod(chol(A) , rnorm(p) ) + mu
			 
			 # Tau_new = 1 / rinvgauss( p , sqrt(lambda^2 * sigma2)/( beta_new^2 )  , lamda^2)
			 # Tau_new = ifelse(is.finite(Tau_new)==F || Tau_new < 1e-5 , 1e-5 , Tau_new)
			 
			 # sigma2_new = 1/rgamma(1 , shape = (n+p-1)/2 , scale = .5*(crossprod( ( y - crossprod(t(X) 
			 # , beta_new) ) )  + t(beta_new)%*% crossprod(diag(1/Tau_new), beta_new ) ))
			
			# tau_save <- matrix(0 , niter, p)
			# tau_save[i,] = Tau_new}
			# return(sum(tau_hat = apply(tau_save , 2 , mean)))
			# }
			# lambda_new = sqrt( 2*P / (gibbs_lambda(sigma2 , Tau , lambda , y , X ,niterEM)) )
			
			# }