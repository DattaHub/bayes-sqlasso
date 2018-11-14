library(glmnet)
library(monomvn)
library(bayesm)
library(plyr)
library(lars)
library(horseshoe)
library("parallel")
library("foreach")
library("doParallel")
library(tcltk)
library(statmod)
library(MASS)
library(GIGrvg)
library(scalreg)
source("DL_sqrt.R")
source("Bayesian_sqrtLasso.R")
setwd("C:/Users/admin/Desktop/Horseshoe")

no_cores <- detectCores()-1
RNGkind("L'Ecuyer-CMRG")
cl <- makeCluster(no_cores)
registerDoParallel(cl)
clusterSetRNGStream(cl, 13)
clusterEvalQ(cl, .libPaths("C:/Users/admin/Documents/R/win-library/Mylib"))
clusterEvalQ(cl , setwd("C:/Users/admin/Desktop/Horseshoe"))
clusterEvalQ(cl , source("DL_sqrt.R"))
clusterEvalQ(cl , source("Bayesian_sqrtLasso.R"))

niter = 100

# Uncomment this if you want a progress bar
# pb <- winProgressBar(title="Progress bar", label="0% done", min=0, max=100, initial=0)

n=100;p=60
p1 = 7
beta <- c(c(7,5,5,4,4,3,3),rep(0,p-p1))
index <- c(rep(1,p1),rep(0,p-p1))

res<- foreach(i = 1:niter, .packages = c("bayesm","horseshoe","glmnet" , "MASS" , "statmod",
							    "GIGrvg" ,"tcltk") , .verbose = TRUE
								)%dopar%{
      if(!exists("pb")) pb <- tkProgressBar("Parallel task", min=1, max=niter)
    setTkProgressBar(pb, i)
    sigma <- rwishart(p,diag(p))
    
    x=matrix(rnorm(n*p),n,p) # generate X matrices
    R <- chol(sigma$W)
    x <- x%*%R
    fx <- x%*% beta
    A <- t(x[,-(1:p1)])%*%x[,(1:p1)]%*%solve((t(x[,(1:p1)])%*%x[,(1:p1)]))%*%rep(1,p1)
    n_inf <- 1-norm(A,'i')
    reps = 50
    mp_hs = rep(0,reps)
    mp_lasso = rep(0, reps)
    mp_Bsql= rep(0, reps)
    mp_DL_sqrt = rep(0, reps)
    mse_hs = rep(0,reps)
    mse_lasso = rep(0, reps)
    mse_Bsql= rep(0, reps)
    mse_DL_sqrt = rep(0, reps)
    #bsq <- matrix(0 , reps , p)
    ptm <- c(proc.time()["elapsed"] , rep(0 ,reps))
    for (j in 1:reps)
    { 
      eps=sqrt(5)*rnorm(n)
      y=drop(fx+eps)
      data = data.frame(y,x)
      ## Using horseshoe
      res1 <- horseshoe(y, x, method.tau = "truncatedCauchy",
                       method.sigma = "Jeffreys", burn = 1000, nmc = 4000)
      hs_ind <- HS.var.select(res1, y, method = "intervals")
      mp_hs[j] <- sum((hs_ind!=index))
	  mse_hs[j]<- sum((res1$BetaHat - beta)^2  )
      

	## Using Bayesian sqrt Lasso
	bsq <- BSqrt_lasso ( y , x , r = 1 ,delta =  1 ,start_values =  1 , nmc = 8000 , nburn = 1000, thin = 2)$BetaHat
	bsq_ind <- class_vector(bsq)
	mp_Bsql[j] = sum((bsq_ind!=index))
	mse_Bsql[j] = sum((bsq - beta)^2)

	## Using DL sqrt
	dl <- DL_square( y , x , a = 1/p  , nmc = 6000 , nburn = 1000 , thin = 2)$BetaHat
	dl_ind <- class_vector(dl)
	mp_DL_sqrt[j] = sum((dl_ind!=index))
	mse_DL_sqrt[j] = sum( (dl - beta)^2 )
    ptm[j+1] <- proc.time()["elapsed"] - sum(ptm)
	
	
	
      lasso.mod = glmnet(x,y,alpha = 1)
      cv.out = cv.glmnet(x,y,alpha = 1, nfolds=10)
      bestlam = cv.out$lambda.min
      lasso.coef = predict(lasso.mod,type ="coefficients",s=bestlam)
      lasso_cv_est = lasso.coef[-1]
      lasso_ind <- ((lasso_cv_est!=0))
      mp_lasso[j] <- sum((lasso_ind!=index))
	  mse_lasso[j] = sum( ( lasso_cv_est - beta )^2 )

	    if(!exists("pb")) pb <- tkProgressBar("Parallel task", min=1, max=reps)
    setTkProgressBar(pb, j)
    }
    prop_true_hs <- sum((mp_hs==0))/reps
    prop_true_lasso <- sum((mp_lasso==0))/reps
    prop_true_Bsqrl <- sum((mp_Bsql==0))/reps
    prop_true_DLsqrt <- sum((mp_DL_sqrt==0))/reps
	MSE_hs <- sum(mse_hs)/reps
	MSE_Bsqrl <- sum(mse_Bsql )/reps
	MSE_DLsqrt <- sum( mse_DL_sqrt )/reps
	MSE_lasso <- sum( mse_lasso ) / reps

     t <- paste("Iteration" , i , sum(ptm[-1]))
    #info <- sprintf("%d%% done", round((i/niter)*niter))
     #setWinProgressBar(pb, i/(niter)*niter, label=info)
    #if(!exists("pb")) pb <- tkProgressBar("Parallel task", min=1, max=n)
   # setTkProgressBar(pb, i)
    cbind(n_inf,prop_true_hs,prop_true_lasso ,prop_true_Bsqrl ,  prop_true_DLsqrt , MSE_hs , MSE_lasso ,  MSE_Bsqrl , MSE_DLsqrt , t)
}

 close(pb)
stopCluster(cl)
registerDoSEQ()


####Trying to come up with a plot


str(res ,0 )
#mp <- lapply(res , function(x) drop(x[,-6]))
dres <- do.call( rbind , res)
#dres <- apply(dres , 2 , as.numeric)
str(dres)

dres <- as.data.frame(dres)

res_mp <- dres[,c(1 ,grep("prop_" , names(dres)))]
res_mse <- dres[,c(1 ,grep("MSE_" , names(dres)))]

library(reshape2)
library(ggplot2)

### Plot for misclassification Probability
data_mp <- melt(res_mp , id = "n_inf")

fig_mp <- ggplot(data_mp , aes(x = n_inf , y = value , by = variable )) + facet_grid(variable~.) +
 geom_point() + theme_bw() 
 
### Plot for MSE
x11()
data_mse <- melt(res_mse , id = "n_inf")

fig_mse <- ggplot(data_mse , aes(x = n_inf , y = value , by = variable )) + facet_grid(variable~.) +
 geom_point() + theme_bw() 

## Recording the execution time
time <- dres[,'t']

#time <- do.call(rbind , time)
