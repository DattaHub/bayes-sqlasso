#setwd("C:/Users/pantho/Desktop/Mohamed/Horseshoe")
#.libPaths("C:/Users/pantho/Desktop/Mohamed/Mylib/Mylib")
setwd("C:/Users/admin/Desktop/Horseshoe")
.libPaths("C:/Users/admin/Documents/R/win-library/Mylib")

library(glmnet)
library(monomvn)
library(bayesm)
library(plyr)
library(lars)
library(horseshoe)
library("parallel")
library("foreach")
library("doParallel")
library(statmod)
library(MASS)
library(GIGrvg)
library(scalreg)
library(tcltk)
source("DL_sqrt.R")
source("Bayesian_sqrtLasso.R")
#setwd("C:/Users/admin/Desktop/Horseshoe")

no_cores <- detectCores()-1 
RNGkind("L'Ecuyer-CMRG")
cl <- makeCluster(no_cores)
registerDoParallel(cl)
clusterSetRNGStream(cl, 123)
#clusterEvalQ(cl, .libPaths("C:/Users/pantho/Desktop/Mohamed/Mylib/Mylib"))
#clusterEvalQ(cl , setwd("C:/Users/pantho/Desktop/Mohamed/Horseshoe"))
clusterEvalQ(cl, .libPaths("C:/Users/admin/Documents/R/win-library/Mylib"))
clusterEvalQ(cl , setwd("C:/Users/admin/Desktop/Horseshoe"))
clusterEvalQ(cl , source("DL_sqrt.R"))
clusterEvalQ(cl , source("Bayesian_sqrtLasso.R"))

niter = 100

# Uncomment this if you want a progress bar
# pb <- winProgressBar(title="Progress bar", label="0% done", min=0, max=100, initial=0)

n=100;p=100
#p1 = 50
#beta <- c(c(rep(7,14),rep(5,14),rep(4,11),rep(3,11)),rep(0,p-p1))
#index <- c(rep(1,p1),rep(0,p-p1))

#Design matrix does not change
 
res<- foreach(i = 1:niter, .verbose = TRUE  ,.packages = c("bayesm","horseshoe","glmnet" , "MASS" , "statmod",
							    "GIGrvg" , "tcltk") )%dopar%{
  
    if(!exists("pb")) pb <- tkProgressBar("Parallel task", min=1, max=niter)
    setTkProgressBar(pb, i)
    #sigma <- rwishart(p,diag(p))
    
    #x=matrix(rnorm(n*p),n,p) # generate X matrices
    #R <- chol(sigma$W)
    #x <- x%*%R
    #fx <- x%*% beta
    #A <- t(x[,-(1:p1)])%*%x[,(1:p1)]%*%solve((t(x[,(1:p1)])%*%x[,(1:p1)]))%*%rep(1,p1)
    #n_inf <- 1-norm(A,'i')
	prop = tail(head(seq(0 , 1 , by = .1) , -1) , -1)
	sparsity = prop
    matbeta = matrix(0 , nrow = p , ncol = 9)
    X = matrix(rnorm(n*p , 0 , sqrt(2)) , ncol = p) 
    for(i in 1:9)  matbeta[,i] = c(rep(5 , prop[i]*p) , rep(0 , ceiling(p-prop[i]*p)) ) 

    L <- as.list(1:9)
	build <- function(X , beta) {
	return(list(beta0 = beta , y = (crossprod(t(X) , beta) + rnorm(n , 0 , sqrt(5)))))
						}

	L <- lapply(L  , function(x) build( X , matbeta[,x] ))
    reps = length(L)
    mp_hs = rep(0,reps)
    mp_lasso = rep(0, reps)
    mp_Bsql= rep(0, reps)
    mp_DL_sqrt = rep(0, reps)
	mse_hs = rep(0,reps)
    mse_lasso = rep(0, reps)
    mse_Bsql= rep(0, reps)
    mse_DL_sqrt = rep(0, reps)
	mse_hs_dec = rep(0,reps)
    #mse_lasso = rep(0, reps)
    mse_Bsql_dec= rep(0, reps)
    mse_DL_sqrt_dec = rep(0, reps)
	
	#mp_hs_median = rep(0,reps)
    #mp_lasso = rep(0, reps)
    mp_Bsql_median= rep(0, reps)
    mp_DL_sqrt_median = rep(0, reps)
	mse_hs_median = rep(0,reps)
    #mse_lasso_median = rep(0, reps)
    mse_Bsql_median= rep(0, reps)
    mse_DL_sqrt_median = rep(0, reps)
	mse_hs_median_dec = rep(0,reps)
    #mse_lasso = rep(0, reps)
    mse_Bsql_median_dec= rep(0, reps)
    mse_DL_sqrt_median_dec = rep(0, reps)
    #bsq <- matrix(0 , reps , p)
    ptm <- c(proc.time()["elapsed"] , rep(0 ,reps))
	



    for (j in 1:reps){ 
     # eps=sqrt(5)*rnorm(n)
      y=drop(L[[j]]$y)
	  beta <- L[[j]]$beta0
	  index <-beta
	  index[(index!=0)]=1
     # data = data.frame(y,x)
      ## Using horseshoe
      res1 <- horseshoe(y, X, method.tau = "truncatedCauchy",
                       method.sigma = "Jeffreys", burn = 1000, nmc = 4000)
      hs_ind <- HS.var.select(res1, y, method = "intervals")
      mp_hs[j] <- sum((hs_ind!=index))
	  mse_hs[j]<- sum((res1$BetaHat - beta)^2  )
	  mse_hs_median[j] <- sum((res1$BetaMedian - beta)^2  )
	  #hs_betahat <- res1$BetaHat*hs_ind
	  mse_hs_median_dec[j] <- sum((res1$BetaMedian*hs_ind - beta)^2  )
	  mse_hs_dec[j] <- sum((res1$BetaHat*hs_ind - beta)^2  )
      

	## Using Bayesian sqrt Lasso
	bsq <- BSqrt_lasso ( y , X , r = 1 ,delta =  1 ,start_values =  1 , nmc = 8000 , nburn = 1000, thin = 2)
	bsq_median <- bsq$BetaMedian
	bsq <- bsq$BetaHat
	ifelse(length(bsq)!=p ,print("Dimension error") ,print("OK"))
	ifelse(length(bsq_median)!=p ,print("Dimension error") , print("OK") )
	bsq_ind <- class_vector(bsq)
	bsq_median_ind <- class_vector(bsq_median)
	mp_Bsql[j] = sum((bsq_ind!=index))
	mse_Bsql[j] = sum((bsq - beta)^2)
	mp_Bsql_median[j] = sum((bsq_median_ind!=index))
	mse_Bsql_median[j] = sum((bsq_median - beta)^2)
	mse_Bsql_dec[j] = sum((bsq*bsq_ind - beta)^2)
	mse_Bsql_median_dec[j] = sum((bsq_median*bsq_median_ind - beta)^2)
	
	## Using DL sqrt
	dl <- DL_square( y , X , a = 1/p  , nmc = 8000 , nburn = 1000 , thin = 2)
	dl_median <- dl$BetaMedian
	dl <- dl$BetaHat
	ifelse(length(dl)!=p ,print("Dimension error") , print("OK"))
	ifelse(length(dl_median)!=p ,print("Dimension error") , print("OK"))
	dl_ind <- class_vector(dl)
	dl_median_ind <- class_vector(dl_median)
	mp_DL_sqrt[j] = sum((dl_ind!=index))
	mp_DL_sqrt_median[j] = sum((dl_median_ind!=index))
	mse_DL_sqrt[j] = sum( (dl - beta)^2 )
	mse_DL_sqrt_median[j] = sum( (dl_median - beta)^2 )
	mse_DL_sqrt_dec[j] = sum((dl*dl_ind - beta)^2)
	mse_DL_sqrt_median_dec[j] = sum((dl_median*dl_median_ind - beta)^2)
	ptm[j+1] <- proc.time()["elapsed"] - sum(ptm)
	
      lasso.mod = glmnet(x,y,alpha = 1)
      cv.out = cv.glmnet(x,y,alpha = 1, nfolds=10)
      bestlam = cv.out$lambda.min
      lasso.coef = predict(lasso.mod,type ="coefficients",s=bestlam)
      lasso_cv_est = lasso.coef[-1]
      lasso_ind <- ((lasso_cv_est!=0))
      mp_lasso[j] <- sum((lasso_ind!=index))
	  mse_lasso[j] = sum( ( lasso_cv_est - beta )^2 )
  
    }
	
    
    # Missclas_hs <- p - sum(mp_hs)/reps
    # Missclas_lasso <- p - sum(mp_lasso)/reps
    # Missclas_Bsqrl <- p - sum(mp_Bsql)/reps
    # Missclas_DLsqrt <- p - sum(mp_DL_sqrt)/reps
    # prop_true_hs <- sum((mp_hs==0))/reps
    # prop_true_lasso <- sum((mp_lasso==0))/reps
    # prop_true_Bsqrl <- sum((mp_Bsql==0))/reps
    # prop_true_DLsqrt <- sum((mp_DL_sqrt==0))/reps
	# MSE_hs <- sum(mse_hs)/reps
	# MSE_Bsqrl <- sum(mse_Bsql )/reps
	# MSE_DLsqrt <- sum( mse_DL_sqrt )/reps
	# MSE_lasso <- sum( mse_lasso ) / reps
	
	 
	
    t <-  ptm[-1]
    #if(!exists("pb")) pb <- tkProgressBar("Parallel task", min=1, max=niter)
    #setTkProgressBar(pb, i)

    cbind(sparsity,mp_hs , mp_lasso , mp_Bsql, mp_DL_sqrt , mse_hs , mse_lasso , mse_Bsql , mse_DL_sqrt , mse_hs_dec , mse_Bsql_dec , mse_DL_sqrt_dec ,mp_Bsql_median ,
    mp_DL_sqrt_median ,	mse_hs_median , mse_Bsql_median , mse_DL_sqrt_median , mse_hs_median_dec , mse_Bsql_median_dec , mse_DL_sqrt_median_dec , t)
    
}

 close(pb)
stopCluster(cl)
registerDoSEQ()
stopImplicitCluster()
m <- paste("replicate" , 1:100)
res_array <- simplify2array(res)
Sim = list(Results = res_array , p = p , n = n , Beta = matbeta)
save(Sim, file = "Sparsity_simulations_n=p.RDATA")

dim(res_array)
dimnames(res_array)
dimnames(res_array)[c(1,3)] <- list( n  = paste( prop) , m)
out <- apply(res_array , c(1,2) , mean)
str(out,1)
dimnames(out)
###
###Now we need the MSE MP and each in a data frame 
###

#mse_dec <- grep("*_dec" , dimnames(out)[[2]])
median <-  grep("*_median_*" , dimnames(out)[[2]])
median <- dimnames(out)[[2]][median]
median <- c("mse_lasso","mp_lasso","mp_hs",median)
mp <- grep("mp_*" , median)
t<-median[mp]
median <- median[-mp]
mp<-t
median_dec <- c("mse_lasso" , median[grep("*_dec" , median)])
median_dec
median_mse <- median[-grep("*_dec" , median)]
median_mse

out_mse <- as.data.frame(out[,median_mse])
out_dec <- as.data.frame(out[,median_dec])
out_mp <- as.data.frame(out[,mp])

out_mp <- out_mp/p
names(out_mp)<- c("Lasso" , "Horseshoe" , "Sqrt-Lasso" , "Sqrt-DL")
names(out_dec)<- c("Lasso" , "Horseshoe" , "Sqrt-Lasso" , "Sqrt-DL")
names(out_mse)<- c("Lasso" , "Horseshoe" , "Sqrt-Lasso" , "Sqrt-DL")





library(reshape2)
#out_mp_median <- melt(as.matrix(out_mp_median))
out_mp <- melt(as.matrix(out_mp))
names(out_mp) <- c("Sparsity" , "Method" ,"Proportion")

library(ggplot2)

 misclas <- ggplot(out_mp , aes(x = Sparsity , y = Proportion , group = Method , color = Method )) + geom_line()  + theme_bw() + geom_point(aes( shape  =Method ), size = 4)
 misclas + scale_x_continuous(breaks = prop)
#misclas


out_mse <- melt(as.matrix(out_mse))
names(out_mse) <- c("Sparsity" , "Method" ,"MSE")

#library(ggplot2)
x11()
 mse <- ggplot(out_mse , aes(x = Sparsity , y = MSE , group = Method , color = Method )) + geom_line()  + theme_bw() + geom_point(aes( shape  =Method ), size = 4)
 mse + scale_x_continuous(breaks = prop)
#mse
 
 #MSE after dec
 out_dec <- melt(as.matrix(out_dec))
names(out_dec) <- c("Sparsity" , "Method" ,"MSE")

#library(ggplot2)
x11()
 dec <- ggplot(out_dec , aes(x = Sparsity , y = MSE , group = Method , color = Method )) + geom_line()  + theme_bw() + geom_point(aes( shape  =Method ), size = 4)
 dec + scale_x_continuous(breaks = prop)
#dec




####Trying to come up with a plot

res <- Out$Output
str(res ,0 )
#mp <- lapply(res , function(x) drop(x[,-6]))
dres <- do.call( rbind , res)
#dres <- apply(dres , 2 , as.numeric)
str(dres)

dres <- as.data.frame(dres)

res_mp <- dres[,c(1 ,grep("prop_" , names(dres)))]
res_mse <- dres[,c(1 ,grep("MSE_" , names(dres)))]
res_Miss <- dres[, c(1 , grep("Missclas_" , names(dres)))]

library(reshape2)
library(ggplot2)

#pdf(  width = 9 ,  height = 13 , "a4" , title = "Beta with different signal strength 7 , 5 , 4 and 3" )

### Plot for misclassification Probability (True Model)
data_mp <- melt(res_mp , id = "n_inf")

fig_mp <- ggplot(data_mp , aes(x = n_inf , y = value , by = variable )) + facet_grid(variable~.) +
 geom_point() + theme_bw() 
 fig_mp
 
### Plot for misclassification Probability (Average number of missclassified coefficients)
data_Miss <- melt(res_Miss , id = "n_inf")
x11()
fig_Miss <- ggplot(data_Miss , aes(x = n_inf , y = value , by = variable )) + facet_grid(variable~.) +
 geom_point() + theme_bw() 
 fig_Miss
 
 
### Plot for MSE
x11()
data_mse <- melt(res_mse , id = "n_inf")

fig_mse <- ggplot(data_mse , aes(x = n_inf , y = value , by = variable )) + facet_grid(variable~.) +
 geom_point() + theme_bw() 
fig_mse


dev.off()
## Recording the execution time
time <- dres[,'t']

#time <- do.call(rbind , time)

Out <- list(Output = dres , Beta0 = beta)

save(Out , file = "mp_n100_p60_q_50_1.Rdata")

