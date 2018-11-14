setwd("C:/Users/admin/Desktop/Horseshoe")
.libPaths("C:/Users/admin/Documents/R/win-library/Mylib")
m <- paste("replicate" , 1:100)
load("Sparsity_simulations_n=p.RDATA")
res_array <- Sim$Results
p<-Sim$p
n<- Sim$n
matbeta <- Sim$Beta
prop = tail(head(seq(0 , 1 , by = .1) , -1) , -1)
sparsity = prop

dim(res_array)
dimnames(res_array)
dimnames(res_array)[c(1,3)] <- list( n  = paste(prop) , m)
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

out_mse[,"id"] <- rep("MSE",9)
out_dec[,"id"] <- rep("Decision_MSE",9)

temp <- rbind(out_mse, out_dec)
temp <- melt(temp ,id.vars = "id")
temp[,"Sparsity"] <- rep(prop , 8)
names(temp)[c(2,3)] <- c("Method" , "MSE")
mmse <- ggplot(temp , aes(x = Sparsity , y = MSE , group = Method , color = Method )) + geom_line()  + theme_bw()+scale_x_continuous(breaks = prop)
 mmse+ geom_point(aes( shape  =Method ), size = 4) + facet_wrap(~id , ncol = 2)
 


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

