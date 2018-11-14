## Ex 1
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

setwd("C:/Users/Jyotishka/OneDrive/Documents/R/sqlasso")
source("bayes.sqrt.lasso.R")

par(mfrow=c(1,1))
theta = c(rep(7,10),rep(0,90))
Y = rnorm(100,theta,1)
ans1 = bayes.sqrt.lasso(Y, method.tau ="fixed", r = 1, delta = 1/2, burn = 1000, nmc = 5000)

## sql.normal.means and eval-sql should be the same. 
# 
par(mfrow=c(1,1))
plot(Y)
lines(ans1$ThetaHat)
points(apply(ans1$ThetaSave,2,median),col="red",pch=15)
# 
# Ex 2
theta = c(rep(7,10),rep(3,10),rep(0,80))
Y = rnorm(100,theta,0.1)s
ans2 = bayes.sqrt.lasso(Y, method.tau ="fixed", burn = 1000, nmc = 5000)
# 
plot(Y)
lines(ans2$ThetaHat)
points(apply(ans2$ThetaSave,2,Mode),col="red",pch=15)

## Posterior for a range of deterministic values 
 
y <- seq(-5, 5, 0.05)

ans3 <- bayes.sqrt.lasso(Y=y, method.tau ="fixed", burn = 1000, nmc = 5000)
plot(y, apply(ans3$ThetaSave,2,mean), type="l",col="blue",ylab = expression(hat(theta)))
lines(y, apply(ans3$ThetaSave,2,median), type="l",col="red")
lines(y, HS.post.mean(y, tau = 0.5, Sigma2 = 1),col="magenta")
lines(y,y)
legend("topleft",c("SQL-Mean","SQL-Median","Horseshoe Mean","Y"),
       col=c("blue","red","magenta","black"),lty=1)
dev.copy2pdf(file="shrinkage_profile.pdf")

# 
# 
# ##############################
# ########## Fancy Plots #######
# ##############################
# 

set.seed(198)
theta_1 = c(rep(7, 10), rep(0, 90))
Y_1 = rnorm(100, theta_1, 1)
Y_2 = rnorm(100, theta_1, 0.5)

theta_2 = c(rep(7, 10), rep(3,10), rep(0, 80))
Y_3 = rnorm(100, theta_2, 1)
Y_4 = rnorm(100, theta_2, 0.5)

prob_set = list(list(name='theta_1, sigma=1', obs=Y_1),
                list(name='theta_1, sigma=2', obs=Y_2),
                list(name='theta_2, sigma=1', obs=Y_3),
                list(name='theta_3, sigma=2', obs=Y_4))

theta.data = lapply(prob_set,
                    function(set)
                      list('settings'=set$name,
                           'samples'=bayes.sqrt.lasso(set$obs,burn = 1000,
                                              nmc=5000)))

theta.data = lapply(theta.data,
                    function(sample.data)
                      cbind(melt(sample.data$samples$ThetaSave),
                            'settings'=sample.data$settings,
                            'var'='theta'))

theta.data = Reduce(rbind, theta.data)

colnames(theta.data)[1:2] = c('sample', 'component')

obs.data = lapply(prob_set,
                  function(set)
                    data.frame(component=seq_along(set$obs),
                               value=set$obs,
                               sample=NA,
                               settings=set$name,
                               var="obs"))
obs.data = Reduce(rbind, obs.data)

all.data = rbind(obs.data, theta.data)

plot.data = all.data %>%
  group_by(component, settings, var) %>%
  summarise(upper = quantile(value, prob=0.95),
            lower = quantile(value, prob=0.05),
            middle = median(value))

library(ggplot2)
theta.plot = ggplot(plot.data,
                    aes(x=component, y=middle, group=component,
                        colour=var)) + theme_bw() +
  geom_pointrange(aes(ymin=lower, ymax=upper), size=0.2, alpha=0.5) +
  facet_grid(settings ~ ., scales="free_y") + #, labeller = label_both) +
  scale_colour_manual(values=c("#D55E00","#0072B2")) +
  xlab("") + ylab("")

print(theta.plot)
#setwd("C:/Users/Jyotishka/OneDrive/Documents/R/sqlasso")
ggsave("sparse-means-sql-1.pdf", plot=theta.plot)

