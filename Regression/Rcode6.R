# setwd("~/0 - Teaching/Computational Statistics/Rcode")

# An example of bootstrap inconsistency

###########################################################
# Example taking from Andrews (2000):
#   https://www.jstor.org/stable/2999432?seq=1#page_scan_tab_contents

# We consider iid data X_1,...,X_n from N(mu,1) distribution, 
#   where we know mu >= 0.
# 
# The MLE for mu: hat mu_n = max(\bar X_n,0).
#   sqrt(n)*(hat mu_n - mu) is asymptotically distributed as:
#      Z if mu>0
#      max(Z,0) if mu=0
#   where Z ~ N(0,1)
# The bootstrap is consistent if mu>0 and inconsistent if mu=0.
#   (Inconsistency can occur if parameter is on the boundary 
#    of the parameter space.)
#
# Note: In this case we know the limiting distribution, 
#       so there is no real need to do the bootstrap. 
#       But since we know the limiting distribution, we can use 
#       this as an example to see when the bootstrap works and 
#       doesn't work. 

library(boot)

# function that samples data set of size n and computes MLE
sim.mle <- function(n, mu){
  data <- rnorm(n, mean=mu, sd=1)
  return(max(mean(data),0))
}

# simulation settings:   
# (run the code once with mu=0.2. Then the bootstrap is consistent.
#  run the code once with mu=0. Then the boostrap is inconsistent.)

#mu <- 0.2
mu <- 0
n <- c(10,50,250,1250)
reps <- 5000

#### First we simulate the limiting distribution:

sim.result <- matrix(numeric(4*reps),ncol=4)
# ith column will contain mles for sample size n[i]

set.seed(123)

par(mfrow=c(2,2))
xlim=c(0,1)
for(i in 1:4){
  sim.result[,i] <- replicate( reps, sim.mle(n[i],mu) )
  hist(sim.result[,i], freq=F, breaks=20, xlim=xlim, main=paste("Simulated distr., n=",n[i],", mu=", mu, sep=""))
  grid <- seq(from=0,to=1,length=100)
  if (mu>0){
    lines(grid, dnorm(grid, mean=mu, sd=1/sqrt(n[i])),col="red")
  }
}
# red curve is density of asymptotic normal distrbution of hat mu_n

# rescaled distribution: sqrt(n)*(hat mu_n - mu)
sim.scaled <- matrix(numeric(4*reps), ncol=4)
for (i in 1:4){
  sim.scaled[,i] <- sqrt(n[i])*(sim.result[,i] - mu)
}

#jpeg(paste("Bootstrap-Simulation-mu",10*mu,".jpg", sep=""))
par(mfrow=c(2,2))
xlim=c(-5,5)
for(i in 1:4){
  hist(sim.scaled[,i], freq=F, breaks=20,xlim=xlim, 
         main=paste("Rescaled sim. distr., n=",n[i],", mu=",mu,sep=""))
  grid <- seq(from=-5,to=5,length=100)
  if (mu>0){
     lines(grid, dnorm(grid), col="red")
  }else{
    lines(grid[grid>=0], dnorm(grid[grid>=0],sd=1), col="red")
  }
  #legend("topleft", c("Limiting distr."), lty=c(1), col=c("red"))
}
#dev.off()


#### Now try bootstrap:

# function that computes median of data[index]
mle.fn <- function(data, index){
  return(max(mean(data[index]),0))
}

# settings for bootstrap:
set.seed(321)
R <- 5000
#n <- 250
n <- 10000
xlim <- c(mu-8/sqrt(n),mu+8/sqrt(n))

boot.result <- matrix(numeric(4*R),ncol=4)
# ith column will contain bootstrapped estimators for data set i

mles <- numeric(4)
# ith entry will contain mle based on data set i

par(mfrow=c(2,2))
for(i in 1:4){
  data <- rnorm(n,mean=mu)
  mles[i] <- max(mean(data),0)
  boot.result[,i] <- boot(data, mle.fn, R)$t
  hist(boot.result[,i], freq=F, breaks=20, xlim=xlim, 
       main=paste("Bootstrap distr., data set ",i),
       xlab=paste("n=",n,", mu=",mu,", sample mean=",
                  round(mean(data),4),sep=""))
  abline(v=mles[i],col="red",lwd=2)
}

# rescaled distribution: sqrt(n)*(hat mu_n* - hat mu_n)
boot.scaled <- matrix(numeric(4*R), ncol=4)
for (i in 1:4){
  boot.scaled[,i] <- sqrt(n)*(boot.result[,i] - mles[i])
}

#jpeg(paste("Bootstrap-mu",10*mu,".jpg", sep=""))
par(mfrow=c(2,2))
xlim=c(-5,5)
ylim=c(0,.7)
for(i in 1:4){
  hist(boot.scaled[,i], freq=F, breaks=20, xlim=xlim, ylim=ylim, 
       main=paste("Rescaled bootstr. distr, data set", i),
       xlab=paste("sd=",round(sd(boot.scaled[,i]),2)))
  grid <- seq(from=-5,to=5,length=1000)
  if (mu>0){
    lines(grid, dnorm(grid, 0, sd=1), col="red")
  }else{
    lines(grid[grid>=0], dnorm(grid[grid>=0],sd=1), col="red")
  }  
  #legend("topleft", c("Limiting distr."), lty=c(1), col=c("red"))
}
#dev.off()

# look at CDFs
#jpeg(paste("Bootstrap-ecdf",10*mu,".jpeg",sep=""))
par(mfrow=c(2,2))
for (i in 1:4){
  plot(ecdf(boot.scaled[,i]), 
       main=paste("Rescaled bootstr. distr., data set",i),
       verticals = TRUE, do.points = FALSE,lwd=3)
  gridl <- seq(from=-5,to=0,length=100)
  gridr <- seq(from=0,to=5,length=100)
  if (mu>0){
    lines(grid, pnorm(grid), col="red",lwd=2)
  }else{
    lines(gridr, pnorm(gridr), col="red",lwd=2)
    lines(gridl, rep(0,length(gridl)), col="red",lwd=2)
  }
  #legend("topleft", c("Emp. CDF bootstrap", "Limiting distr."), lty=c(1,1), col=c(1,"red"))
}
#dev.off()



##############################################################
# Estimating the Accuracy of a Linear Regression Model

# Simulate from model with heteroscedastic errors 
# (non-constant variance)

library(boot)
set.seed(999)

n <- 100

generate.data <- function(){
   x <- runif(n,0,100)
   eps <- rnorm(n,0,abs(x-50))
   y <- 50 + 0.3*x + eps
   return(data.frame(cbind(x,y)))
}

dat <- generate.data()

par(mfrow=c(1,1))
plot(y~x,data=dat)
abline(50,0.3)

fit1 <- lm(y~x,data=dat)
summary(fit1)
confint(fit1)

# We cannot trust the standard error estimates, and hence
#   we cannot trust the p-values and confidence intervals
# Since we generated the data ourselves, we can approximate
#   the standard errors by simulation

simul <- function(){
  data.temp <- generate.data()
  return(coef(lm(y~x,data=data.temp)))
}

res.coef <- replicate(1000, simul())
(se.int.sim <- sd(res.coef[1,]))
(se.slope.sim <- sd(res.coef[2,]))

# This is a situation where the bootstrap should work.
# We resample pairs of (x,y)
# Why is it wrong to resample only the y's?

boot.fn=function(data,index){
  return(coef(lm(y~x,data=data,subset=index)))
}

boot(dat, boot.fn, 5000)

# compare to simulation values
se.int.sim
se.slope.sim

# compare to output of lm()
summary(fit1)$coef


