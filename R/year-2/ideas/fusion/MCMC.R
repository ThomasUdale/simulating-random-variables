library(ggplot2)
library(mvtnorm)
library(progress)

f <- function(x) {
  dnorm(x,4,sqrt(2/3))
}


mu = c(10,1)
sigma = c(sqrt(2),1)

fc <- function(x,c) {
  dnorm(x,mu[c],sigma[c])
}

samplefc <- function(c) {
  rnorm(1,mu[c],sigma[c])
}

s1 <- rnorm(1e4,mu[1],sigma[1])
s2 <- rnorm(1e4,mu[2],sigma[2])

d1 <- dnorm(s1,mu[1],sigma[1])
d2 <- dnorm(s2,mu[2],sigma[2])

pi <- function(x) {
  prod(dnorm(x,mu,sigma))
}


mcmc <- function(burn,ss){
  pb <- progress_bar$new(
    format = "sampling [:bar] :current/:total (:percent) eta: :eta, elapsed: :elapsed",
    clear = FALSE, total = (ss+burn), width = 80)
  y <- c(0)
  xc <- rnorm(2,mu,sigma)
  mean_xc <- mean(xc)
  sd_xc <- exp(-sd(xc))
  out <- c()
  for (t in 1:(burn+ss)){
    pb$tick()
    new_xc <- rnorm(2,mu,sigma)
    sd_new_xc <- exp(-sd(new_xc))
    mean_new_xc <- mean(new_xc)
    new_y <- rnorm(1,mean_new_xc*sd_new_xc+y*(1-sd_new_xc),(1-sd_new_xc))
    r <- pi(new_y)/pi(y)*dnorm(y,mean_xc*sd_xc+new_y*(1-sd_xc),(1-sd_xc))/dnorm(new_y,mean_new_xc*sd_new_xc+y*(1-sd_new_xc),(1-sd_new_xc))
    # new_y <- rnorm(1,y)
    # r <- pi(new_y)/pi(y)
    if (runif(1)<min(r,1)){
      xc <- new_xc
      mean_xc <- mean_new_xc
      sd_xc <- sd_new_xc
      y <- new_y
    }
    if (t>burn){
      out <- rbind(out,y)
    }
  }
  pb$terminate()
  out
}

out <- mcmc(1e4,1e5)

plot <- ggplot()+
  xlim(0,11)+
  geom_density(aes(x=s1),col='blue',alpha=0.5)+
  geom_density(aes(x=s2),col='blue',alpha=0.5)+
  geom_density(aes(x=out),col='red',alpha=0.4)+
  geom_function(fun=f)

print(plot)