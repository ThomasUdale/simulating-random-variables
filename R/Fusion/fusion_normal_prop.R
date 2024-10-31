library(mcreplicate)
library(parallel)
library("ggplot2")

C = 5.0
temp = 1.0

x0 <- seq(-5, 5, by = 0.01)

f <- function(x) {
  exp(-(x)^4/2.0)/2.1558
}

fc <- function(x) {
  exp(-(x)^4/(2.0*C))/(2.1558)^(1/C)
}

qc_x <- function(x,y) {
  exp(-((y-x)^2)/(2.0*temp))
}

qc_unif <- function(x,y) {
  if (abs(x-y)<0.5) {
    1
  }
  else {
    0
  }
}

pc_x <- function(x,y) {
  fc(y)/fc(x)
}


g <- function(x,y) {
  prod(fc(x))*prod(sapply(x,pc_x,y=y))
}

h <- function(x,y) {
  prod(fc(x))*exp(-(C*(y-mean(x))^2)/(2.0*temp))
}

r <- function(x,y) {
  g(x,y)/h(x,y)
}

samplefc <- function() {
  while (TRUE) {
    x <- rnorm(1)
    if (runif(1)<(fc(x)/exp(-x^2/2+C/8))) {
      return(x)
    }
  }
}


# Full algo

sample_fusion <- function(){

  while (TRUE) {
    x <- replicate(C,samplefc())
    x_bar <- mean(x)
    y <- rnorm(1,x_bar,sqrt(temp/C))
    alpha = prod(fc(y)/fc(x))
    
    h_star = exp((C*(y-x_bar)^2)/(2.0*temp))
    
    sigma_2 = var(x)*(C-1)/C
    rho = exp(-0.5*C*sigma_2/temp)
    
    if (runif(1) < alpha*h_star) {
      return(y)
    }
  }
}

sample = pbapply::pbreplicate(10000,sample_fusion())


curve(f,-3,3,ylim=c(0,0.8))
# hist(sample,add=T,freq=F,breaks=30)
dens = density(sample)
lines(dens,col='red')


# gs = c()
# 
# for (i in x0) {
#   gs = append(gs,pc_x(0,i))
# }





