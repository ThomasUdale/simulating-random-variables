set.seed(1)
library(progress)
library(layeredBB)
library(glue)
library(mvtnorm)
library(ggplot2)

f <- function(x) {
  dmvnorm(x=x,mean=c(0,0),sigma=diag(2))
}

mu=c(-0.5,0.5)

fc <- function(x,c) {
  dmvnorm(x=x,mean=rep(mu[c],2),sigma=diag(rep(sqrt(2),2)))
}

x <- seq(-3,3,length.out=100)
y <- seq(-3,3,length.out=100)
out <- expand.grid(x,y)

out$z <- apply(out,1,f)
print(ggplot(data=out,aes(x=Var1,y=Var2,z=z)) + geom_contour())

phi <- function(x,c) {
  
}

find_bound <- function(bes_layer,c) {
  
}

ea3 <- function() {
  
}

mcf <- function() {
  
}

bf <- function() {
  
}

po_est <- function(x_in,x_out) {
  
}

pmbf <- function() {
  
}


## run mcf
## run bf
## run pmbf

