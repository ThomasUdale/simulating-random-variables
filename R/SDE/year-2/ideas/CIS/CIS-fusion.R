library(layeredBB)
library(ggplot2)
library(mvtnorm)
library(progress)

f <- function(x) {
  dnorm(x,0,sqrt(0.5))
}

mu_0 = 3
mu = c(0.5*mu_0,-0.5*mu_0)
C = length(mu)

fc <- function(x,c) {
  dnorm(x,mu[c],1)
}

samplefc <- function(c) {
  rnorm(1,mu[c],1)
}

phi <- function(x,c) {
  0.5 * ((x-mu[c])^2-1) + 0.5
}

alpha <- 0.3
delta <- 1
lambda <- function(s) {
  delta * s^(alpha-1)
}
gamma <- function(x) {
  2
}
gamma_p <- function(x) {
  0
}
gamma_pp <- function(x) {
  0
}
b <- function(x,c) {
  -(x-mu[c])
}
b_p <- function(x,c){
  -1
}

rho <- function(x,y,u,c){
  1 + 1/lambda(u)*(1/(gamma(x)*u)*((gamma(y)-gamma(x))/2*(((y-x-u*b(x,c))^2)/(gamma(x)*u)-1)+(b(y,c)-b(x,c)-gamma_p(y))*(y-x-u*b(x,c))) + gamma_pp(y)/2 - b_p(y,c))
}

g_prop <- function(z,y,s,tt,t) {
  rnorm(1,z + (t-s)/(tt-s)*(y-z),sqrt((tt-t)*(t-s)/(tt-s)))
}

g_den <- function(x,z,y,s,tt,t) {
  dnorm(x,z + (t-s)/(tt-s)*(y-z),sqrt((tt-t)*(t-s)/(tt-s)))
}

q_den <- function(x,y,event_time,c) {
  dnorm(y,x + event_time*b(x,c),sqrt(event_time)*sqrt(gamma(x)))
}

cts_is <- function(z,y,tt,c) {
  tau <- 0
  w <- 1
  x <- z
  k <- 0
  while(TRUE){
    event_time <- (-alpha*log(runif(1))/delta)^(1/alpha)
    new_tau <- tau + event_time
    if (new_tau > tt) {
      return(w*q_den(x,y,tt-tau,c))
    } else {
      new_x <- g_prop(x,y,tau,tt,new_tau)
      w <- w*rho(x,new_x,event_time,c)*q_den(x,new_x,event_time,c)/g_den(new_x,x,y,tau,tt,new_tau)
      tau <- new_tau
      x <- new_x
      k <- k+1
    }
  }
}

test <- function() {
  y <- rnorm(1,0,4)
  x1 <- samplefc(1)
  x2 <- samplefc(2)
  w1 <- cts_is(x1,y,1,1)
  w2 <- cts_is(x2,y,1,2)
  
  w <- (w1*w2)/(dnorm(y,0,4))
  return(c(y,w))
}
