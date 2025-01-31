library(ggplot2)
library(layeredBB)

tt <- pi

s0 = 0

l = -0.5
r = 9/8
phi <- function(x) {
  (sin(x)^2 + cos(x))/2
}

sim_end <- function(t) {
  while (T) {
    u <- rnorm(1,0,sqrt(t))
    if (runif(1) < exp(-cos(u)-1)) {
      return(u)
    }
  }
}

ea1 <- function(t){
  while (T){
    y <- sim_end(t)
    k <- rpois(1,r*t)
    if (k==0){
      return(y)
    }
    skel_t <- sort(runif(k,0,t))
    skel_u <- runif(k,0,r)
    bb <- Brownian_bridge_path_sampler(
      x=s0,
      y=y,
      s=0,
      t=t,
      times=skel_t
    )
    if (all((phi(bb$simulated_path[1,])-l)<skel_u)){
      return(y)
    }
  }
}

euler_scheme <- function() {
  N <- 100
  h <- tt/N
  z <- s0
  for (i in 1:N){
    z <- z + sin(z)*h+rnorm(1,0,sqrt(h))
  }
  return(z)
}

alpha <- 0.5
delta <- 1
lambda <- function(s) {
  delta * s^(alpha-1)
}


gamma <- function(x) {
  1
}
gamma_p <- function(x) {
  0
}
gamma_pp <- function(x) {
  0
}
b <- function(x) {
  sin(x)
}
b_p <- function(x){
  cos(x)
}

rho <- function(theta,x,y,u){
  1 + 1/lambda(u)*(1/(gamma(x)*u)*((gamma(y)-gamma(x))/2*(((y-x-u*b(x))^2)/(gamma(x)*u)-1)+(b(y)-b(x)-gamma_p(y))*(y-x-u*b(x))) + gamma_pp(y)/2 - b_p(y))
}



cts_is <- function() {
  tau <- 0
  w <- 1
  x <- s0
  k <- 0
  theta <- c(b(x),gamma(x))
  while(TRUE){
    event_time <- (-alpha*log(runif(1))/delta)^(1/alpha)
    new_tau <- tau + event_time
    if (new_tau > tt) {
      x <- rnorm(1, x + (tt-tau)*theta[1],sqrt(tt-tau)*sqrt(theta[2]))
      return(c(x,w))
    } else {
      tau <- new_tau
      new_x <- rnorm(1, x + event_time*theta[1],sqrt(event_time)*sqrt(theta[2]))
      w <- w*rho(theta,x,new_x,event_time)
      x <- new_x
      theta <- c(b(x),gamma(x))
      k <- k+1
    }
  }
}

cts_is_adapt <- function() {
  tau <- 0
  w <- 1
  x <- s0
  k <- 0
  theta <- c(b(x),gamma(x))
  while(TRUE){
    event_time <- (-alpha*log(runif(1))/(delta*(abs(theta[1]+1)/(abs(theta[1])+2))))^(1/alpha)
    new_tau <- tau + event_time
    if (new_tau > tt) {
      x <- rnorm(1, x + (tt-tau)*theta[1],sqrt(tt-tau)*sqrt(theta[2]))
      return(c(x,w))
    } else {
      tau <- new_tau
      new_x <- rnorm(1, x + event_time*theta[1],sqrt(event_time)*sqrt(theta[2]))
      w <- w*rho(theta,x,new_x,event_time)
      x <- new_x
      theta <- c(b(x),gamma(x))
      k <- k+1
    }
  }
}

# ea_out <- pbapply::pbreplicate(1e4,ea1(tt))
# print(summary(ea_out))
# print(summary(ea_out^2))
# eu_out <- pbapply::pbreplicate(1e4,euler_scheme())
# print(summary(eu_out))
# print(summary(eu_out^2))
cts_out <- pbapply::pbreplicate(1e4,cts_is())
print(summary(cts_out[1,]*cts_out[2,]/mean(cts_out[2,])))
print(summary(cts_out[1,]^2*cts_out[2,]/mean(cts_out[2,])))
cts_out_adapt <- pbapply::pbreplicate(1e4,cts_is_adapt())
print(summary(cts_out_adapt[1,]*cts_out_adapt[2,]/mean(cts_out_adapt[2,])))
print(summary(cts_out_adapt[1,]^2*cts_out_adapt[2,]/mean(cts_out_adapt[2,])))



# plot <- ggplot()+
#   xlim(-5,5)+
#   geom_density(aes(x=ea_out))+
#   geom_density(aes(x=cts_out[1,],weight=cts_out[2,]))
# 
# print(plot)