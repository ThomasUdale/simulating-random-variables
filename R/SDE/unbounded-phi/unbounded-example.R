library(ggplot2)
# devtools::load_all("~/Documents/simulating-random-variables/R-packages/cts.smc")

f <- function(x) {
  exp(-(x)^4/2.0)/2.1558
}

phi <- function(x) {
  0.5*(4*x^6 - 6*x^2)
}

l = -sqrt(2)

tdash <- function(){
  while(TRUE){
    u <- runif(1)
    v <- runif(1)
    x <- floor(1/u)
    if (2*v*x <= x+1){
      return(x)
    }
  }
}

est_exp <- function(x,y,s,t,n){
  
}

sim_ber <- function(x,y,s,t,n,l){
  as.integer(est_exp(x,y,s,t,n)>=l)
}

linear_bernoulli <- function(x,y,s,t,n,l) {
  C <- l^2/(2*u_n)
  if (C<=1){
    return((runif(1)<C && sim_ber(x,y,s,t,n,l)))
  }
  eps <- 1/2
  gamma <- 1/2
  k <- 2.3/(gamma*eps)
  i <- 1
  R <- TRUE
  
  while((i!=0) && R){
    while ((i>0) && (i<k)) {
      B <- sim_ber(x,y,s,t,n,l)
      G <- rgeom(1,(C-1)/C)+1
      i <- i-1 + (1-B)*G
    }
    if (i>=k){
      R <- runif(1)<((1/(1+gamma*eps))^i)
      C <- C*(1+gamma*eps)
      eps <- eps*(1-gamma)
      k <- k/(1-gamma)
    }
  }
  
  return(as.integer(i==0))
}

sim_bridge_prob<- function(x,y,s,t){
  while(TRUE){
    k <- tdash()
    n <- rgeom(1,1-exp(-1))
    if(linear_bernoulli(x,y,s,t,n,k)){
      return(k==1)
    }
  }
}