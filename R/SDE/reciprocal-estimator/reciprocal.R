library(ggplot2)

theta = 4
lam <- exp(-theta)
u_n = 1e32

unbiased <- function(){
  k <- rgeom(1,1-exp(-1))
  if(k==0){
    return(1)
  }
  x <- 1
  for (i in 1:k) {
    x <-x + prod(rexp(i,1/theta))*exp(i)/factorial(i)
  }
  return(x)
}

sim_ber <- function(l){
  as.integer(unbiased()>=l)
}

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

tdash3 <- function(){
  while(TRUE){
    x <- tdash()
    if (runif(1)<1/x){
      return(x)
    }
  }
}

linear_bernoulli <- function(l) {
  C <- l^2/(2*u_n)
  if (C<=1){
    return((runif(1)<C && sim_ber(l)))
  }
  eps <- 1/2
  gamma <- 1/2
  k <- 2.3/(gamma*eps)
  i <- 1
  R <- TRUE
  
  while((i!=0) && R){
    while ((i>0) && (i<k)) {
      B <- sim_ber(l)
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

tsim <- function(vm){
  while(TRUE){
    k <- tdash()
    a <- linear_bernoulli(k)
    if (a){
      return(k==1)
    }
  }
}

# tsim_out <- pbmcapply::pbmclapply((1:1e3),tsim,mc.cores=10)
# 
# plot = ggplot()+
#   xlim(0,1)+
#   ylim(0,1)+
#   geom_vline(xintercept=1/theta)+
#   geom_vline(xintercept=mean(bias_est),col='red')+
#   geom_vline(xintercept=mean(unlist(tsim_out)),col='green')
# 
# print(plot)


