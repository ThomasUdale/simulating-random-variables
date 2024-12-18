library(MASS)
library(ggplot2)
library(progress)
library(parallel)
library(layeredBB)
# devtools::load_all("~/Documents/simulating-random-variables/R-packages/cts.smc")

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

ea1 <- function(x,y,t){
  while (T){
    k <- rpois(1,r*t)
    if (k==0){
      return(1)
    }
    skel_t <- sort(runif(k,0,t))
    skel_u <- runif(k,0,r)
    bb <- Brownian_bridge_path_sampler(
      x=x,
      y=y,
      s=0,
      t=t,
      times=skel_t
    )
    if (all((phi(bb$simulated_path[1,])-l)<skel_u)){
      return(1)
    } else {
      return(0)
    }
  }
}

pint <- function(x) {
  0.5*(0.5*(x-sin(x)*cos(x))+sin(x))-l*x
}

phiint <- function(x,y,s,t) {
  (t-s)/(y-x) * (pint(y)-pint(x))
}

phispline <- function(x,t) {
  n = 10
  while(TRUE){
    y <- sim_end(t)
    test_time <- rexp(1)
    times <- seq(0,t,length.out=(n+2))
    bb <- Brownian_bridge_path_sampler(x,y,0,t,times=times[2:n+1])
    m_t <- 0
    i <- 0
    resolved <- FALSE
    while(!resolved){
      i <- i+1
      m_t <- m_t + phiint(bb$full_path[1,i],bb$full_path[1,i+1],bb$full_path[2,i],bb$full_path[2,i+1])
      if (m_t>test_time){
        resolved <- TRUE
      } else if (i==n) {
        return(y)
      }
    }
  }
}

po_bound_alt <- function(x,y,t) {
  while(TRUE){
    y <-sim_end(t)
    p = 0
    n=1
    while((p/n)*(1-p/n)<=0){
      k <- rpois(1,t)
      skel_t <- runif(k,0,t)
      bb <- Brownian_bridge_path_sampler(
        x=x,
        y=y,
        s=0,
        t=t,
        times=skel_t
      )
      p = p+prod((-phi(bb$simulated_path[1,])+l))*exp(t)
      n=n+1
    }
    p <- p/n
    if(runif(1)<p){
      return(y)
    }
  }
}


ss=1e5
target_time = 2

start_time <- proc.time()
cts_smc_sample_r <- cts_smc_rust(ss,target_time)
print(proc.time()-start_time)
cts_smc_sample_r <- data.frame(w=cts_smc_sample_r)

start_time <- proc.time()
spline_output <- mcreplicate::mc_replicate(1e4,po_bound_alt(0,1,target_time))
print(proc.time()-start_time)
spline_output <- data.frame(w=spline_output)

plot <-
  ggplot()+
  geom_density(data=cts_smc_sample_r,aes(x=w),col='blue')+
  geom_density(data=spline_output,aes(x=w),col='red')+
  xlim(-12,12)
print(plot)