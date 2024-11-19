# set.seed(1)
library(MASS)
library(ggplot2)
library(progress)
library(parallel)
library(layeredBB)
library(RcppZiggurat)
# devtools::build("~/Documents/simulating-random-variables/R-packages/cts.smc")
devtools::load_all("~/Documents/simulating-random-variables/R-packages/cts.smc")

cts_smc_rust(1,1)


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

ea1 <- function(id,x,t,return_path=F){
  while (T){
    y <- sim_end(t)
    k <- rpois(1,r*t)
    if (k==0){
      if (return_path){
        return(data.frame(t=c(0,t),w=c(x,y)))
      } else {
        return(y)
      }
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
      if (return_path){
        return(data.frame(t=bb$full_path[2,],w=bb$full_path[1,]))
      } else {
        return(y)
      }
    }
  }
}

cts_smc <- function(id,t_target,n_particles){
  t <- 0
  x <- array(0,dim=c(n_particles))
  last_sample_time <- array(0,dim=c(n_particles))
  while (t < t_target) {
    k_plus <- rexp(1,n_particles*(5/8))
    k_minus <- rexp(1,n_particles*(1/2))
    new_t = min(t_target,t+k_plus,t+k_minus)
    test_particle <- dqrng::dqsample.int(n_particles,1)
    x[test_particle] <- dqrng::dqrnorm(1,x[test_particle],sqrt(new_t-last_sample_time[test_particle])) 
    last_sample_time[test_particle] <- new_t
    if (k_plus<k_minus) {
      if (runif(1) < max(phi(x[test_particle]),0)/(5/8)) {
        dup_particle <- dqrng::dqsample.int(n_particles,1)
        if (dup_particle!=test_particle){
          x[dup_particle] <- dqrng::dqrnorm(1,x[dup_particle],sqrt(new_t-last_sample_time[dup_particle])) 
          last_sample_time[dup_particle] <- new_t
          x[test_particle] <- x[dup_particle]
        }
      }
    } else {
      if (runif(1) < max(-phi(x[test_particle]),0)/(1/2)) {
        killed_particle <- dqrng::dqsample.int(n_particles,1)
        if (killed_particle!=test_particle){
          last_sample_time[killed_particle] <- new_t
          x[killed_particle] <- x[test_particle]
        }
      }
    }
    t <- new_t
  }
  x <- rnorm(n_particles,x,sqrt(t-last_sample_time))
  x <- x[sample(n_particles,n_particles,replace=T,prob=exp(-cos(x-pi)+cos(0-pi)))]
  return(x)
}

euler_approx <- function(id,t_target,n_intervals) {
  eps <- t_target/n_intervals
  x <- 0
  for (i in 1:n_intervals) {
    x <- x + sin(x)*eps + rnorm(1,0,sqrt(eps))
  }
  return(x)
}

ss=1e5
target_time = 10*pi

# start_time <- proc.time()
# exact = pbmcapply::pbmclapply(c(1:ss),ea1,x=0,t=target_time,mc.cores=8)
# print(proc.time()-start_time)
# exact <- data.frame(w=unlist(exact))

# start_time <- proc.time()
# euler = pbmcapply::pbmclapply(c(1:ss),euler_approx,t_target=target_time,n_intervals=1e4,mc.cores=8)
# print(proc.time()-start_time)
# euler <- data.frame(w=unlist(euler))

start_time <- proc.time()
cts_smc_sample_r <- cts_smc_rust(ss,target_time)
print(proc.time()-start_time)
cts_smc_sample_r <- data.frame(w=cts_smc_sample_r)

start_time <- proc.time()
euler_r <- euler_approx_rust(ss,target_time,1e4)
print(proc.time()-start_time)
euler_r <- data.frame(w=euler_r)


plot <-
  ggplot()+
  # geom_density(data=exact,aes(x=w),col='green')+
  geom_density(data=euler_r,aes(x=w),col='red')+
  geom_density(data=cts_smc_sample_r,aes(x=w),col='blue')+
  xlim(-12,12)
print(plot)