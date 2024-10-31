# set.seed(1)
library(MASS)
library(ggplot2)
library(progress)
library(parallel)
library(layeredBB)
library(RcppZiggurat)

ss = 100000

interval_time = 8
target_time = 4

l = -0.5
r = 9/8
phi <- function(x) {
  (sin(x-pi)^2 + cos(x-pi))/2
}

sim_end <- function(t) {
  while (T) {
    u <- rnorm(1,0,sqrt(t))
    if (runif(1) < exp(-cos(u-pi)-1)) {
      return(u)
    }
  }
}

ea1 <- function(id,t,return_path=F){
  while (T){
    y <- sim_end(t)
    k <- rpois(1,r*t)
    if (k==0){
      if (return_path){
        return(data.frame(t=c(0,t),w=c(0,y)))
      } else {
        return(y)
      }
    }
    skel_t <- sort(runif(k,0,t))
    skel_u <- runif(k,0,r)
    bb <- Brownian_bridge_path_sampler(
      x=0,
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
print(Sys.time())
exact <- data.frame(w=unlist(parallel::mclapply(c(1:ss),ea1,t=target_time,mc.cores=8)))

ea1_inter <- function(id,t_target,t_interval,return_path=F){
  skeleton <- ea1(NaN,t_interval,T)
  end_ind <- which(skeleton$t>=t_target,arr.ind=T)[1]
  bb <- Brownian_bridge_path_sampler(
    x=skeleton$w[end_ind-1],
    y=skeleton$w[end_ind],
    s=skeleton$t[end_ind-1],
    t=skeleton$t[end_ind],
    times=t_target
  )
  
  if (!return_path){
    return(bb$simulated_path[1])
  } else {
    return(rbind(skeleton,c(t_target,b <- bb$simulated_path[1])))
  }
}

print(Sys.time())
# exact_i <- data.frame(w=unlist(pbmcapply::pbmclapply(c(1:ss),ea1_inter,t_target=target_time,t_interval=interval_time,mc.cores = 8)))

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
          x[killed_particle] <- dqrng::dqrnorm(1,x[killed_particle],sqrt(new_t-last_sample_time[killed_particle])) 
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

print(Sys.time())
cts_smc_sample=cts_smc(1,target_time,ss)
print(Sys.time())
cts_smc_sample <- data.frame(w=cts_smc_sample)
# cts_smc_sample <- data.frame(w=unlist(pbmcapply::pbmclapply(c(1:8),cts_smc,t_target=target_time,n_particles=ss/8,mc.cores = 8)))

print(Sys.time())
plot <-
  ggplot()+
  # geom_histogram(data=exact,aes(x=w,y=..density..),binwidth=0.1)+
  geom_density(data=exact,aes(x=w),col='green')+
  geom_density(data=exact_i,aes(x=w),col='blue')+
  geom_density(data=cts_smc_sample,aes(x=w),col='red')
print(plot)