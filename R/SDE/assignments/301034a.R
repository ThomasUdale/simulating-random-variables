library(MASS)
library(ggplot2)
library(progress)
library(parallel)
library(layeredBB)
devtools::load_all("~/Documents/simulating-random-variables/R-packages/cts.smc")
mu =1 

f <- function(x) {
  exp(-(x)^4/2.0)/2.1558
}

phi <- function(x) {
  0.5*(4*x^6 - 6*x^2)
}

l = -sqrt(2)

find_bound <- function(bes_layer) {
  m <- max(c(abs(bes_layer$L),bes_layer$U))
  if (m < (3/2)^(1/4)) {
    phi(0)
  } else {
    phi(m)
  }
}

sim_h <- function(t) {
  while(T) {
    x <- rnorm(1,0,sqrt(t))
    if (runif(1)<exp(-x^4/2)) {
      return(x)
    }
  }
}

ea3 <- function(x,t) {
  while(T){
    y <- sim_h(t)
    bes_layer <- bessel_layer_simulation(x = x, y = y, s = 0, t = t, mult = 0.2)
    rbound <- find_bound(bes_layer)
    k <- rpois(1,t*rbound)
    skeleton_u <- runif(k,0,rbound)
    skeleton_v <- runif(k,0,t)
    layered_bb <- layered_brownian_bridge(x = x,
                                          y = y,
                                          s = 0,
                                          t = t,
                                          bessel_layer = bes_layer,
                                          times = skeleton_v)
    if (all((sapply(layered_bb$simulated_path[1,],phi)) < skeleton_u+l)) {
      return(y)
    }
  }
}

cts_smc <- function(n_particles,target_time,theta) {
  t <- 0
  x <- array(0,dim=c(n_particles))
  last_sample_time <- array(0,dim=c(n_particles))
  while(t<target_time){
    print(t)
    layer_ends <- x + sample(c(-theta,theta),n_particles,replace=T)
    layer_info <- list(U=max(layer_ends),L=min(layer_ends))
    passage_times <- brownian_motion_fpt(n_particles)*theta^2
    min_pt_pid <- which.min(passage_times)
    min_pt <- passage_times[min_pt_pid]
    r <- find_bound(layer_info)
    k_plus = if (r>0){
      rexp(1,r*n_particles)
    } else {
      Inf
    }
    k_minus = if (-l>0){
      rexp(1,-l*n_particles)
    } else {
      Inf
    }
    k <- min(k_plus,k_minus)
    kid <- sample(n_particles,1)
    new_t = min(t+k,t+min_pt,target_time)
    for (id in 1:n_particles) {
      x[id] <- algo5(x[id],layer_ends[id],t,passage_times[id],new_t)
    }
    if (new_t==target_time){
      return(x)
    }
    
    if (new_t==t+k) {
      if (k_plus<k_minus) {
        if (runif(1)<max(phi(x[kid]),0)/r) {
          x[kid] <- x[sample(n_particles,1)]
        }
      } else {
        if (runif(1)<max(-phi(x[kid]),0)/(-l)) {
          x[sample(n_particles,1)] <- x[kid]
        }
      }
    }
    t <- new_t
  }
}

ss= 1000
target_time=1

# start_time <- proc.time()
# exact <- replicate(ss,ea3(0,target_time))
# print(proc.time()-start_time)


start_time <- proc.time()
output_l <- cts_smc(ss,target_time,0.5)
print(proc.time()-start_time)


# plot <-
#   ggplot()+
#   geom_density(data=data.frame(w=exact),aes(x=w),col='green')+
#   geom_density(data=data.frame(w=output_l),aes(x=w),col='red')+
#   xlim(-4,4)
# print(plot)
