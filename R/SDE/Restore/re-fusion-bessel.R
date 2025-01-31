set.seed(1)

library(glue)
library(ggplot2)
library(layeredBB)

C = 2


f <- function(x) {
  exp(-x^4/2)/2.1558
}

fc <- function(x) {
  exp(-x^4/(2*C))
}

samplefc <- function(c) {
  while (TRUE) {
    x <- rnorm(1)
    if (runif(1)<(fc(x)/exp(-x^2/2+C/8))) {
      return(x)
    }
  }
}

k_lower <- 0.1


k <- function(x) {
  0.5*(-6*x^2+4*x^6) + sqrt(2)*fc(x)/(fc(x)^C)
}

k_prime <- function(x) {
  k(x)-k_lower
}

upper_bound_k <- function(bes_layer) {
  if (max(abs(bes_layer$U),abs(bes_layer$L))<1.15){
    k_prime(1.15)
  } else {
    k_prime(max(abs(bes_layer$U),abs(bes_layer$L)))
  }
}

lower_bound_k <- function(bes_layer) {
  if ((0.81<(bes_layer$L)) | (-0.81>(bes_layer$U))){
    k_prime(min(abs(bes_layer$U),abs(bes_layer$L)))
  } else {
    0
  }
}

sim_pi <- function(z0,time_step) {
  sample_time <- rexp(1,k_lower)
  t <- 0
  while (TRUE) {
    x <- samplefc(sample(c(1:C),1))
    regen <- FALSE
    while(!regen) {
      y <- rnorm(1,x,sqrt(min(time_step,sample_time-t)))
      layer_t <- min(t+time_step,sample_time)
      bes_layer <- bessel_layer_simulation(x = x, y = y, s = t, t = layer_t, mult = 0.2)
      U_x <- upper_bound_k(bes_layer)
      L_x <- lower_bound_k(bes_layer)
      cur_time <- t
      k_times <- c()
      while (cur_time<layer_t){
        new_kill_time <- rexp(1,U_x-L_x)
        if (cur_time+new_kill_time<layer_t){
          cur_time <- cur_time+new_kill_time
        } else {
          cur_time <- layer_t
        }
        k_times <- append(k_times,cur_time)
      }
      
      new_x <- layered_brownian_bridge(x = x,
                                   y = y,
                                   s = t,
                                   t = layer_t,
                                   bessel_layer = bes_layer,
                                   times = k_times)$simulated_path[1,]
      in_layer <- TRUE
      stop_time_id <- 1
      while (in_layer) {
        if (k_times[stop_time_id]==sample_time){
          return(new_x[stop_time_id])
        } else if (k_times[stop_time_id]==layer_t){
          in_layer <- FALSE
          t <- k_times[stop_time_id]
          x <- new_x[stop_time_id]
        } else {
          t <- k_times[stop_time_id]
          p <- (U_x-k_prime(new_x[stop_time_id]))/(U_x - 0)
          if (runif(1)>p){
            regen <- TRUE
            in_layer <- FALSE
          }
        }
        stop_time_id = stop_time_id+1
      }
    }
  }
}

output <- sim_pi(1,0.05)
print(output)

samples <- pbmcapply::pbmclapply(seq(10000),sim_pi,time_step=0.05,mc.cores=10)
hist(unlist(samples),freq=F,nclass=50,ylim=c(0,0.6))
curve(f,xlim=c(-3,3),ylim=c(0,0.5),add=T)
lines(density(unlist(samples)),col='red')
curve(k_prime,add=T)
