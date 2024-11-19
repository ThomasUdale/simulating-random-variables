# set.seed(1)
library(MASS)
library(ggplot2)
library(progress)
library(parallel)
library(layeredBB)
library(RcppZiggurat)
# devtools::build("~/Documents/simulating-random-variables/R-packages/cts.smc")
# devtools::load_all("~/Documents/simulating-random-variables/R-packages/cts.smc")

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

est_int<- function(n,path) {
  return(prod((phi(path)-l)/r))
}

sim_path <- function(path,new_points) {
  
}

series_sample <- function(t) {
  y <- sim_end(t)
  path <- data.frame(t=c(0,t),w=c(0,y))
  n <- 0
  f_u <- 1
  new_points <- runif(n+1,0,t)
  new_path <- sim_path(path,new_points)
  path <- update_path(path,new_points,new_path)
  f_d <- 1 - est_int(new_path)
  u <- runif(1)
  while ((u-f_d)*(f_u-u)>0) {
    n=n+2
    new_points <- runif(n,0,t)
    new_path <- sim_path(path,new_points)
    path <- update_path(path,new_points,new_path)
    f_u <- f_d + est_int(new_path)
    new_points <- runif(n+1,0,t)
    new_path <- sim_path(path,new_points)
    path <- update_path(path,new_points,new_path)
    f_d <- f_u + est_int(new_path)
  }
  if (u<f_d){
    return(path)
  }
}

series_sample(pi)

# ss=1e5
# target_time = pi
# 
# start_time <- proc.time()
# cts_smc_sample_r <- cts_smc_rust(ss,target_time)
# print(proc.time()-start_time)
# cts_smc_sample_r <- data.frame(w=cts_smc_sample_r)
# 
# plot <-
#   ggplot()+
#   # geom_density(data=exact,aes(x=w),col='green')+
#   # geom_density(data=euler_r,aes(x=w),col='red')+
#   geom_density(data=cts_smc_sample_r,aes(x=w),col='blue')+
#   xlim(-12,12)
# print(plot)
