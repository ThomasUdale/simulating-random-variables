set.seed(10)

library(glue)
library(ggplot2)
library(layeredBB)
library(Rmpfr)
library(mcreplicate)
library(pbmcapply)

C = 4


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

k_lower <- 1

k <- function(x) {
  0.5*(-6*x^2+4*x^6) + (k_lower +sqrt(2))*fc(x)^(1-C)
}


k_prime <- function(x) {
  k(x)-k_lower
}

upper_bound_k <- function(x,x_extreme) {
  if (max(abs(x),abs(x_extreme))<1){
    k(1)
  } else {
    k(max(abs(x),abs(x_extreme)))
  }
}

lower_bound_k <- function(x,x_extreme) {
  large_x <- max(x,x_extreme)
  small_x <- min(x,x_extreme)
  if ((0.8<(small_x)) | (-0.8>(large_x))){
    k(min(abs(small_x),abs(large_x)))
  } else {
    0
  }
}

samples <- c()

M1 = 0.42599
M2 = 0.578103
M = M1+M2
t_star = 0.64


g_bar <- function(tau) {
  if (tau<t_star){
    2/(pi*tau^(3/2))*exp(-1/(2*tau))/M
  } else {
    pi/2*exp(-(pi^2*tau)/8)/M
  }
}

f_up <- function(tau,n) {
  min(pi*sum((-1)^seq(from=0,to=2*n)*(seq(from=0,to=2*n)+1/2)*exp(-0.5*(seq(from=0,to=2*n)+0.5)^2*pi^2*tau)),g_bar(tau))
}

f_down <- function(tau,n) {
  max(pi*sum((-1)^seq(from=0,to=2*n+1)*(2/(pi*tau))^3/2*(seq(from=0,to=2*n+1)+0.5)*exp((-2/tau)*(seq(from=0,to=2*n+1)+0.5)^2)),0)
}



algo4 <- function(theta,w0){
  while(TRUE) {
    if (runif(1)<M1/M){
      X <- rexp(1)
      tau_bar <- t_star + 8*X/pi^2
    } else {
      found_x <- FALSE
      X <- rexp(1)
      i <- 1
      while (!found_x){
        next_x <- rexp(1)
        if (X^2<=2*next_x/t_star & (i %% 2)==1) {
          found_x <- TRUE
        } else {
          X <- next_x
          i <- i+1
        }
      }
      tau_bar <- t_star/(1+t_star*X)^2
    }
    u <- runif(1)
    n <- 0
    gbar <- g_bar(tau_bar)
    while (gbar>f_down(tau_bar,n) & gbar<f_up(tau_bar,n)){
      n <- n + 1
    }
    if (gbar <= f_down(tau_bar,n)) {
      tau <- theta^2*tau_bar
      if(runif(1)<0.5){
        return(list(tau=tau,w=w0+theta))
      } else {
        return(list(tau=tau,w=w0-theta))
      }
    }
  }
}

psi1<- function(t_change,j,w_change,theta){
  exp(-(2*theta^2*(2*j-1)^2)/(t_change)+(2*(2*j-1)*theta*w_change)/t_change) + exp(-(2*theta^2*(2*j-1)^2)/(t_change)-(2*(2*j-1)*theta*w_change)/t_change)
}

psi2<- function(t_change,j,w_change,theta){
  exp(-(8*theta^2*j^2)/(t_change)+(4*j*theta*w_change)/t_change) + exp(-(8*theta^2*j^2)/(t_change)-(4*j*theta*w_change)/t_change)
}

psi3<- function(t_change,j,w_change,theta,m){
  (4*theta*j + m*w_change)/(m*w_change)*exp(-(4*theta*j*(2*theta*j + m*w_change))/t_change)
}

pn_down <- function(n,w_start,w_end,w_sim,theta,t_end,t_sim,t_start,m){
  n0 <- ceiling(sqrt((t_end-t_sim)+4*theta^2)/(4*theta))
  out <- if(n==1){
    (1-psi1(t_sim-t_start,1,w_start-w_sim,theta))/(1 - exp(-(2*theta*(m*(w_start-w_sim)+theta))/(t_sim-t_start)))*(1+sum(psi3(t_end-t_sim,seq(n0+1),w_sim-w_end,theta,m))+sum(psi3(t_end-t_sim,seq(n0),w_sim-w_end,theta,-m)))
  } else {
    (1-sum(psi1(t_sim-t_start,seq(n),w_start-w_sim,theta))+sum(psi2(t_sim-t_start,seq(n-1),w_start-w_sim,theta)))/(1 - exp(-(2*theta*(m*(w_start-w_sim)+theta))/(t_sim-t_start)))*(1+sum(psi3(t_end-t_sim,seq(n0+n),w_sim-w_end,theta,m))+sum(psi3(t_end-t_sim,seq(n0+n-1),w_sim-w_end,theta,-m)))
  }
  if (is.na(out)){
    print('')
    print('---')
    print(c(out,n,n0))
    
  }
  out
}

pn_down2 <- function(n,w_start,w_end,w_sim,theta,t_end,t_sim,t_start,m) {
  n0 <- ceiling(sqrt((t_end-t_sim)+4*theta^2)/(4*theta))
  out <- if (n==1){
    (
      1/(1 - exp(-(2*theta*(m*(w_start-w_sim)+theta))/(t_sim-t_start)))
      -(exp(-(2*theta^2*(2*1-1)^2)/(t_sim-t_start)+(2*(2*1-1)*theta*(w_start-w_sim))/(t_sim-t_start)+(2*theta*(m*(w_start-w_sim)+theta))/(t_sim-t_start)) + exp(-(2*theta^2*(2*1-1)^2)/((t_sim-t_start))-(2*(2*1-1)*theta*(w_start-w_sim))/(t_sim-t_start)+(2*theta*(m*(w_start-w_sim)+theta))/(t_sim-t_start)))/(exp((2*theta*(m*(w_start-w_sim)+theta))/(t_sim-t_start))-1)
    )*(
       1+sum(psi3(t_end-t_sim,seq(n0+1),w_sim-w_end,theta,m))+sum(psi3(t_end-t_sim,seq(n0),w_sim-w_end,theta,-m))
    )
  } else {
    (
      1/(1 - exp(-(2*theta*(m*(w_start-w_sim)+theta))/(t_sim-t_start)))
      -sum((exp(-(2*theta^2*(2*seq(n)-1)^2)/(t_sim-t_start)+(2*(2*seq(n)-1)*theta*(w_start-w_sim))/(t_sim-t_start)+(2*theta*(m*(w_start-w_sim)+theta))/(t_sim-t_start)) + exp(-(2*theta^2*(2*seq(n)-1)^2)/((t_sim-t_start))-(2*(2*seq(n)-1)*theta*(w_start-w_sim))/(t_sim-t_start)+(2*theta*(m*(w_start-w_sim)+theta))/(t_sim-t_start)))/(exp((2*theta*(m*(w_start-w_sim)+theta))/(t_sim-t_start))-1))
      + sum(exp(-(8*theta^2*seq(n-1)^2)/(t_sim-t_start)+(4*seq(n-1)*theta*(w_start-w_sim))/(t_sim-t_start)+(2*theta*(m*(w_start-w_sim)+theta))/(t_sim-t_start)) + exp(-(8*theta^2*seq(n-1)^2)/(t_sim-t_start)-(4*seq(n-1)*theta*(w_start-w_sim))/(t_sim-t_start)+(2*theta*(m*(w_start-w_sim)+theta))/(t_sim-t_start)))/(exp((2*theta*(m*(w_start-w_sim)+theta))/(t_sim-t_start))-1)
    )*(
      1+sum(psi3(t_end-t_sim,seq(n0+n),w_sim-w_end,theta,m))+sum(psi3(t_end-t_sim,seq(n0+n-1),w_sim-w_end,theta,-m))
    )
  }
  if (is.na(out)){
    return(1)
  } else {
    out
  }
}

pn_up <- function(n,w_start,w_end,w_sim,theta,t_end,t_sim,t_start,m){
  n0 <- ceiling(sqrt((t_end-t_sim)+4*theta^2)/(4*theta))
  out <- (1-sum(psi1(t_sim-t_start,seq(n),w_start-w_sim,theta))+sum(psi2(t_sim-t_start,seq(n),w_start-w_sim,theta)))/(1 - exp(-(2*theta*(m*(w_start-w_sim)+theta))/(t_sim-t_start)))*(1+sum(psi3(t_end-t_sim,seq(n0+n),w_sim-w_end,theta,m))+sum(psi3(t_end-t_sim,seq(n0+n),w_sim-w_end,theta,-m)))
  out
}

algo5 <- function(w_start,w_end,theta,t_start,t_end,t_sim){
  if (t_end==t_sim){
    return(w_end)
  }
  
  if (t_sim==t_start){
    return(w_start)
  }
  
  while (TRUE){
    b <- rnorm(3,0,sqrt(((t_end-t_sim)*(t_sim-t_start))/(t_end-t_start)^2))
    w_sim <- w_start + (-1)^as.integer(w_end<w_start)*sqrt((t_end-t_start)*((theta*(t_end-t_sim)/(t_end-t_start)^(3/2)+b[1])^2 + b[2]^2 +b[3]^2))
    u <- runif(1)
    n <- 1
    m <- as.integer(w_end>w_start)-as.integer(w_end<w_start)
    while (u>pn_down2(n,w_start,w_end,w_sim,theta,t_end,t_sim,t_start,m) & u<pn_up(n,w_start,w_end,w_sim,theta,t_end,t_sim,t_start,m)){
      n <- n + 1
    }
    if (u<pn_down2(n,w_start,w_end,w_sim,theta,t_end,t_sim,t_start,m)){
      return(w_sim)
    }
  }
  
}

sim_pi <- function(z,theta) {
  t <- 0
  i <- 0
  j <- -1
  while (TRUE) {
    x <- samplefc(sample(c(1:C),1))
    regen<-FALSE
    while (!regen) {
      j <- j + 1
      exit_info <- algo4(theta,x)
      U_x <- upper_bound_k(x,exit_info$w)
      # L_x <- lower_bound_k(x,exit_info$w)
      hit_out <- FALSE
      while(!hit_out) {
        i <- i+1
        E <- rexp(1,U_x-k_lower)
        sample_time <- rexp(1,k_lower)
        new_t <- min(t+E,t+exit_info$tau,t+sample_time)
        new_z <- algo5(x,exit_info$w,theta,t,min(t+exit_info$tau,t+sample_time),new_t)
        if (new_t==t+exit_info$tau){
          x <- new_z
          hit_out<-TRUE
        } else if (new_t==t+sample_time){
          x <- new_z
          return(x)
        } else {
          p <- (U_x - k(new_z))/(U_x-k_lower)
          if (runif(1)<p) {
            x <- new_z
          } else {
            hit_out<-TRUE
            regen<-TRUE
          }
        }
        t <- new_t
      }
    }
  }
}

# output <- sim_pi(1,0.1)

samples <- pbmcapply::pbmclapply(seq(10000),sim_pi,theta=0.1,mc.cores=12)

# samples <- pbapply::pbreplicate(1000,sim_pi(1,0.1))
# samples <- mcreplicate::mc_replicate(50000,sim_pi(0.1))

curve(f,xlim=c(-3,3))
hist(unlist(samples),freq=F,add=T,nclass=50)

