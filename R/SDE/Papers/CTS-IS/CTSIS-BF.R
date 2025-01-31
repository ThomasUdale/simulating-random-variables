library(layeredBB)
library(ggplot2)
library(mvtnorm)
library(progress)

f <- function(x) {
  dnorm(x,0,sqrt(0.5))
}

mu_0 = 1
mu = c(0.5*mu_0,-0.5*mu_0)
C = length(mu)

fc <- function(x,c) {
  dnorm(x,mu[c],1)
}

samplefc <- function(c) {
  rnorm(1,mu[c],1)
}

phi <- function(x,c) {
  0.5 * ((x-mu[c])^2-1) + 0.5
}



find_bound <- function(bes_layer,c) {
  if ((mu[c]-bes_layer$l)<(bes_layer$u-mu[c])) {
    phi(bes_layer$u,c)
  } else {
    phi(bes_layer$l,c)
  }
}

ea3 <- function(x,y,t,phi,r,c) {
  bes_layer <- bessel_layer_simulation(x = x, y = y, s = 0, t = t, mult = 0.2)
  rbound <- r(bes_layer,c)
  k <- rpois(1,t*rbound)
  skeleton_u <- runif(k,0,rbound)
  skeleton_v <- runif(k,0,t)
  layered_bb <- layered_brownian_bridge(x = x,
                                        y = y,
                                        s = 0,
                                        t = t,
                                        bessel_layer = bes_layer,
                                        times = skeleton_v)
  all(sapply(layered_bb$simulated_path[1,],phi,c=c) < skeleton_u)
}

mcf <- function(phi,r,t) {
  while (TRUE) {
    x <- sapply(c(1:C),samplefc)
    xbar <- mean(x)
    xvar <- var(x) * (C-1)/C
    y = rnorm(1,xbar,sqrt(t/C))
    if (log(runif(1)) < -C*xvar/(2*t)) {
      if (all(mapply(ea3,x,c=c(1:C),MoreArgs=c(y=y,t=t,r=r,phi=phi)))) {
        return(y)
      }
    }
  }
}

alpha <- function(x,c){
  -(x-mu[c])
}
lambda <- function(s) {
  1 * s^(0.5-1)
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
b <- function(x,c) {
  -(x-mu[c])
}
b_p <- function(x,c){
 -1
}

rho <- function(x,y,u,c){
  1 + 1/lambda(u)*(1/(gamma(x)*u)*((gamma(y)-gamma(x))/2*(((y-x-u*b(x,c))^2)/(gamma(x)*u)-1)+(b(y,c)-b(x,c)-gamma_p(y))*(y-x-u*b(x,c))) + gamma_pp(y)/2 - b_p(y,c))
}

g_prop <- function(x,s,tt,t,C){
  rmvnorm(
    n=1,
    mean=(tt-t)/(tt-s)*x + (t-s)/(tt-s)*mean(x),
    sigma=diag((t-s)*(tt-t)/(tt-s),nrow=C)+matrix(rep((t-s)^2/(C*(tt-s)),C^2),nrow=C)
  )
}

g_den <- function(y,x,s,tt,t,C) {
  dmvnorm(
    x=y,
    mean=(tt-t)/(tt-s)*x + (t-s)/(tt-s)*mean(x),
    sigma=diag((t-s)*(tt-t)/(tt-s),nrow=C)+matrix(rep((t-s)^2/(C*(tt-s)),C^2),nrow=C)
  )
}

q_den <- function(theta,x,y,event_time) {
  dnorm(
    x=y,
    mean=x+event_time*theta[1],
    sd=sqrt(event_time*theta[2])
  )
}

ctsis_bf <- function(t){
  tau <- 0
  x <- sapply(1:C, samplefc)
  theta <- cbind(mapply(b,x=x,c=c(1:C)),gamma(x))
  w <- exp(-sum((x-mean(x))^2)/(2*t))
  while(TRUE) {
    event_time <- (-0.5*log(runif(1))/1)^(1/0.5)
    new_tau <- tau + event_time
    if (new_tau>t){
      x <- g_prop(x,tau,t,t,C)
      return(c(x[1],w))
    } else {
      new_x <- g_prop(x,tau,t,new_tau,C)
      w <- w/g_den(new_x,x,tau,t,new_tau,C)
      for (c in 1:C){
        w <- w*rho(x[c],new_x[c],event_time,c)*q_den(theta[c,],x[c],new_x[c],event_time)
      }
      x <- new_x
      theta <- cbind(mapply(b,x=x,c=c(1:C)),gamma(x))
      tau <- new_tau
    }
  }
}

ctsis_bf_resample <- function(t,n) {
  pb <- progress_bar$new(
    format = "sampling [:bar] (:percent) eta: :eta, elapsed: :elapsed",
    clear = FALSE, total = t, width = 80,force = TRUE)
  tau <- 0
  x <- matrix(sapply(rep(1:C,times=n),samplefc),nrow=C)
  theta_b <- cbind(sapply(c(1:C),function(y) b(x[y,],y)))
  w <- exp(-colSums(scale(x,scale=F)^2)/2*t)
  while(tau<t) {
    event_time <- (-0.5*log(runif(1))/1)^(1/0.5)
    new_tau <- tau + event_time
    if (new_tau>t){
      x <- apply(x,2,g_prop,s=tau,tt=t,t=t,C)
      pb$terminate()
      return(rbind(x[1,],w))
    } else {
      new_x <- apply(x,2,g_prop,s=tau,tt=t,t=new_tau,C)
      for(i in 1:n){
        w[i] <- w[i]/g_den(new_x[,i],x[,i],tau,t,new_tau,C)
      }
      for (c in 1:C){
        for( i in 1:n){
          w[i] <- w[i]*rho(x[c,i],new_x[c,i],event_time,c)*q_den(c(theta_b[i,c],1),x[c,i],new_x[c,i],event_time)
        }
      }
      x <- new_x
      theta_b <- cbind(sapply(c(1:C),function(y) b(x[y,],y)))
      tau <- new_tau
      pb$update(tau/t)
    }
  }
}

# ptm <- proc.time()
# bf_out <- pbapply::pbreplicate(1e4,ctsis_bf(1))
# print(proc.time()-ptm)
# print(mean(bf_out[1,]*bf_out[2,]))

# ptm <- proc.time()
# mcf_out <- pbapply::pbreplicate(1e4,mcf(phi,find_bound,1))
# print(proc.time()-ptm)
# print(mean(mcf_out))

# ptm <- proc.time()
# out <- ctsis_bf_resample(1,1e4)
# print(proc.time()-ptm)
# print(mean(out[1,]*out[2,]))
# print(ggplot()+geom_point(aes(x=out[1,],y=out[2,])))

# 
# plot <- ggplot()+
#   geom_histogram(aes(x=mcf_out,y=after_stat(density)))+
#   geom_vline(xintercept=mean(bf_out[1,]*bf_out[2,]))+
#   geom_function(fun=f)
# print(plot)


