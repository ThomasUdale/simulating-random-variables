library(layeredBB)
library(ggplot2)
library(mvtnorm)
library(progress)

f <- function(x) {
  dnorm(x,4,sqrt(2/3))
}


mu = c(10,1)
sigma = c(sqrt(2),1)

C = length(mu)

fc <- function(x,c) {
  dnorm(x,mu[c],sigma[c])
}

samplefc <- function(c) {
  rnorm(1,mu[c],sigma[c])
}

phi <- function(x,c) {
  0.5 * ((x-mu[c])^2/sigma[c]^2-1) + 0.5
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
  -(x-mu[c])/sigma[c]^2
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
  -(x-mu[c])/sigma[c]^2
}
b_p <- function(x,c){
 -1/sigma[c]^2
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

ctsis_bf <- function(id,t){
  tau <- 0
  x <- sapply(1:C, samplefc)
  theta <- cbind(mapply(b,x=x,c=c(1:C)),gamma(x))
  w <- exp(-sum((x-mean(x))^2)/(2*t))
  while(TRUE) {
    event_time <- (-0.2*log(runif(1))/1)^(1/0.2)
    new_tau <- tau + event_time
    print(new_tau)
    if (new_tau>t){
      for (c in 1:C){
        w <- w*q_den(theta[c,],x[c],mean(x),t-tau)
      }
      return(c(mean(x),w))
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
  tau <- 0
  x <- matrix(sapply(rep(1:C,times=n),samplefc),nrow=C)
  theta_b <- cbind(sapply(c(1:C),function(y) b(x[y,],y)))
  w <- exp(-colSums(scale(x,scale=F)^2)/2*t)
  while(tau<t) {
    event_time <- (-0.3*log(runif(1))/1)^(1/0.3)
    new_tau <- tau + event_time
    if (new_tau>t){
      new_x <- colMeans(x)
      for (c in 1:C){
        for( i in 1:n){
          w[i] <- w[i]*q_den(c(theta_b[i,c],1),x[c,i],new_x[c],t-tau)
        }
      }
      return(rbind(new_x,w))
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
      # resampling
      # ess <- sum(abs(w))/(sum(w^2))
      # print(c(tau,ess))
      # new_ids <- sample(n,replace=T,prob=abs(w))
      # w <- sign(w[new_ids])*sum(abs(w))/n
      # x <- x[,new_ids]
      print(c(tau,sum(abs(w))))
      
      theta_b <- cbind(sapply(c(1:C),function(y) b(x[y,],y)))
      tau <- new_tau
      
    }
  }
}

ptm <- proc.time()
bf_out <- do.call(rbind,pbmcapply::pbmclapply(c(1:1e4),ctsis_bf,t=1,mc.cores=8))
print(proc.time()-ptm)
print(mean(bf_out[1,]*bf_out[2,]))
ggplot()+xlim(0,10)+geom_point(aes(x=bf_out[1,],y=bf_out[2,]))+geom_density(aes(bf_out[1,]))+geom_function(fun=f)
# ptm <- proc.time()
# mcf_out <- pbapply::pbreplicate(1e4,mcf(phi,find_bound,1))
# print(proc.time()-ptm)
# print(mean(mcf_out))
# 
# ptm <- proc.time()
# out <- t(ctsis_bf_resample(1,1e4))
# print(proc.time()-ptm)
# print(mean(out[,1]*out[,2]))


#re-weighting to eliminate negative weights. 
# 
# out2 <- out[order(out[,1]),]
# 
# while(min(out2[,2])<0) {
#   print(min(out2[,2]))
#   cur_id <- which.min(out2[,2])
#   resolved <- FALSE
#   idx <- 0
#   l <- cur_id
#   u <- cur_id
#   while (!resolved){
#     idx <- -idx
#     if (idx>=0){
#       idx <- idx+1
#       u <- min(1e5,cur_id+idx)
#     } else {
#       l <- max(1,cur_id+idx)
#     }
#     s <- sum(out2[l:u,2])
#     if (s>0) {
#       out2[l:u,2] <- s/(u-l)
#       resolved <- TRUE
#     }
#   }
# }
# 
# plot <- ggplot()+
#   geom_histogram(aes(x=mcf_out,y=after_stat(density)))+
#   geom_vline(xintercept=mean(bf_out[1,]*bf_out[2,]))+
#   geom_density(aes(x=out2[,1],weight=out2[,2]),col='red')+
#   geom_function(fun=f)
# print(plot)