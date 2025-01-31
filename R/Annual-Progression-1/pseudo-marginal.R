# devtools::load_all("~/Documents/simulating-random-variables/R-packages/bm.rust")
library(layeredBB)
library(progress)
library(ggplot2)



temp = 1
x = 0
y = 1

phi <- function(x) {
  0.5*(4*x^6 - 6*x^2) + sqrt(2)
}

find_bound <- function(bes_layer) {
  m <- max(c(abs(bes_layer$L),bes_layer$U))
  if (m < (3/2)^(1/4)) {
    phi(0)
  } else {
    phi(m)
  }
}

find_bound_k<-function(k) {
  if (k<0) {
    return(0)
  }
  if (k < (3/2)^(1/4)) {
    phi(0)
  } else {
    phi(k)
  }
}

po_est_bound <- function(x,y,temp) {
  bes_layer <- bessel_layer_simulation(x = x, y = y, s = 0, t = temp, mult = 0.2)
  rbound <- find_bound(bes_layer)
  k <- rpois(1,rbound*temp)
  
  skeleton = runif(k,0,temp)
  
  layered_bb <- layered_brownian_bridge(x = x,
                                        y = y,
                                        s = 0,
                                        t = temp,
                                        bessel_layer = bes_layer,
                                        times = skeleton)
  
  prod((rbound - phi(layered_bb$simulated_path[1,]))/rbound)
}

po_est_cut <- function(x,y,t) {
  N <- rgeom(1,1-exp(-1))
  gN <- find_bound_k(N)
  
  k <- rpois(1,gN*t)
  if (k==0){
    return(1)
  }
  skel_t <- runif(k,0,t)
  bb <- brownian_bridge(
    x=x,
    y=y,
    s=0,
    t=t,
    times=sort(skel_t)
  )
  
  phis <- phi(bb)
  g <- find_bound_k(0)
  z <- prod((gN - pmin(phis,g))/gN)
  if (N<=1){
    return(z)
  }
  po_old <- prod((gN - pmin(phis,g))/gN)
  for (i in 2:N) {
    g <- find_bound_k(i)
    po_new <- prod((gN - pmin(phis,g))/gN)
    if (po_new==po_old){
      break
    }
    z <- z+(po_new-po_old)*exp(i)
    po_old <- po_new
  }
  return(z)
}

bb_density <- function(path){
  len_path <- ncol(path)-1
  start_t <- path[2,1]
  end_t <- path[2,ncol(path)]
  mu <- path[1,1] + (path[2,2:len_path] - start_t)/(end_t-start_t)*(path[1,ncol(path)]-path[1,1])
  cov <- outer(
    path[2,2:len_path],
    path[2,2:len_path],
    Vectorize(function(x,y,start_t,end_t){(end_t-max(x,y))*(min(x,y)-start_t)/(end_t-start_t)}),start_t=start_t,end_t=end_t
  )
  mvtnorm::dmvnorm(path[1,2:len_path],mean=mu,sigma=cov)
}

pseudo_marginal <- function(x,y,t,N,burn){
  pb = progressBar(min = 0, max = N+burn, initial = 0) 
  bes_layer <- bessel_layer_simulation(x = x, y = y, s = 0, t = temp, mult = 0.2)
  rbound <- find_bound(bes_layer)
  k <- rpois(1,rbound*temp)
  skeleton = runif(k,0,temp)
  layered_bb <- layered_brownian_bridge(x = x,
                                        y = y,
                                        s = 0,
                                        t = temp,
                                        bessel_layer = bes_layer,
                                        times = skeleton)
  path <- layered_bb$full_path
  rho <- prod((rbound - phi(layered_bb$simulated_path[1,]))/rbound)
  q <- bb_density(path)
  output <- list()
  for (i in 1:(N+burn)){
    setTxtProgressBar(pb,i)
    bes_layer <- bessel_layer_simulation(x = x, y = y, s = 0, t = temp, mult = 0.2)
    rbound <- find_bound(bes_layer)
    k <- rpois(1,rbound*temp)
    skeleton = runif(k,0,temp)
    layered_bb <- layered_brownian_bridge(x = x,
                                          y = y,
                                          s = 0,
                                          t = temp,
                                          bessel_layer = bes_layer,
                                          times = skeleton)
    prop_path <- layered_bb$full_path
    prop_rho <- prod((rbound - phi(layered_bb$simulated_path[1,]))/rbound)
    prop_q <- bb_density(prop_path)
    accept <- min(1,(prop_rho*q)/(rho*prop_q))
    if (runif(1)<accept){
      path <- prop_path
      rho <- prop_rho
      q <- prop_q
    }
    
    if (i>burn){
      output[[i-burn]] <- path
    }
  }
  close(pb)
  return(output)
}

pseudo_marginal_neg <- function(x,y,t,N,burn){
  pb = progressBar(min = 0, max = N+burn, initial = 0) 
}



N = 1e4
# print('bounded poisson estimator')
# po_bound_out = pbapply::pbreplicate(N,po_est_bound(x,y,temp))
# print(summary(po_bound_out))
# print(sd(po_bound_out))
# print('fixed bound')
# po_bound_est_cut = pbapply::pbreplicate(N,po_est_cut(x,y,temp))
# print(summary(po_bound_est_cut))
# print(sd(po_bound_est_cut))

# output <- pseudo_marginal(x,y,temp,N,N)
# df <- NULL
# for (i in 1:N){
#   df <- rbind(df,data.frame(t(output[[i]]),col=i))
# }
# 
# plot <- ggplot(df)+
#   geom_line(aes(x=time,y=X,group=col),alpha=0.2)
# print(plot)



