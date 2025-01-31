# set.seed(1)

library(layeredBB)
library(mcreplicate)
library(parallel)
library(ggplot2)

C = 4

f <- function(x) {
  exp(-x^4/2)/2.1558
}



fc <- function(x) {
  exp(-x^4/(2*C))
}

samplefc <- function() {
  while (TRUE) {
    x <- rnorm(1)
    if (runif(1)<(fc(x)/exp(-x^2/2+C/8))) {
      return(x)
    }
  }
}

phi <- function(x) {
  0.5*(4*x^6/C^2 - 6*x^2/C) + (2/C)^(1/2)
}

find_bound <- function(bes_layer) {
  m <- max(c(abs(bes_layer$L),bes_layer$U))
  if (m < (3*C/2)^(1/4)) {
    phi(0)
  } else {
    phi(m)
  }
}

ea3 <- function(x,y,t,phi,r) {
  bes_layer <- bessel_layer_simulation(x = x, y = y, s = 0, t = t, mult = 0.2)
  rbound <- r(bes_layer)
  k <- rpois(1,t*rbound)
  skeleton_u <- runif(k,0,rbound)
  skeleton_v <- runif(k,0,t)
  layered_bb <- layered_brownian_bridge(x = x,
                                        y = y,
                                        s = 0,
                                        t = t,
                                        bessel_layer = bes_layer,
                                        times = skeleton_v)
  all(phi(layered_bb$simulated_path[1,]) < skeleton_u)
}

po_est_bound <- function(x,y,temp) {
  bes_layer <- bessel_layer_simulation(x = x, y = y, s = 0, t = temp, mult = 0.2)
  rbound <- find_bound(bes_layer)
  k <- rpois(1,rbound*temp)
  # print(c(rbound,k))
  skeleton = runif(k,0,temp)
  
  layered_bb <- layered_brownian_bridge(x = x,
                                        y = y,
                                        s = 0,
                                        t = temp,
                                        bessel_layer = bes_layer,
                                        times = skeleton)
  
  prod((rbound - phi(layered_bb$simulated_path[1,]))/rbound)
}

mcf <- function(x,phi,r,t,c) {
  xbar <- mean(x[c(1:c)])
  xvar <- var(x[c(1:c)]) * (c-1)/c
  while (TRUE) {
    y = rnorm(1,xbar,sqrt(t/c))
    if (log(runif(1)) < -c*xvar/(2*t)) {
      if (all(mapply(ea3,x[c(1:c)],MoreArgs=c(y=y,t=1,r=r,phi=phi)))) {
        return(y)
      }
    }
  }
}


smc <- function(n_particles,temp){
  x <- array(replicate(n_particles*C,samplefc()),dim=c(n_particles,C))
  y <- rnorm(n_particles)*sqrt(temp) + x[,1]
  w_unnorm <- pbmcapply::pbmcmapply(po_est_bound,x=x[,1],y=y,MoreArgs=list(temp=temp),mc.cores=12)
  w_norm <- w_unnorm/sum(w_unnorm)
  z <- w_unnorm
  for (c in 2:C){
    print(c)
    if (sum(1/w_norm^2)<C/2) {
      A_t <- sample(c(1:n_particles),n_particles,replace=TRUE,prob=w_norm)
      w_hat_norm <- 1
    } else {
      A_t <- c(1:n_particles)
      w_hat_norm <- w_norm
    }
    x_new <- array(replicate(n_particles*C,samplefc()),dim=c(n_particles,C))[A_t,]
    y_new <- rnorm(n_particles)*sqrt(temp/c) + rowMeans(x_new[,c(1:c)])
    y_new <- array(replicate(c,y_new),dim=c(n_particles,c))
    z_new <- exp(-0.5*rowSums((x_new[,c(1:c)]-rowMeans(x_new[,c(1:c)]))^2)/temp)*apply(array(mapply(po_est_bound,x=x_new[,c(1:c)],y=y_new,MoreArgs=list(temp=temp)),dim=c(n_particles,c)),1,prod)
    for (i in 1:n_particles){
      if (runif(1)<z_new[i]/z[i]) {
        y[i] <- y_new[i]
        z[i] <- z_new[i]
        x[i,] <- x_new[i,]
      }
    }
    w_unnorm <- mapply(po_est_bound,x=x[,c],y=y,MoreArgs=list(temp=temp)) * w_hat_norm
    w_norm <- w_unnorm/sum(w_unnorm)
  }
  return(list(x=x,w_norm=w_norm,w_unnorm=w_unnorm,y=y))
}

output <- data.frame(smc(20000,1))

plot <- ggplot()+
  geom_histogram(data=output,aes(x=y,weight=w_norm,y=..density..),binwidth=0.1)+
  geom_function(fun = f, colour = "red")+
  xlim(-3, 3)+
  ylim(0,0.7)
print(plot)