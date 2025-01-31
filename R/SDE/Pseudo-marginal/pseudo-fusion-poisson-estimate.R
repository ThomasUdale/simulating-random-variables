library(mcreplicate)
library(parallel)
library("ggplot2")
library(progress)
library(layeredBB)

C = 5.0
temp = 1.0

lambda = 1
gam =  1

x0 <- seq(-5, 5, by = 0.01)

f <- function(x) {
  exp(-(x)^4/2.0)/2.1558
}

fc <- function(x) {
  exp(-(x)^4/(2.0*C))/(2.1558)^(1/C)
}

samplefc <- function() {
  while (TRUE) {
    x <- rnorm(1)
    if (runif(1)<(fc(x)/exp(-x^2/2+C/8))) {
      return(x)
    }
  }
}

prop_y <- function(y) {
  rnorm(1,mean=y,sd=0.5)
}

prop_x <- function() {
  replicate(C,samplefc())
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

find_bound_approx <- function(bes_layer) {
  opt_func <- function(x) {
    -phi(x)
  }
  
  -optim(0.5*(bes_layer$L+bes_layer$U),opt_func,lower=bes_layer$L,upper=bes_layer$U,method='L-BFGS-B')$value
}

po_est_bound <- function(x,y) {
  bes_layer <- bessel_layer_simulation(x = x, y = y, s = 0, t = temp, mult = 0.2)
  rbound <- find_bound_approx(bes_layer)
  k <- rpois(1,temp*rbound)
  
  skeleton = runif(k,0,temp)
  
  layered_bb <- layered_brownian_bridge(x = x,
                                        y = y,
                                        s = 0,
                                        t = temp,
                                        bessel_layer = bes_layer,
                                        times = skeleton)

  prod((rbound - phi(layered_bb$simulated_path[1,]))/rbound)
}

po_est_bound_mean <- function(x,y,N) {
  mean(replicate(N,po_est_bound(x,y)))
}


x = replicate(C,samplefc())
y = 0
g_tilde = exp(-sum((y-x)^2/(2*T))) * prod(fc(x)) * exp(C*(lambda-gam)*temp)

sample = c()
jumps = 0

n_sample = 5000
n_burn = 2500



pb <- progress_bar$new(
  format = "sampling [:bar] :current/:total (:percent) eta: :eta, elapsed: :elapsed, jumps :jumps",
  clear = FALSE, total = n_sample+n_burn, width = 80)
for (i in 1:(n_sample+n_burn)) {
  pb$tick(tokens = list(jumps = jumps))
  y_new = prop_y(y)
  x_new = replicate(C,samplefc())
  
  g_tilde_new = exp(-sum((y_new-x_new)^2/(2*temp)))/dnorm(y_new,mean=y,sd=0.5)
  
  g2 = g_tilde_new
  
  pc_tildes = c()
  for (c in 1:C) {
    pc_tilde = po_est_bound_mean(x_new[c],y_new,10)
    
    pc_tildes = append(pc_tildes,pc_tilde)
  }
  
  g_tilde_new = g_tilde_new * prod(pc_tildes)
  
  alpha = min(1,g_tilde_new/g_tilde)
  
  
  if (runif(1)<alpha) {
    jumps = jumps + 1
    g_tilde = g_tilde_new
    x = x_new
    y = y_new
  }
  
  
  if (i>n_burn) {
    sample = append(sample,y)
    
  }
}

pb$terminate()

curve(f,-3,3)
hist(sample,freq=F,add=T,breaks=40)