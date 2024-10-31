library(layeredBB)


C = 2.0
c = 1
temp = 1.0

mu = c(0.5,-0.5)

x = 0
y = 0

fc <- function(x,c) {
  dnorm(x,mu[c],1)
}

phi <- function(x,c) {
  0.5 * ((x-mu[c])^2-1) + 0.5
}

find_bound <- function(bes_layer,c) {
  if ((mu[c]-bes_layer$L)<(bes_layer$U-mu[c])) {
    phi(bes_layer$U,c)
  } else {
    phi(bes_layer$L,c)
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

N = 500

ea3_out = mean(replicate(N,ea3(x,y,temp,phi,find_bound,c)))

po_est_bound <- function(x,y,temp,c) {
  bes_layer <- bessel_layer_simulation(x = x, y = y, s = 0, t = temp, mult = 0.2)
  rbound <- find_bound(bes_layer,c)
  k <- rpois(1,rbound*temp)
  
  skeleton = runif(k,0,temp)
  
  layered_bb <- layered_brownian_bridge(x = x,
                                        y = y,
                                        s = 0,
                                        t = temp,
                                        bessel_layer = bes_layer,
                                        times = skeleton)
  
  prod((rbound - phi(layered_bb$simulated_path[1,],c))/rbound)
}

po_est_bound_mean <- function(x,c,y,N,temp) {
  mean(replicate(N,po_est_bound(x,y,temp,c)))
}

po_bound_out = pbapply::pbreplicate(N,po_est_bound_mean(x,c,y,1,temp))


geo_sample <- function(x,c,y,temp) {
  bes_layer <- bessel_layer_simulation(x = x, y = y, s = 0, t = temp, mult = 0.2)
  rbound <- 10*find_bound(bes_layer,c)

  out <- 1
  for (i in 1:ceiling(rbound)*temp) {
    geo_x <- 0
    geo_sim <- TRUE
    while(geo_sim){
      geo_x <- geo_x +1
      layered_bb <- layered_brownian_bridge(x = x,
                                            y = y,
                                            s = 0,
                                            t = temp,
                                            bessel_layer = bes_layer,
                                            times = runif(1,0,temp))
      p_hat <- phi(layered_bb$simulated_path[1,],c)/ceiling(rbound)
      
      if (runif(1)<p_hat) {
        geo_sim <- FALSE
      }
    }
    if (sum(rexp(geo_x,1))<1) {
      out <- 0
    }
  }
  out
}

geo_sample_mean <- mean(pbapply::pbreplicate(N,geo_sample(x,c,y,temp)))

plot.new()
axis(1,xlim=c(-0.1,1.1))
axis(2,ylim=c(0,1.5))

hist(po_bound_out,freq=F)

abline(v=ea3_out,col='red')
text(x=ea3_out,y=3,ea3_out,col='red')
abline(v=mean(po_bound_out),col='green')
text(x=mean(po_bound_out),y=2,round(mean(po_bound_out),3),col='green')
abline(v=geo_sample_mean,col='blue')
text(x=geo_sample_mean,y=4,geo_sample_mean,col='blue')
