### https://github.com/rchan26/layeredBB/tree/master?tab=readme-ov-file

library(layeredBB)
library(mcreplicate)
library(parallel)

par(mfrow = c(2, 2)) 

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
  0.5*(4*x^6/C^2 - 6*x^2/C)
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

mcf <- function(phi,r,t) {
  while (TRUE) {
    x <- replicate(C,samplefc())
    xbar <- mean(x)
    xvar <- var(x) * (C-1)/C
    y = rnorm(1,xbar,sqrt(t/C))
    if (log(runif(1)) < -C*xvar/(2*t)) {
      if (all(mapply(ea3,x,MoreArgs=c(y=y,t=1,r=r,phi=phi)))) {
        return(y)
      }
    }
  }
}

x<- pbapply::pbreplicate(20000,mcf(phi,find_bound,1))
curve(f,xlim = c(-3,3),ylim=c(0,0.8))
hist(x,add=T,freq=F,nclass = 50)

# po_est_bound <- function(x,y,temp) {
#   bes_layer <- bessel_layer_simulation(x = x, y = y, s = 0, t = temp, mult = 0.2)
#   rbound <- find_bound(bes_layer)
#   k <- rpois(1,temp*rbound)
#   
#   skeleton = runif(k,0,temp)
#   
#   layered_bb <- layered_brownian_bridge(x = x,
#                                         y = y,
#                                         s = 0,
#                                         t = temp,
#                                         bessel_layer = bes_layer,
#                                         times = skeleton)
#   
#   prod((rbound - phi(layered_bb$simulated_path[1,]))/rbound)
# }
# 
# po_est_bound_mean <- function(x,y,N,temp) {
#   mean(replicate(N,po_est_bound(x,y,temp)))
# }
# 
# pomf <- function(ss,phi,find_bound,temp,N,y0,g0,y_sd) {
#   sample <- c(y0)
#   y_current <- y0
#   g_current <- g0
#   jumps = 0
#   pb <- progress_bar$new(
#     format = "sampling [:bar] :current/:total (:percent) eta: :eta, elapsed: :elapsed, jumps :jumps",
#     clear = FALSE, total = ss, width = 80)
#   while (length(sample)<ss) {
#     pb$tick(tokens = list(jumps = jumps))
#     x <- replicate(C,samplefc())
#     y_new <- rnorm(1,mean(x),y_sd)
#     pc <- sapply(x,po_est_bound_mean,y=y_new,N=N,temp=temp)
#     g_new = exp(-0.5*(C-1)*var(x)/temp) * prod(pc)
#     
#     if (runif(1) < g_new/g_current) {
#       jumps = jumps + 1
#       y_current <- y_new
#       g_current <- g_new
#     }
#     sample <- append(sample,y_current)
#   }
#   pb$terminate()
#   return(sample)
# }
# 
# z3 <- pomf(25000,phi,find_bound,1,5,0,0.0001,sqrt(1/C))
# curve(f,xlim = c(-3,3),ylim=c(0,0.8))
# hist(z3[5000:25000],freq=F,add=T,col="red",nclass=50)
# 
