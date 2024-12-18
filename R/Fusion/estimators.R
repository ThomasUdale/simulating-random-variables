library(layeredBB)

C = 2.0
c = 1
temp = 1

mu = c(0.5,-0.5)

x = 0
y = 1

# phi <- function(x,c) {
#   0.5 * ((x-mu[c])^2-1) + 0.5
# }
# 
# find_bound_k <- function(k) {
#   if (k<0) {return(0)}
#   if ((mu[c]+k)<(k-mu[c])) {
#     phi(k,c)
#   } else {
#     phi(-k,c)
#   }
# }
# 
# find_bound <- function(bes_layer,c) {
#   if ((mu[c]-bes_layer$L)<(bes_layer$U-mu[c])) {
#     phi(bes_layer$U,c)
#   } else {
#     phi(bes_layer$L,c)
#   }
# }

phi <- function(x,c) {
  0.5*(4*x^6 - 6*x^2) + sqrt(2)
}

find_bound <- function(bes_layer,c) {
  m <- max(c(abs(bes_layer$L),bes_layer$U))
  if (m < (3/2)^(1/4)) {
    phi(0,c)
  } else {
    phi(m,c)
  }
}

find_bound_k<-function(k) {
  if (k<0) {
    return(0)
  }
  if (k < (3/2)^(1/4)) {
    phi(0,c)
  } else {
    phi(k,c)
  }
}

# phi <- function(x,c) {
#   (sin(x)^2 + cos(x))/2 + 0.5
# }
# 
# find_bound <- function(bes_layer,c){
#   9/8
# }
# 
# find_bound_k <- function(k){
#   if(k<0){return(0)}
#   9/8
# }

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

po_est_cut <- function(x,y,t,c) {
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
  
  phis <- phi(bb,c)
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
  if(abs(z)>10){
    print(c(N,gN,k,z))
  }
  return(z)
}

compound_est <- function(x,y,t,c) {
  N <- rgeom(1,1-exp(-1))
  zc <- exp(N)/(1-exp(-1))
  gN <- find_bound_k(N)
  gN1 <- find_bound_k(N-1)
  k1 <- rpois(1,gN1*t)
  k <- rpois(1,(gN-gN1)*t)
  
  if (k+k1==0) {
    if (N==0){
      return(1*zc)
    }
    return(0*zc)
  }
  
  skel_t1 <- sort(runif(k1,0,t))
  skel_t <- append(skel_t1,runif(k,0,t))
  skel_0 <- order(skel_t)
  bb <- brownian_bridge(
    x=x,
    y=y,
    s=0,
    t=t,
    times=sort(skel_t)
  )
  
  phis <- phi(bb,c)
  
  y <- prod((gN-pmin(phis,gN))/gN)
  y1 <- prod((gN1-pmin(phis[skel_0<=k1],gN1))/gN1)
  
  if (N==0){
    y1 <- 0
  }
  
  z <- (y-y1)*zc
  
  return(z)
}

lim_est <- function(x,y,t,c){
  N <- rgeom(1,1-exp(-1))
  
  bb <- brownian_bridge(
    x=x,
    y=y,
    s=0,
    t=t,
    times=seq(0,t,length.out=(2^N+1))
  )
  
  phis = phi(bb,c)
  old_x <- 0
  z <- 0
  for (i in 0:N){
    if (i==N)
      I = t * mean(phis)
    else {
      I = t * mean(phis[seq_along(bb) %% 2^(N-i)==1])
    }
    
    xn <- (1-I/i)^i
    z <- z + (xn-old_x)*exp(i)
    old_x <- xn
  }
  return(z)
}


N = 1e4
print('bounded poisson estimator')
po_bound_out = pbapply::pbreplicate(N,po_est_bound_mean(x,c,y,1,temp))
print(summary(po_bound_out))
print(sd(po_bound_out))

print('fixed bound biased')
po_bound_est_cut = pbapply::pbreplicate(N,po_est_cut(x,y,temp,c))
print(summary(po_bound_est_cut))
print(sd(po_bound_est_cut))

print('single diff fixed bounds')
c_out = pbapply::pbreplicate(N,compound_est(x,y,temp,c))
print(summary(c_out))
print(sd(c_out))

print('exact:')
ea3_out = pbapply::pbreplicate(N,ea3(x,y,temp,phi,find_bound,c))
print(summary(as.integer(ea3_out)))
print(sd(ea3_out))

plot.new()
axis(1,xlim=c(-0.1,1.1))
axis(2,ylim=c(0,1.5))

hist(as.integer(ea3_out),freq=F,nclass=50,ylim=c(0,1),xlim=c(0,1))

abline(v=mean(ea3_out),col='red')
text(x=mean(ea3_out),y=0.4,mean(ea3_out),col='red')
abline(v=mean(po_bound_out),col='green')
text(x=mean(po_bound_out),y=0.5,round(mean(po_bound_out),3),col='green')
abline(v=mean(po_bound_est_cut),col='blue')
text(x=mean(po_bound_est_cut),y=0.3,round(mean(po_bound_est_cut),3),col='blue')
abline(v=mean(c_out),col='black')
text(x=mean(c_out),y=0.2,round(mean(c_out),3),col='black')