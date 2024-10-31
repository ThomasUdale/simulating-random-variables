# set.seed(1)
library(progress)
library(layeredBB)
library(glue)

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

find_bound_approx <- function(bes_layer,c) {
  opt_func <- function(x) {
    -phi(x,c)
  }
  
  -optim(mu[c],opt_func,lower=bes_layer$L,upper=bes_layer$U,method='L-BFGS-B')$value
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

ea_est<- function(x,c,y,N,temp) {
  mean(replicate(N,ea3(x,y,temp,phi,find_bound,c)))
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

pomf <- function(ss,phi,find_bound,temp,N,y0,g0) {
  sample <- c()
  y_current <- y0
  g_current <- g0
  pb <- progress_bar$new(
    format = "sampling [:bar] :current/:total (:percent) eta: :eta, elapsed: :elapsed",
    clear = FALSE, total = ss, width = 80)
  while (length(sample)<ss) {
    pb$tick()
    x <- sapply(c(1:C),samplefc)
    y_new <- rnorm(1,mean(x),sqrt(temp/C))
    pc <- mapply(po_est_bound_mean,x=x,c=c(1:C),MoreArgs=list(y=y_new,N=N,temp=temp))
    g_new = exp(-0.5*(C-1)*var(x)/temp) * prod(pc)
    if (runif(1) < g_new/g_current) {
      y_current <- y_new
      g_current <- g_new
    }
    sample <- append(sample,y_current)
  }
  pb$terminate()
  return(sample)
}

pmbf <- function(ss,phi,find_bound,temp,n_partition,N,y0,g0) {
  sample <- c()
  y_current <- y0
  g_current <- g0
  pb <- progress_bar$new(
    format = "sampling [:bar] :current/:total (:percent) eta: :eta, elapsed: :elapsed",
    clear = FALSE, total = ss, width = 80)
  while (length(sample)<ss) {
    pb$update(length(sample)/ss)
    x <- sapply(c(1:C),samplefc)
    part_time <- seq(from=0,to=temp,length.out=n_partition)
    x_part_old <- x
    pc <- 1
    for (i in 2:length(part_time)){
      x_part_new = mvrnorm(
        mu=(x_part_old*(temp-part_time[i])/(temp-part_time[i-1]) + mean(x_part_old)*(part_time[i]-part_time[i-1])/(temp-part_time[i-1])),
        Sigma=matrix(rep((part_time[i]-part_time[i-1])^2/(C*(temp-part_time[i-1])),C^2),nrow=C)
      )
      for (c in (1:C)){
        pc <- pc * po_est_bound_mean(x=x_part_old[c],c=c,y=x_part_new[c],N=N,temp=(part_time[i]-part_time[i-1]))
      }
      x_part_old <- x_part_new
    }
    g_new = exp(-0.5*(C-1)*var(x)/temp) * pc
    if (runif(1) < g_new/g_current) {
      y_current <- x_part_new[1]
      g_current <- g_new
    }
    sample <- append(sample,y_current)
  }
  pb$terminate()
  return(sample)
}

dev.new()
par(mfrow = c(3, 2))

# for (i in c(1,2,3,4,5,6)) {
#   mu = c(0.5*6,-0.5*6)
#   C = 2
#   N=16
#   n_partition=4
#   print(i)
#   a <- Sys.time()
#   z<- pmbf(ss=25000,phi=phi,find_bound=find_bound,temp=2^(i-1),n_partition=n_partition,N=N,y0=0,g0=1e-20)
#   b <- Sys.time()
#   den_z <- density(z)
#   error <- mean(abs(f(den_z$x)-den_z$y))
#   pmf_times = append(pmf_times,round(log(as.numeric(difftime(time1 = b, time2 = a, units = "secs"))),1))
#   curve(f,xlim = c(-3,3),ylim=c(0,0.8),main=glue('mean diff:{6},N:{N},log time:{round(log(as.numeric(difftime(time1 = b, time2 = a, units = "secs"))),1)},error:{round(error,3)},n_par={n_partition}'))
#   hist(z[5001:25000],freq=F,add=T,col="red",nclass=50)
#   print('-----')
# }

mu_diff <- 1
mu = c(0.5*mu_diff,-0.5*mu_diff)

for (i in c(1,2,3,4,5,6)) {
  z <- pomf(
    ss=10^i,
    phi=phi,
    find_bound = find_bound,
    temp = 1,
    N=1,
    y0=0,
    g0=1e-20
  )
  den_z <- density(z)
  error <- mean(abs(f(den_z$x)-den_z$y))
  curve(f,xlim = c(-3,3),ylim=c(0,0.8),main=glue('diff={mu_diff},N={1},temp={1},ss={10^i}: error:{error}'))
  hist(z[(10^i/2+1):(10^i)],freq=F,add=T,col="red",nclass=50)
}


# z2<- pmbf(ss=10000,phi=phi,find_bound=find_bound,temp=2,n_partition=4,N=5,y0=0,g0=1e-20)
# den_z2 <- density(z2)
# error2 <- mean(abs(f(den_z2$x)-den_z2$y))
# curve(f,xlim = c(-3,3),ylim=c(0,0.8),main=glue('diff={mu_diff},N={1},temp={1},ss={10^6}: error:{error2}'))
# hist(z2,freq=F,add=T,col="red",nclass=50)


