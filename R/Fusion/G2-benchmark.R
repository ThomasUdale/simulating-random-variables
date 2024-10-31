set.seed(3)

library('MASS')
library('progress')
library(layeredBB)
library(parallel)

C = 2

true_beta = c(-1,0.5)

ss = 2000
burn = 2000

p <- function(X,beta) {
  z <- exp(X %*% beta)
  z/(1+z)
}



log_reg_posterior <- function(X,y,ss,burn,C) {
  kappa = y-0.5
  true_post <- c()
  omega <- rep(0,nrow(X))
  beta <- c(1,1)
  pb <- progress_bar$new(
    format = "sampling [:bar] :current/:total (:percent) eta: :eta, elapsed: :elapsed",
    clear = FALSE, total = (ss+burn), width = 80)
  for (i in 1:(ss+burn)) {
    pb$tick()
    omega <- BayesLogit::rpg(nrow(X),1,beta[1] + beta[2]*x)
    var <-  solve(t(X) %*% diag(omega) %*% X + solve(diag(rep(C*10,2))))
    mu <- var %*% (t(X) %*% kappa)
    beta <- mvrnorm(n = 1, mu, var)
    if (i > burn) {
      true_post <- rbind(true_post,beta)
    }
  }
  pb$terminate()
  true_post
}

sub_post <- function(c) {
  log_reg_posterior(X[c((m/C*(c-1)+1):(m/C*c)),],y[c((m/C*(c-1)+1):(m/C*c))],ss,burn,C)
}

# m = 1000
# x <- rnorm(m,0.2,1)
# X <- cbind(rep(1,m),x)
# y <- rbinom(m,1,p(X,true_beta))
# 
# print(sum(y))
# 
# true_post <- log_reg_posterior(X,y,ss,burn,1)
# true_dens1 <- density(true_post[,1])
# true_dens2 <- density(true_post[,2])
# 
# sub_posts <- rep(NaN, C*ss*2)
# sub_posts <- array(sub_posts,dim=c(C,ss,2))
# for (c in 1:C) {
#   sub_posts[c,,] = sub_post(c)
# }
# 
# par(mfrow=c(1, 2))
# plot(true_dens1)
# for (c in 1:C) {
#   lines(density(sub_posts[c,,1]),col=rgb(1,0,0, 0.5))
# }
# plot(true_dens2)
# for (c in 1:C) {
#   lines(density(sub_posts[c,,2]),col=rgb(1,0,0, 0.5))
# }

samplefc <- function(c) {
  total <- length(sub_posts[c,,1])
  sub_posts[c,sample(total,1),]
}

phi <- function(omega_beta,c) {
  y <- y[c((m/C*(c-1)+1):(m/C*c))]
  x <- X[c((m/C*(c-1)+1):(m/C*c)),]
  z <- exp(omega_beta[1] + x*omega_beta[2])
  grad <- (sum(y - z/(1+z)) - 1/(C*10)*omega_beta[1])^2 + (sum((y - z/(1+z))*x) - 1/(C*10)*omega_beta[2])^2
  div <- -sum((z/(z+1)^2)*(1+x^2)) - 1/(C*10)
  0.5*(grad+div)
}

# global_lb <- sapply(c(1:C),function(c) NMOF::gridSearch(phi,lower=c(-100,-100),upper=c(30,30),npar=2,c=c,printDetail=F,n=30)$minfun)

find_lower_bound_approx <- function(bes_layer,c) {
  lower <- unlist(lapply(bes_layer, function(x) x$L))
  upper <- unlist(lapply(bes_layer, function(x) x$U))
  NMOF::gridSearch(phi,lower=lower,upper=upper,npar=2,c=c,printDetail=F)$minfun
}

find_upper_bound_approx <- function(bes_layer,c) {
  
  opt_func <- function(omega_beta,c){
    -phi(omega_beta,c)
  }

  lower <- unlist(lapply(bes_layer, function(x) x$L))
  upper <- unlist(lapply(bes_layer, function(x) x$U))
  
  -NMOF::gridSearch(opt_func,lower=lower,upper=upper,npar=2,c=c,printDetail=F,n=5)$minfun
}

po_est_bound <- function(start_x,end_y,temp,c) {
  bes_layer <- multi_bessel_layer_simulation(dim=2,x = start_x, y = end_y, s = 0, t = temp, mult = 0.2)
  lower_bound <- find_lower_bound_approx(bes_layer,c)
  upper_bound <- find_upper_bound_approx(bes_layer,c) - lower_bound
  k <- rpois(1,upper_bound*temp)
  if (is.na(k)){
    print(upper_bound)
    print(lower_bound)
    print(bes_layer)
  }
  if (k==0) {
    return(1)
  }
  skeleton = runif(k,0,temp)

  layered_bb <- multi_layered_brownian_bridge(dim=2,x = start_x,
                                        y = end_y,
                                        s = 0,
                                        t = temp,
                                        bessel_layer = bes_layer,
                                        times = skeleton)
  if (prod((upper_bound - (apply(layered_bb$simulated_path,2,phi,c=c)-lower_bound))/upper_bound)<0){
    print("---")
    print(k)
    print(c)
    print(upper_bound)
    print(apply(layered_bb$simulated_path,2,phi,c=c)-lower_bound)
    print(layered_bb$simulated_path[c(1:2),])
    print(mapply(phi,asplit(layered_bb$simulated_path[c(1:2),],2),MoreArgs=list(c=c)))
    print(bes_layer)
  }

  prod((upper_bound - (apply(layered_bb$simulated_path,2,phi,c=c)-lower_bound))/upper_bound)
}

po_est_bound_mean <- function(start_x,c,end_y,N,temp) {
  mean(replicate(N,po_est_bound(start_x,end_y,temp,c)))
}

pomf <- function(ss,temp,y0,g0) {
  dev.new()
  sample <- array(rep(NaN, ss*2),dim=c(ss,2))
  tick <- 0
  y_current <- y0
  g_current <- g0
  pb <- progress_bar$new(
    format = "sampling [:bar] :current/:total (:percent) eta: :eta, elapsed: :elapsed",
    clear = FALSE, total = ss, width = 80)
  while (tick<ss) {
    tick <- tick + 1
    pb$tick()
    x <- sapply(c(1:C),samplefc)
    y_new <- mvrnorm(1,mu=rowMeans(x),sigma=diag(2)*(temp/C))
    pc <- mapply(po_est_bound,start_x=asplit(x,2),c=c(1:C),MoreArgs=list(end_y=y_new,temp=temp))
    g_new = exp(-0.5*sum((x-rowMeans(x))^2)/temp) * prod(pc)
    if (runif(1) < g_new/g_current) {
      y_current <- y_new
      g_current <- g_new
    }
    sample[tick,] <- y_current
    
    if (tick%%100==0){
      par(mfrow=c(1, 2))
      plot(true_dens1,xlim=c(-8,1))
      hist(sample[c(1:tick),1],freq=F,add=T)
      plot(true_dens2,xlim=c(-8,1))
      hist(sample[c(1:tick),2],freq=F,add=T)
    }
  }
  pb$terminate()
  return(sample)
}

pobf <- function(ss,dims,N,temp,n_partition,y0,g0) {
  
  sample <- array(rep(NaN, ss*dims),dim=c(ss,dims))
  tick <- 0
  y_current <- y0
  g_current <- g0
  exp_current <- 1e-200
  
  pb <- progress_bar$new(
    format = "sampling [:bar] :current/:total (:percent) eta: :eta, elapsed: :elapsed",
    clear = FALSE, total = ss, width = 80)
  while (tick<ss) {
    tick <- tick + 1
    pb$tick()
    print(1)
    x <- sapply(c(1:C),samplefc)
    part_time <- seq(from=0,to=temp,length.out=n_partition)
    # x_part <- array(rep(NaN, length(part_time*dims),dim=c(ss,dims))
    x_part_old <- c(x)
    pc <- 0 
    
    for (i in 2:length(part_time)){
      x_part_new = mvrnorm(
        mu=(x_part_old*(temp-part_time[i])/(temp-part_time[i-1]) + rep(rowMeans(matrix(x_part_old,nrow=dims)),C)*(part_time[i]-part_time[i-1])/(temp-part_time[i-1])),
        Sigma=matrix(rep((part_time[i]-part_time[i-1])^2/(C*(temp-part_time[i-1])),(C*dims)^2),nrow=(C*dims))
      )
      for (c in (1:C)){
        pc <- pc + log(po_est_bound_mean(start_x=matrix(x_part_old,nrow=dims)[,c],c=c,end_y=matrix(x_part_new,nrow=dims)[,c],temp=(part_time[i]-part_time[i-1]),N=N))
      }
      x_part_old <- x_part_new
    }
    
    g_new = -0.5*sum((x-x_part_new[c(1:dims)])^2)/temp + pc
    print(g_new)
    if (tick==1){
      y_current <- x_part_new[c(1:dims)]
      g_current <- g_new
    }
    if (log(runif(1)) < g_new/g_current) {
      y_current <- x_part_new[c(1:dims)]
      g_current <- g_new
    }
    sample[tick,] <- y_current
    
    if (tick%%100==0){
      par(mfrow=c(1, 2))
      plot(true_dens1,xlim=c(-8,1))
      for (c in 1:C) {
        lines(density(sub_posts[c,,1]),col=rgb(1,0,0, 0.5))
      }
      lines(density(sample[c(tick/2:tick),1]),col=rgb(0,1,0, 0.5))
      plot(true_dens2,xlim=c(-8,1))
      for (c in 1:C) {
        lines(density(sub_posts[c,,2]),col=rgb(1,0,0, 0.5))
      }
      lines(density(sample[c(tick/2:tick),2]),col=rgb(0,1,0, 0.5))
    }
  }
  pb$terminate()
  return(sample)
  
  
}

# out_s <- pomf(10000,1,c(-1,-1),1e-100)

out_bf <- pobf(20000,2,1,0.5,2,c(-1,-1),1e-200)



