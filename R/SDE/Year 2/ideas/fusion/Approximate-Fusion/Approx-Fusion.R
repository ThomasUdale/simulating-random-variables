library(progress)
library(layeredBB)
library(glue)
library(ggplot2)
library(mvnfast)
library(ggpmisc)
# devtools::load_all("~/Documents/simulating-random-variables/R-packages/bm.rust")

f <- function(x) {
  dnorm(x,0,sqrt(0.5))
}

mu_0 = 2
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


ea3 <- function(x,y,t,c) {
  bes_layer <- bessel_layer_simulation(x = x, y = y, s = 0, t = t, mult = 0.2)
  rbound <- find_bound(bes_layer,c)
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

mcf <- function(t) {
  while (TRUE) {
    x <- sapply(c(1:C),samplefc)
    xbar <- mean(x)
    xvar <- var(x) * (C-1)/C
    y = rnorm(1,xbar,sqrt(t/C))
    if (log(runif(1)) < -C*xvar/(2*t)) {
      if (all(mapply(ea3,x,c=c(1:C),MoreArgs=list(y=y,t=t)))) {
        return(y)
      }
    }
  }
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
  sum(log(pmax((rbound - phi(layered_bb$simulated_path[1,],c))/rbound,0)))
}

bf <- function(n,t,m) {
  x <- t(replicate(n,sapply(c(1:C),samplefc)))
  w <- apply(x,1,function(xc) {exp(-0.5*(C-1)*var(xc)/t)})
  
  for (j in 1:m) {
    print(j)
    
    x <- x[sample(n,n,replace=T,prob=w),]
    w <- rep(1/n,n)
    
    cov <- diag(rep(t/m*1*(m-j)/(m-j+1),C))+
      matrix(t/m*1/(C*(m-j+1)),nrow=C,ncol=C)
    for (i in 1:n){
      if (j<m){
        new_x <- rmvn(
          1,
          x[i,]*(m-j)/(m-j+1) + mean(x[i,])*(1/(m-j+1)),
          cov
        )
      } else {
        new_x <- rmvnorm(
          1,
          x[i,]*(m-j)/(m-j+1) + mean(x[i,])*(1/(m-j+1)),
          cov
        )
      }
      
      w[i] <- w[i]*exp(
        sum(
          mapply(
            po_est_bound,
            x=x[i,],
            y=new_x,
            c=c(1:C),
            MoreArgs=list(temp=t/m)
          )
        )
      )
      x[i,] <- new_x
    }
  } 
  return(list(x=x,w=w))
}

approx_integral <- function(x,y,s,t,eps,c) {
  cur_length = 1
  cur_eps = eps+1
  phis = c(phi(x,c),phi(y,c))
  cur_I = (t-s)/(cur_length+1)*sum(phis)
  times <- (0:cur_length)/cur_length*(t-s) + s
  path = rbind(brownian_bridge(x,y,s,t,times),times)
  while(cur_eps>eps){
    new_points <- matrix(nrow=2,ncol=cur_length)
    for (i in 1:cur_length) {
      # sim new point
      new_points[,i] <- c(brownian_bridge(
        path[1,i],
        path[1,i+1],
        path[2,i],
        path[2,i+1],
        (path[2,i]+path[2,i+1])/2
      ),(path[2,i]+path[2,i+1])/2)
    }
    phis <- append(phis,phi(new_points[1,],c))
    new_path <- matrix(nrow=2,ncol=(cur_length*2+1))
    new_path[,(1:(cur_length*2+1))[(1:(cur_length*2+1))%%2==1]] = path
    new_path[,(1:(cur_length*2+1))[(1:(cur_length*2+1))%%2==0]] = new_points
    path <- new_path
    cur_length = cur_length*2
    new_I <- 0.5*(t-s)/cur_length*(2*sum(phis)-(phis[1]+phis[2]))
    cur_eps <- abs(new_I-cur_I)
    cur_I <- new_I
  }
  return(cur_I)
}



abf <- function(n,t,m,eps) {
  x <- t(replicate(n,sapply(c(1:C),samplefc)))
  w <- apply(x,1,function(xc) {exp(-0.5*(C-1)*var(xc)/t)})
  
  for (j in 1:m) {
    print(j)
    
    x <- x[sample(n,n,replace=T,prob=w),]
    w <- rep(1/n,n)
    
    cov <- diag(rep(t/m*1*(m-j)/(m-j+1),C))+
      matrix(t/m*1/(C*(m-j+1)),nrow=C,ncol=C)
    for (i in 1:n){
      if (j<m){
        new_x <- rmvn(
          1,
          x[i,]*(m-j)/(m-j+1) + mean(x[i,])*(1/(m-j+1)),
          cov
        )
      } else {
        new_x <- rmvnorm(
          1,
          x[i,]*(m-j)/(m-j+1) + mean(x[i,])*(1/(m-j+1)),
          cov,
          checkSymmetry=FALSE
        )
      }
      
      w[i] <- w[i]*exp(-
        sum(
          mapply(
            approx_integral,
            x=x[i,],
            y=new_x,
            c=c(1:C),
            MoreArgs=list(
              s=0,
              t=t/m,
              eps=eps
            )
          )
        )
      )
      x[i,] <- new_x
    }
  } 
  return(list(x=x,w=w))
}

iad <- function(x,w=NaN){
  if (any(is.nan(w))) {
    dens <- density(x,from=-4,to=4,n=1024)
  } else {
    dens <- density(x,weights=w,from=-4,to=4,n=1024)
  }
  sum(abs(dens$y-f(seq(-4,4,length.out=1024))))*8/1024
}

##### Running Simulations

N=1e4
t=2
# mcf_out <- pbapply::pbreplicate(N,mcf(t))
# bf_out <- bf(N,t,10)
abf_out <- abf(N,t,5,0.05)


##### Graphs

iad_df <- data.frame(rbind(
  c("mcf",iad(mcf_out)),
  c("bf",iad(bf_out$x[,1],bf_out$w/sum(bf_out$w))),
  c("abf",iad(abf_out$x[,1],abf_out$w/sum(abf_out$w)))
))
colnames(iad_df) <- c("algo","IAD")
iad_df$IAD <- round(as.numeric(iad_df$IAD),4)


plot <- ggplot()+
  xlim(mu[2]-3,mu[1]+3)+
  geom_function(fun=f,aes(color='fusion-density'),key_glyph='path')+
  geom_function(fun=fc,args=list(c=1),aes(color='sub-post'),alpha=0.4,key_glyph='path')+
  geom_function(fun=fc,args=list(c=2),aes(color='sub-post'),alpha=0.4,key_glyph='path')+
  geom_density(aes(x=mcf_out,color='mcf'),alpha=0.3,key_glyph='path')+
  geom_density(aes(x=bf_out$x[,1],weight=bf_out$w,color='bf'),alpha=0.3,key_glyph='path')+
  geom_density(aes(x=abf_out$x[,1],weight=abf_out$w,color='abf'),alpha=0.3,key_glyph='path')+
  scale_color_manual(
    values=c(
      'sub-post'='green',
      'fusion-density'='black',
      'mcf'='red',
      'bf'='blue',
      'abf'='purple'
    )
  )+
  guides(color = guide_legend(title = "Densities"))+
  annotate(
    geom='table',
    x=4,
    y=0.6,
    label=list(iad_df)
  )
print(plot)

