library('MASS')
library('progress')
library(layeredBB)
library(parallel)
library(ggplot2)
library(patchwork)
library(reshape2)

C = 10

true_beta = c(-4,-2)

ss = 10000
burn = 2000

p <- function(X,beta) {
  z <- exp(X %*% beta)
  z/(1+z)
}


log_likelihood_full <- function(beta){
  z <- X %*% t(beta)
  return(sum(z*y - log(1+exp(z))) - sum((beta)^2/(2*C*10)))
}


log_likelihood <- function(beta,c){
  z <- X[c((m/C*(c-1)+1):(m/C*c)),] %*% t(beta)
  return(sum(z*y[c((m/C*(c-1)+1):(m/C*c))] - log(1+exp(z))) - sum((beta)^2/(2*C*10)))
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

phi <- function(omega_beta,c) {
  y <- y[c((m/C*(c-1)+1):(m/C*c))]
  x <- X[c((m/C*(c-1)+1):(m/C*c)),]
  z <- exp(x %*% omega_beta)
  p <- as.vector(z/(1+z))
  alpha <- t(x) %*% (y-p)-1/(C*10)*omega_beta
  div <- -t(x) %*% diag(p*(1-p)) %*% x - 1/(C*10)
  0.5*(sum(alpha^2)+sum(diag(div)))
}

approx_integral <- function(x,y,s,t,eps,c,d) {
  cur_length = 1
  cur_eps = eps+1
  phis = c(phi(x,c),phi(y,c))
  cur_I = (t-s)/(cur_length+1)*sum(phis)
  times <- (0:cur_length)/cur_length*(t-s) + s
  path <- matrix(nrow=d+1,ncol=cur_length+1)
  path[d+1,] <- times
  for (d_theta in 1:d){
    path[d_theta,] <- brownian_bridge(x[d_theta],y[d_theta],s,t,times)
  }
  while(cur_eps>eps) {
    new_points <- matrix(nrow=d+1,ncol=cur_length)
    for (i in 1:cur_length) {
      for (d_theta in 1:d) {
        new_points[d_theta,i] <- brownian_bridge(
          path[d_theta,i],
          path[d_theta,i+1],
          path[d+1,i],
          path[d+1,i+1],
          (path[d+1,i]+path[d+1,i+1])/2
        )
      }
      new_points[d+1,i] <- (path[d+1,i]+path[d+1,i+1])/2
    }
    phis <- append(phis,apply(new_points[1:d,,drop=F],2,phi,c=c))
    new_path <- matrix(nrow=d+1,ncol=(cur_length*2+1))
    new_path[,(1:(cur_length*2+1))[(1:(cur_length*2+1))%%2==1]] = path
    new_path[,(1:(cur_length*2+1))[(1:(cur_length*2+1))%%2==0]] = new_points
    path <- new_path
    cur_length = cur_length*2
    new_I <- (t-s)/(cur_length+1)*sum(phis)
    cur_eps <- abs(new_I-cur_I)
    cur_I <- new_I
  }
  return(cur_I)
}

abf <- function(n,t,m,d) {
  x <- array(dim=c(n,C,d))
  j <- sample.int(ss,C*n,replace=T)
  for (i in 0:(n-1)){
    x[i+1,,] <- t(do.call(cbind, lapply(c(1:C),function(c){sub_post_sample[j[c+i*C],,c]})))
  }
  
  w <- apply(0.5*(C-1)*apply(x,c(1,3),var)/t,1,sum)
  
  for (j in 1:m) {
    print(j)
    x <- x[sample(n,n,replace=T,prob=exp(-w)),,]
    w <- rep(0,n)
    
    cov <- diag(rep(t/m*1*(m-j)/(m-j+1),C*d))+
          matrix(t/m*1/(C*(m-j+1)),nrow=C*d,ncol=C*d)
    for (i in 1:n) {
      if (j<m){
        new_x <- matrix(rmvn(
          1,
          as.vector(x[i,,])*(m-j)/(m-j+1) + rep(colMeans(x[i,,]),each=C)*(1/(m-j+1)),
          cov
        ),ncol=d)
      } else {
        new_x <- matrix(rmvnorm(
          1,
          as.vector(x[i,,])*(m-j)/(m-j+1) + rep(colMeans(x[i,,]),each=C)*(1/(m-j+1)),
          cov
        ),ncol=d)
      }
      for (c in 1:C) {
        w[i] <- w[i]+approx_integral(
          x[i,c,],
          new_x[c,],
          s=0,
          t=t/m,
          eps=10,
          c=c,
          d=d
        )
      }
      x[i,,] <- new_x 
    }
  }
  return(list(x=x[,1,],w=exp(-w)))
}

# m = 1000
# x <- rnorm(m,0.7,1)
# X <- cbind(rep(1,m),x)
# y <- rbinom(m,1,p(X,true_beta))
# print(sum(y))
# 
# full_post_sample <- log_reg_posterior(X,y,ss,burn,C)
# sub_post_sample <- simplify2array(pbmcapply::pbmclapply(c(1:C),sub_post,mc.cores=4))
# dimnames(sub_post_sample) <- list(c(1:ss),c(1:2),c(1:C))
# sub_post_sample_df <- melt(sub_post_sample)
# 
abf_out <- abf(1e4,1,5,2)

plot1 <- ggplot()+
  xlim(-6,-1)+
  geom_density(data=subset(sub_post_sample_df,Var2==1),aes(x=value,group=Var3),col='blue',alpha=0.05)+
  geom_density(aes(x=full_post_sample[,1]),col='black',alpha=0.8)+
  geom_density(aes(x=abf_out$x[,1],weight=abf_out$w),col='red',alpha=0.5)

plot2 <- ggplot()+
  xlim(-4,1)+
  geom_density(data=subset(sub_post_sample_df,Var2==2),aes(x=value,group=Var3),col='blue',alpha=0.05)+
  geom_density(aes(x=full_post_sample[,2]),col='black',alpha=0.8)+
  geom_density(aes(x=abf_out$x[,2],weight=abf_out$w),col='red',alpha=0.5)


print(plot1+plot2)



