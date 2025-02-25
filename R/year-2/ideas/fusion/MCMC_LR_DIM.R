library('MASS')
library('progress')
library(parallel)
library(ggplot2)
library(patchwork)
library(reshape2)
library(mvtnorm)



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

log_reg_posterior <- function(X,y,ss,burn,C,d) {
  kappa = y-0.5
  true_post <- c()
  omega <- rep(0,nrow(X))
  beta <- rep(1,d)
  pb <- progress_bar$new(
    format = "sampling [:bar] :current/:total (:percent) eta: :eta, elapsed: :elapsed",
    clear = FALSE, total = (ss+burn), width = 80)
  for (i in 1:(ss+burn)) {
    pb$tick()
    omega <- BayesLogit::rpg(nrow(X),1,X %*% beta)
    var <-  solve(t(X) %*% diag(omega) %*% X + solve(diag(rep(C*10,d))))
    mu <- var %*% (t(X) %*% kappa)
    beta <- mvrnorm(n = 1, mu, var)
    if (i > burn) {
      true_post <- rbind(true_post,beta)
    }
  }
  pb$terminate()
  true_post
}

sub_post <- function(c,d) {
  log_reg_posterior(X[c((m/C*(c-1)+1):(m/C*c)),],y[c((m/C*(c-1)+1):(m/C*c))],ss,burn,C,d)
}

mcmc <- function(burn_mcmc,ss_mcmc,d){
  pb <- progress_bar$new(
    format = "sampling [:bar] :current/:total (:percent) eta: :eta, elapsed: :elapsed",
    clear = FALSE, total = (ss_mcmc+burn_mcmc), width = 80)
  y <- rep(0,d+1)
  j <- sample.int(ss,C,replace=T)
  xc <- do.call(cbind, lapply(c(1:C),function(c){sub_post_sample[j[c],,c]}))
  mean_xc <- rowMeans(xc)
  sd_xc <- exp(-apply(xc,1,sd))
  pi_y <- exp(sum(unlist(lapply(c(1:C),log_likelihood,beta=t(y)))))
  out <- matrix(,nrow=ss_mcmc,ncol=d+1)
  for (t in 1:(burn_mcmc+ss_mcmc)){
    pb$tick()
    j <- sample.int(ss,C,replace=T)
    new_xc <- do.call(cbind, lapply(c(1:C),function(c){sub_post_sample[j[c],,c]}))
    sd_new_xc <- exp(-apply(new_xc,1,sd))
    mean_new_xc <- rowMeans(new_xc)
    
    new_y <-rmvnorm(1,mean_new_xc*sd_new_xc+y*(1-sd_new_xc),diag((1-sd_new_xc)/4))
    pi_new_y <- exp(sum(unlist(lapply(c(1:C),log_likelihood,beta=new_y))))
    
    r <- pi_new_y/pi_y * dmvnorm(y,mean_xc*sd_xc+new_y*(1-sd_xc),diag((1-sd_xc)/4))/dmvnorm(new_y,mean_new_xc*sd_new_xc+y*(1-sd_new_xc),diag((1-sd_new_xc)/4))
    if (runif(1)<min(r,1)){
      xc <- new_xc
      mean_xc <- mean_new_xc
      sd_xc <- sd_new_xc
      y <- new_y
      pi_y <- pi_new_y
    }
    if (t>burn_mcmc){
      out[t-burn_mcmc,] <- y
    }
  }
  pb$terminate()
  out
}

C = 40

set.seed(1)

# d=8
# true_beta = c(-2,rnorm(d,-1,2))
# print(true_beta)
# 
# ss = 10000
# burn = 2000
# 
# m = 1000
# 
# x <- matrix(rnorm(m*d,1,1),ncol=d)
# X <- cbind(rep(1,m),x)
# y <- rbinom(m,1,p(X,true_beta))
# print(sum(y))
# 
# full_post_sample <- log_reg_posterior(X,y,ss,burn,C,d+1)
# dimnames(full_post_sample) <- list(c(1:ss),c(1:(d+1)))
# full_post_sample_df <- melt(full_post_sample)
# sub_post_sample <- simplify2array(pbmcapply::pbmclapply(c(1:C),sub_post,mc.cores=8,d=(d+1)))
# dimnames(sub_post_sample) <- list(c(1:ss),c(1:(d+1)),c(1:C))
# sub_post_sample_df <- melt(sub_post_sample)

for (n in c(1e6)) {
  out <- mcmc(1e4,n,d)
  dimnames(out) <- list(c(1:n),c(1:(d+1)))
  out_df <- melt(out[,c(1:(d+1))])

  calc_iad <- function() {
    iad <- 0
    for (var_id in c(1:(d+1))){
      iad <- iad + sum(abs(density(out[,var_id],from=-10,to=10,n=2000)$y - density(full_post_sample[,var_id],from=-10,to=10,n=2000)$y))/2000*20
    }
    iad/(d+1)
  }

  print(calc_iad())
  print(c(length(unique(out[,1])),length(unique(out[,1]))/n))

  l <- vector("list",d+1)

  for (pltid in c(1:(d+1))){
    l[[pltid]] <- ggplot()+
      xlim(true_beta[pltid]-2,true_beta[pltid]+2)+
      ylim(0,6)+
      geom_density(data=subset(sub_post_sample_df,Var2==pltid),aes(x=value,group=Var3),col='blue',alpha=0.05)+
      geom_density(data=subset(full_post_sample_df,Var2==pltid),aes(x=value),col='black',alpha=0.8)+
      geom_density(data=subset(out_df,Var2==pltid),aes(x=value),col='red',alpha=0.4)+
      geom_vline(xintercept=true_beta[pltid],col='green')
  }

  print(wrap_plots(l))

}



  



