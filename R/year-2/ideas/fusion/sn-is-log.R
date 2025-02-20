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

grad_log_likelihood <- function(beta,c){
  t(X[c((m/C*(c-1)+1):(m/C*c)),]) %*% (y[c((m/C*(c-1)+1):(m/C*c))] - p(X[c((m/C*(c-1)+1):(m/C*c)),],beta))
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


snis_single<- function(i) {
  j <- sample.int(ss,C,replace=T)
  sj <- do.call(cbind, lapply(c(1:C),function(c){sub_post_sample[j[c],,c]}))
  mu <- rowMeans(sj)
  # mu <- apply(sj,1,median)
  sigma <- diag((apply(sj,1,sd)^2))
  # sigma <- diag(rowMeans(abs(sapply(c(1:C),grad_log_likelihood,beta=mu)))^2)
  y <- rmvnorm(1,mu,sigma)
  w <- sum(unlist(lapply(c(1:C),log_likelihood,beta=y))) - log(dmvnorm(y,mu,sigma))
  
  # uniform proposal:
   # max_sj <- apply(sj,1,max)
   # min_sj <- apply(sj,1,min)
   # y <- runif(2,min_sj,max_sj)
   # w <- sum(unlist(lapply(c(1:C),log_likelihood,beta=t(y)))) - log(prod(dunif(y,min_sj,max_sj)))
  return(c(y,exp(w)))
}

snis_par<-function(n) {
  out <- pbmcapply::pbmclapply(c(1:n),snis_single,mc.cores=10)
  do.call(rbind,out)
}

calc_iad <- function() {
  iad <- 0
  for (var_id in c(1:2)){
    iad <- iad + sum(abs(density(out[,var_id],weights=out[,3]/sum(out[,3]),from=-10,to=10,n=2000)$y - density(full_post_sample[,var_id],from=-10,to=10,n=2000)$y))/2000*20
  }
  iad/2
}

# resampling without replacement

resample <- function(s,d,m){
  remaining_weight <- sum(s[,d+1])
  remaining_ids <- c(1:nrow(s))
  running_max_w <- max(s[remaining_ids,d+1])
  out_s <- c()
  for (k in 1:m){
    print(k)
    id_found <- FALSE
    while(!id_found){
      id <- sample(remaining_ids,1,prob=s[remaining_ids,d+1])
      id_weight <-s[id,d+1]
      if ((remaining_weight - k*running_max_w)<0){
        id_found=TRUE
      } else if (runif(1)<(remaining_weight-k*running_max_w)/(remaining_weight-k*id_weight)) {
        id_found=TRUE
      }
    }
    remaining_ids <- remaining_ids[-c(id)]
    remaining_weight = remaining_weight-s[id,d+1]
    out_s <- rbind(out_s,s[id,1:d])
  }
  return(out_s)
}

C = 40

true_beta = c(-4,-2)

ss = 10000
burn = 2000

# m = 1000
# x <- rnorm(m,0.7,1)
# X <- cbind(rep(1,m),x)
# y <- rbinom(m,1,p(X,true_beta))
# print(sum(y))
# 
# full_post_sample <- log_reg_posterior(X,y,ss,burn,1)
# sub_post_sample <- simplify2array(pbmcapply::pbmclapply(c(1:C),sub_post,mc.cores=8))
# dimnames(sub_post_sample) <- list(c(1:ss),c(1:2),c(1:C))
# sub_post_sample_df <- melt(sub_post_sample)

for (n in c(1e6)){
  # out <- rbind(out,snis_par(n))
  # out <- snis_par(n)
  # out_resample <- resample(out,2,1e2)
  print(calc_iad())

  plot1 <- ggplot()+
    xlim(-8,1)+
    geom_density(data=subset(sub_post_sample_df,Var2==1),aes(x=value,group=Var3),col='blue',alpha=0.05)+
    geom_density(aes(x=full_post_sample[,1]),col='black',alpha=0.8)+
    geom_density(aes(x=out[,1],weight=out[,3]),col='red',alpha=0.4)+
    geom_density(aes(x=out_resample[,1]),col='orange',alpha=0.4)


  plot2 <- ggplot()+
    xlim(-4,3)+
    geom_density(data=subset(sub_post_sample_df,Var2==2),aes(x=value,group=Var3),col='blue',alpha=0.05)+
    geom_density(aes(x=full_post_sample[,2]),col='black',alpha=0.8)+
    geom_density(aes(x=out[,2],weight=out[,3]),col='red',alpha=0.4)+
    geom_density(aes(x=out_resample[,2]),col='orange',alpha=0.4)

  print(plot1+plot2)
}






