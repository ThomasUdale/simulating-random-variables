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


snis_single<- function(i) {
  j <- sample.int(ss,C,replace=T)
  sj <- do.call(cbind, lapply(c(1:C),function(c){sub_post_sample[j[c],,c]}))
  mu <- rowMeans(sj)
  sigma <- diag((apply(sj,1,sd)^2))
  y <- rmvnorm(1,mu,sigma)
  w <- sum(unlist(lapply(c(1:C),log_likelihood,beta=y))) - log(dmvnorm(y,mu,sigma))
  return(c(y,exp(w)))
}

snis_par<-function(n) {
  out <- pbmcapply::pbmclapply(c(1:n),snis_single,mc.cores=12)
  do.call(rbind,out)
}

calc_iad <- function(d) {
  iad <- 0
  for (var_id in c(1:(d+1))){
    iad <- iad + sum(abs(density(out[,var_id],weights=out[,d+2]/sum(out[,d+2]),from=-10,to=10,n=2000)$y - density(full_post_sample[,var_id],from=-10,to=10,n=2000)$y))/2000*20
  }
  iad/(d+1)
}

resample <- function(s,d,m){
  remaining_weight <- sum(s[,d+1])
  remaining_ids <- c(1:nrow(s))
  running_max_w <- max(s[remaining_ids,d+1])
  out_s <- c()
  pb <- progress_bar$new(
    format = "sampling [:bar] :current/:total (:percent) eta: :eta, elapsed: :elapsed",
    clear = FALSE, total = (m), width = 80)
  for (k in 1:m){
    pb$tick()
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
  pb$terminate()
  return(out_s)
}

C = 40

set.seed(1)

d=3
true_beta = c(-0,rnorm(d,1,2))
print(true_beta)

ss = 10000
burn = 2000

m = 1000

x <- matrix(rnorm(m*d,1,1),ncol=d)
X <- cbind(rep(1,m),x)
y <- rbinom(m,1,p(X,true_beta))
print(sum(y))

full_post_sample <- log_reg_posterior(X,y,ss,burn,C,d+1)
dimnames(full_post_sample) <- list(c(1:ss),c(1:(d+1)))
full_post_sample_df <- melt(full_post_sample)
sub_post_sample <- simplify2array(pbmcapply::pbmclapply(c(1:C),sub_post,mc.cores=8,d=(d+1)))
dimnames(sub_post_sample) <- list(c(1:ss),c(1:(d+1)),c(1:C))
sub_post_sample_df <- melt(sub_post_sample)

for (n in c(1e4)) {
  out <- snis_par(n)
  dimnames(out) <- list(c(1:n),c(1:(d+2)))
  out_resample <- melt(resample(out[out[,d+2]!=0,],d+1,1e3))
  out_df <- melt(out[,c(1:(d+1))])
  out_df$weight <- rep(out[,d+2],d+1)
  print(calc_iad(d))

  l <- vector("list",d+1)

  for (pltid in c(1:(d+1))){
    l[[pltid]] <- ggplot()+
      xlim(true_beta[pltid]-5,true_beta[pltid]+5)+
      geom_density(data=subset(sub_post_sample_df,Var2==pltid),aes(x=value,group=Var3),col='blue',alpha=0.05)+
      geom_density(data=subset(full_post_sample_df,Var2==pltid),aes(x=value),col='black',alpha=0.8)+
      geom_density(data=subset(out_df,Var2==pltid),aes(x=value,weight=weight),col='red',alpha=0.4)+
      geom_density(data=subset(out_resample,Var2==pltid),aes(x=value),col='orange',alpha=0.4)+
      geom_vline(xintercept=true_beta[pltid],col='green')
  }

  print(wrap_plots(l))

}


  



