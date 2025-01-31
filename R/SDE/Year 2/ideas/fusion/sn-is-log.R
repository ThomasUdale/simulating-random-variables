library('MASS')
library('progress')
library(layeredBB)
library(parallel)
library(ggplot2)
library(patchwork)
library(reshape2)

C = 40

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




snis <- function(n) {
  
  pb <- progress_bar$new(
    format = "sampling [:bar] :current/:total (:percent) eta: :eta, elapsed: :elapsed",
    clear = FALSE, total = (n), width = 80)
  
  out <- c()
  
  for(i in 1:n){
    pb$tick()
    j <- sample.int(ss,C,replace=T)
    sj <- do.call(cbind, lapply(c(1:C),function(c){sub_post_sample[j[c],,c]}))
    mu <- rowMeans(sj)
    sigma <- diag(2)
    y <- rmvnorm(1,mu,sigma)
    w <- sum(unlist(lapply(c(1:C),log_likelihood,beta=y))) - log(dmvnorm(y,mu,sigma))
    out <- rbind(out,c(y,exp(w)))
  }
  pb$terminate()
  return(out)
}

m = 1000
x <- rnorm(m,0.7,1)
X <- cbind(rep(1,m),x)
y <- rbinom(m,1,p(X,true_beta))
print(sum(y))

full_post_sample <- log_reg_posterior(X,y,ss,burn,C)
sub_post_sample <- simplify2array(pbmcapply::pbmclapply(c(1:C),sub_post,mc.cores=4))
dimnames(sub_post_sample) <- list(c(1:ss),c(1:2),c(1:C))
sub_post_sample_df <- melt(sub_post_sample)


out1 <- snis(1e4)
out2 <- snis(1e5)
out3 <- snis(1e6)

plot1 <- ggplot()+
  xlim(-5,0)+
  geom_density(data=subset(sub_post_sample_df,Var2==1),aes(x=value,group=Var3),col='blue',alpha=0.05)+
  geom_density(aes(x=full_post_sample[,1]),col='black',alpha=0.8)+
  geom_density(aes(x=out1[,1],weight=out1[,3]),col='red')+
  geom_density(aes(x=out2[,1],weight=out2[,3]),col='orange')+
  geom_density(aes(x=out3[,1],weight=out3[,3]),col='green')
  

plot2 <- ggplot()+
  xlim(-3,1.5)+
  geom_density(data=subset(sub_post_sample_df,Var2==2),aes(x=value,group=Var3),col='blue',alpha=0.05)+
  geom_density(aes(x=full_post_sample[,2]),col='black',alpha=0.8)+
  geom_density(aes(x=out1[,2],weight=out1[,3]),col='red')+
  geom_density(aes(x=out2[,2],weight=out2[,3]),col='orange')+
  geom_density(aes(x=out3[,2],weight=out3[,3]),col='green')

print(plot1+plot2)



