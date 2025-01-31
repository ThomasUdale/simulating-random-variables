set.seed(1)
library(ggplot2)
library(progress)
library('MASS')
library(mvtnorm)

C = 2

true_beta = c(-4,-2)

ss = 200000
burn = 10000

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

f <- function(beta) {
  prod(p(X,beta)^y*(1-p(X,beta))^(1-y))*dmvnorm(beta,mean=rep(0,2),sigma=diag(10,2))
}

mu <- function(x) {
  dmvnorm(x,mean=rep(0,2),sigma=diag(10,2))
}

sample_mu <- function() {
  rmvnorm(1,mean=rep(0,2),sigma=diag(10,2))[1,]
}

lambda_pt <- function(x) {
  y <- rmvnorm(1,mean=rep(0,2),sigma=diag(1,2))[1,]
  a <- dmvnorm(y,mean=rep(0,2),sigma=diag(1,2))/dmvnorm(x,mean=rep(0,2),sigma=diag(1,2))
  if (runif(1)<a){
    return(y)
  } else {
    return(x)
  }
}

lambda <- function(x) {
  dmvnorm(x,mean=rep(0,2),sigma=diag(1,2))
}



jump_sample <- function(ss){
  samples <- data.frame()
  cur_x = sample_mu()
  i = 0
  
  pb <- progress_bar$new(
    format = "sampling [:bar] :current/:total (:percent) eta: :eta, elapsed: :elapsed",
    clear = FALSE, total = ss, width = 80)
  
  while(nrow(samples)<ss){
    pb$update(nrow(samples)/ss)
    i <- i + 1
    f <- f(cur_x)
    tau_step <- rexp(1,lambda(cur_x)/f)
    tau_regen <- rexp(1,mu(cur_x)/f)
    tau <- min(tau_step,tau_regen)
    if (tau>0){
      samples <- rbind(samples,c(i,tau,cur_x[1],cur_x[2],tau_step<tau_regen))
    }
    
    if (tau_step<tau_regen){
      cur_x <- lambda_pt(cur_x)
    } else {
      cur_x <- sample_mu()
    }
  }
  colnames(samples) <- c('index','sample_time','beta1','beta2','step')
  return(samples)
}

output <- jump_sample(ss)
colnames(true_post) <- c('beta1','beta2')
plot <- ggplot() + geom_histogram(data=output,aes(x=beta1,y=..density..,weight=sample_time),binwidth=0.02)+geom_density(data=true_post,aes(x=beta1),col='green')
print(plot)
