
library(ggplot2)
library(progress)
library('MASS')
library(mvtnorm)

C = 2

true_beta = c(-1,0.5)

ss = 500
burn = 500

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

m = 1000
x <- rnorm(m,0.2,1)
X <- cbind(rep(1,m),x)
y <- rbinom(m,1,p(X,true_beta))
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

sample_pt <- function(x) {
  y <- rmvnorm(1,mean=x,sigma=diag(1,2))[1,]
  a <- dmvnorm(y,mean=rep(0,2),sigma=diag(10,2))/dmvnorm(x,mean=rep(0,2),sigma=diag(10,2))
  if (runif(1)<a){
    return(y)
  } else {
    return(x)
  }
}

lambda <- function(x) {
  dmvnorm(x,mean=rep(0,2),sigma=diag(10,2))
}

k_lower <-0.8

jump_sample <- function(index){
  x = sample_mu()
  i = 0
  while(TRUE){
    i <- i + 1
    p_step <- lambda(x)
    p_regen <- 0.5^(m)*mu(x)-k_lower*f(x)
    p_exit <- k_lower*f(x)
    print(c(p_step,p_regen,p_exit))
    tau <- sample(c(1,2,3),1,prob=c(p_step,p_regen,p_exit))
    
    if (tau==1){
      x <- sample_pt(x)
    } else if (tau==2){
      x <- sample_mu()
    } else {
      return(x)
    }
  }
}

output <- jump_sample()
print(output)

# samples <- pbmcapply::pbmclapply(seq(ss),jump_sample,mc.cores=12)
# samples <- do.call(rbind,samples)
# hist(samples[,1],freq=F,nclass=50,xlim=c(-1.5,-0.5))
# lines(density(samples[,1]),col='blue')
# lines(true_dens1,col='red')
