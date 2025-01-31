library(MASS)
library(ggplot2)
library(progress)
library(parallel)
library(Matrix)

C = 2
M_plus = 30
M_minus = 1
N = 50000
total = 1/max(M_plus,M_minus)

f <- function(x) {
  exp(-x^4/2)/2.1558
}

fc <- function(x) {
  exp(-x^4/(2*C))
}

samplefc <- function() {
  while (TRUE) {
    x <- rnorm(1)
    if (runif(1)<(fc(x)/exp(-x^2/2+C/8))) {
      return(x)
    }
  }
}

phi <- function(x) {
  0.5*(4*x^6/C^2 - 6*x^2/C)
}

v_plus <- function(x) {
  min(max(-sum(phi(x)),0),M_minus)
}

v_minus <- function(x) {
  min(max(sum(phi(x)),0),M_plus)
}

x <- array(replicate(N*C,samplefc()),dim=c(N,C))
rho_0 <- exp(-rowSums((x - rowMeans(x))^2)/(2*total))
x = x[sample(N,N,replace=T,prob=rho_0),]

sim_forward_all <- function(x,t,new_t){
  mean <- (total-new_t)/(total-t)*x + (new_t-t)/(total-t) * rowMeans(x)
  cov <- diag((new_t-t)*(total-new_t)/(total-t),nrow=C) + matrix((new_t-t)^2/(C*(total-t)),nrow = C,ncol=C)
  if (cov[1,1]==cov[1,2]) {
    mean
  } else {
    t(chol(cov) %*% matrix(rnorm(N*C),nrow=C)) + mean
  }
  
}

t <- 0
k <- 0

pb <- progress_bar$new(
  format = "sampling [:bar] (:percent) eta: :eta, elapsed: :elapsed, k: :k",
  clear = FALSE, total = 10000, width = 80)
while (t<total){
  pb$update(t/total,tokens =list(k=k))
  k <- k+1
  tau_plus <- rexp(N,M_minus)
  tau_minus <- rexp(N,M_plus)
  tau <- cbind(tau_plus,tau_minus)
  id <- arrayInd(which.min(tau), dim(tau))
  new_t <- min(t+min(tau),total)
  x <- sim_forward_all(x,t,new_t)
  if (id[2]==1){
    if (runif(1) < v_plus(x[id[1],])/M_minus) {
      x[sample(N,1),] <- x[id[1],]
    }

  } else {
    if (runif(1) < v_minus(x[id[1],])/M_plus) {
      x[id[1],] <- x[sample(N,1),]
    }
  }
  t <- new_t
}
pb$terminate()

output <- data.frame(x[,1])
colnames(output) <- c("x")

plot <-
  ggplot()+
  geom_histogram(data=output,aes(x=x,y=..density..),binwidth=0.1)+
  xlim(-3,3)+
  geom_function(fun=f,col='red')+
  geom_function(fun=phi,col='blue')+
  ylim(-0.7,1.25)
  
print(plot)

