library(ggplot2)

mu <- 10
sigma <- 2

f <- function(x) {
  dnorm(x,mu,sqrt(sigma))
}

grad_log_f <- function(x) {
  -(x-mu)/sigma
}

pdmp <- function(t) {
  v <- sample(c(-1,1),1)
  x <- rnorm(1,mu,sqrt(sigma))
  s <- 0
  out <- c(0,x,v,0)
  while (s<t) {
    e <- rexp(1)
    jt <- rexp(1,10/t)
    if (v*(x-mu)<0){
      event_time <- abs(x-mu) + sqrt(2*e*sigma)
    } else {
      event_time <- -v*(x-mu) + sqrt((x-mu)^2 + 2*e*sigma)
    }
    if (s+min(event_time,jt)>t) {
      out <- rbind(out,c(t,x+(t-s)*v,v,0))
      return(out)
    }
    if (event_time<jt){
      x <- x + event_time*v
      v <- -v
      s <- s + event_time
      out <- rbind(out,c(s,x,v,1))
    } else {
      x <- rnorm(1,mu,sqrt(sigma))
      v <- sample(c(-1,1),1)
      s <- s+jt
      out <- rbind(out,c(s,x,v,2))
    }
  }
}

interpolate <- function(out,t,m) {
  h <- t/m
  s <- h
  id <- 1
  idx <- 1
  int_out <- matrix(nrow=m,ncol=2)
  while (idx<m) {
    if (out[id,1]>=s) {
      int_out[idx,] <- c(s,out[id-1,2]+(s-out[id-1,1])*out[id-1,3])
      s <- s+h
      idx <- idx+1
    } else {
      id <- id+1
    }
  }
  return(int_out)
}

out <- pdmp(1000)
int_out <- interpolate(out,1000,1000)

plot <- ggplot()+
  xlim(mu-(sqrt(sigma)*5),mu+(sqrt(sigma)*5))+
  geom_function(fun=f)+
  geom_density(aes(x=int_out[,2]),alpha=0.5,col='red')

print(plot)