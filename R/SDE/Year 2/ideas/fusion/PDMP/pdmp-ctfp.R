library(ggplot2)

mu <- 0
sigma <- 1

f <- function(x) {
  dnorm(x,mu,sqrt(sigma))
}

grad_log_f <- function(x) {
  -(x-mu)/sigma
}

mu_time <- function(x,v){
  cur_x <- x
  total <- 0
  delta <-1
  while(TRUE) {
    if (v*cur_x>0){
      s <- rexp(1,(exp((cur_x+v*delta)^2)/4-1))
      if (s>delta){
        total <- total+delta
        cur_x <- cur_x+v*delta
      } else {
        if (runif(1)<(exp((cur_x+v*s)^2/4)-1)/(exp((cur_x+v)^2)/4-1)) {
          return(s+total)
        } else {
          total <- total+s
          cur_x <- cur_x+v*s
        }
      }
    }
  }
}

mu_time_0 <- function(x,v) {
  if (v*x>0){
    total <- 0
    while(TRUE){
      s <- rexp(1,exp(-x^2)/2)
      if (runif(1)<exp(-(x+v*s)^2/2)/exp(-x^2)/2) {
        return(total+s)
      } else {
        total <- total+s
      }
    }
    
  } else {
    total <- 0
    cur_x <- x
    while(TRUE){
      s <- rexp(1,exp(-(cur_x+v*1)^2/2))
      if (s>1){
        total <- total+1
        cur_x <- x+v*1
      } else {
        if (runif(1)<exp(-(cur_x)^2/2)/exp(-(cur_x+v*s)^2)/2) {
          return(s+total)
        } else {
          total <- total+s
          cur_x <- cur_x+v*s
        }
      }
      
    }
  }
}

pdmp <- function(t) {
  v <- sample(c(-1,1),1)
  x <- rnorm(1,mu,sqrt(sigma))
  s <- 0
  out <- c(0,x,v,0)
  while (s<t) {
    e <- rexp(1)
    # jt <- rexp(1,10/t)
    jt <- mu_time_0(x,v)
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
      x <- rnorm(1,0,1/sqrt(2))
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