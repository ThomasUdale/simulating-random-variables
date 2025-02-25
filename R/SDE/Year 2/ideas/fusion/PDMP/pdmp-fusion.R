library(ggplot2)
library(mvtnorm)

f <- function(x) {
  dnorm(x)
}


grad_log_f <- function(x) {
  -x
}

target_time <- 2

pdmp <- function(n) {
  mu_events <- cumsum(rexp(n,0.1))
  exit_time <- mu_events[n]
  v <- sample(c(-1,1),1)
  x <- rnorm(1)
  s <- 0
  next_event <- 1
  w <- 0
  out <- c(s,x,v,w,2)
  while (s<exit_time) {
    e <- rexp(1)
    if (v*x<0){
      event_time <- abs(x) + sqrt(2*e)
    } else {
      event_time <- -v*x + sqrt(x^2 + 2*e)
    }
    
    if (s+event_time>mu_events[next_event]){
      new_x <- (target_time*next_event/n + x*(1-next_event/n))
      w <- w+log(f(new_x))-log(dnorm(new_x,x,1))
      s <- mu_events[next_event]
      next_event <- next_event+1
      x <- new_x
      v <- sample(c(-1,1),1)
      out <- rbind(out,c(s,x,v,w,1))
    } else {
      x <- x+event_time*v
      v <- -v
      s <- s+event_time
      out <- rbind(out,c(s,x,v,w,0))
    }
    
  }
  return(out)
}

# out <- pdmp(25)
# print(out[nrow(out),4])
# 
# plot <- ggplot()+
#   geom_point(aes(x=out[,1],y=out[,2],color=as.factor(out[,5])))+
#   geom_line(aes(x=out[,1],y=out[,2]),alpha=0.1)
# 
# print(plot)

mu = c(10,1)
sigma = c(sqrt(2),1)
C=2



pf_single <- function(id,n,paths=FALSE) {
  mu_ts <- cumsum(rexp(n,0.1))
  exit_time <- mu_ts[n]
  v <- sample(c(-1,1),C,replace=TRUE)
  x <- rnorm(2,mu,sigma)
  s <- rep(0,C)
  w <- 0
  out <- c(s,x,v,w,2)
  out <- cbind(s,x,v,w,2,c(1:C))
  for (i in 1:n) {
    for( j in 1:C) {
      while(s[j]<mu_ts[i]) {
        e <- rexp(1)
        if (v[j]*(x[j]-mu[j])<0){
          et <- abs(x[j]-mu[j]) + sqrt(2*e*sigma[j]^2)
        } else {
          et <- -v[j]*(x[j]-mu[j]) + sqrt((x[j]-mu[j])^2 + 2*e*sigma[j]^2)
        }
        
        if (s[j]+et>mu_ts[i]) {
          x[j] <- x[j]+v[j]*(mu_ts[i]-s[j])
          s[j] <- mu_ts[i]
        } else {
          x[j] <- x[j]+et*v[j]
          v[j] <- -v[j]
          s[j] <- s[j]+et
          out <- rbind(out,c(s[j],x[j],v[j],w,0,j))
        }
      }
      
    }
    new_x <- rnorm(C,x*(1-mu_ts[i]/exit_time) + mean(x)*(mu_ts[i]/exit_time))
    new_x <- rmvnorm(
      1,
      x*(1-mu_ts[i]/exit_time) + mean(x)*(mu_ts[i]/exit_time),
      diag(rep(mu_ts[i]*(exit_time-mu_ts[i])/exit_time,C))+
        matrix(mu_ts[i]^2/(2*exit_time^2),nrow=C,ncol=C)
    )
    w <-w+
      sum(log(dnorm(new_x,mu,sigma)))-
      dmvnorm(
        new_x,
        x*(1-mu_ts[i]/exit_time) + mean(x)*(mu_ts[i]/exit_time),
        diag(rep(mu_ts[i]*(exit_time-mu_ts[i])/exit_time,C))+
          matrix(mu_ts[i]^2/(2*exit_time^2),nrow=C,ncol=C),
        log=TRUE
      )
    x <- new_x[1,]
    v <- rep(sample(c(-1,1),1),C)
    out <- rbind(out,cbind(s,x,v,w,1,c(1:C)))
  }
  if (!paths) {
    return(c(out[nrow(out),2],out[nrow(out),4]))
  } else {
    return(out)
  }
  
}

# sample <- do.call(rbind,
#   pbmcapply::pbmclapply(c(1:1e4),pf_single,n=15,mc.cores=10))
# 
# print(summary(exp(sample[,2])))
# print(summary(sample[,1]))
# 
# plot <- ggplot()+
#   xlim(min(sample[,1]),max(sample[,1]))+
#   geom_function(fun=dnorm,args=list(mean=4,sd=sqrt(2/3)))+
#   geom_density(aes(x=sample[,1],weight=exp(sample[,2])),alpha=0.4,col='blue')
# 
# print(plot)

out <- pf_single(1,20,paths=T)

plot <- ggplot()+
  geom_point(aes(x=out[,1],y=out[,2],group=out[,6],color=as.factor(out[,5])))+
  geom_line(aes(x=out[,1],y=out[,2],group=out[,6]),alpha=0.1)

print(plot)