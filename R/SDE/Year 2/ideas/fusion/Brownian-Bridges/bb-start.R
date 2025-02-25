library(ggplot2)
library(progress)
library(data.table)
library(mvtnorm)

f <- function(x) {
  dnorm(x,4,sqrt(2/3))
}

C = 2
mu = c(10,1)
sigma = c(sqrt(2),1)

proposal <- function(C,t,n) {
  x <- c(10,1)
  weight <- 0
  s <- 0
  out <- matrix(nrow=n+1,ncol=C+2)
  out[1,] <- c(s,x,weight)
  for (i in 1:n) {
    new_x <- rmvnorm(
      1,
      x*(1-i/n) + mean(x)*(i/n),
      diag(rep(t*i*(n-i)/n,C))+
        matrix(t*i^2/(C*n),nrow=C,ncol=C)
    )
    x_den <- dmvnorm(
      new_x,
      x*(1-i/n) + mean(x)*(i/n),
      diag(rep(t*i*(n-i)/n,C))+
        matrix(t*i^2/(C*n),nrow=C,ncol=C),
      log=TRUE
    )
    x <- new_x
    s <- s+t/n
    if (i<n){
      weight <- weight+sum(dnorm(x,mu,sigma,log=TRUE))-x_den
    } else {
      weight <- weight+sum(dnorm(x,mu,sigma,log=TRUE))
    }
    
    out[(i+1),] <- c(s,x,weight)
  }
  return(out)
}



out <- melt(
  data.table(proposal(C,1,1000)),
  id.vars=c(1,C+2),
  # measure.vars=list(c(2:(C+1))),
  value.name=c("value"),
)

plot <- ggplot()+
  geom_line(data=out,aes(x=V1,y=value,group=variable))

print(plot)

prop_2 <- function(C,t,n,np){
  x <- t(replicate(np,rnorm(C,mu,sigma)))
  weight <- rep(log(1/np),np)
  s <- 0
  pb <- progress_bar$new(
    format = "sampling [:bar] :current/:total (:percent) eta: :eta, elapsed: :elapsed",
    clear = FALSE, total = n*np, width = 80)
  for (i in 1:n){
    cov <- diag(rep(t*i*(n-i)/n,C))+
      matrix(t*i^2/(C*n),nrow=C,ncol=C)
    for (j in 1:np) {
      pb$tick()
      new_x <- rmvnorm(
        1,
        x[j,]*(1-i/n) + mean(x[j,])*(i/n),
        cov
      )
      x_den <- dmvnorm(
        new_x,
        x[j,]*(1-i/n) + mean(x[j,])*(i/n),
        cov,
        log=TRUE
      )
      
      
      if (i<n){
        weight[j] <- weight[j]+
          sum(dnorm(new_x,mu,sigma,log=TRUE))-
          sum(dnorm(x[j,],mu,sigma,log=TRUE))+
          dmvnorm(new_x,x[j,],log=TRUE)-
          x_den
      } else {
        weight[j] <- weight[j]+
          sum(dnorm(new_x,mu,sigma,log=TRUE))-
          sum(dnorm(x[j,],mu,sigma,log=TRUE))+
          dmvnorm(new_x,x[j,],log=TRUE)
      }
      x[j,] <- new_x
    }
    # if (i<n){
    #   print(1/sum(exp(weight)^2))
    #   new_id <- sample(np,np,replace=T,prob=exp(weight))
    #   x <- x[new_id,]
    #   weight <- rep(log(1/np),np)
    #   print(length(unique(new_id)))
    # }
    s <- s+t/n
  }
  pb$terminate()
  return(cbind(x[,1],weight))
}
# 
# for(i in c(1e4)){
#   out2 <- prop_2(2,10,10,i)
#   
#   plot <- ggplot()+
#     xlim(0,8)+
#     geom_function(fun=f)+
#     geom_density(aes(x=out2[,1],weight=exp(out2[,2])),col='blue',alpha=0.5)
#   
#   print(plot)
# }
