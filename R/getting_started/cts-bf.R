library(MASS)
library(ggplot2)
set.seed(4)

tp <- seq(0,1,0.01)
n <- length(tp)
C <- 2
N <- 2

p <- 0.1

g <- rep(1,N)

x <- array(dim=c(n,C,N))
table_x <- data.frame()
swaps <- data.frame()

add_rows <- function(table_x,times,particles=1:N,j="0",g){
  for (t in times){
    for (p in particles) {
      add_table <- reshape2::melt(x[t,,p])
      add_table$t <- t
      add_table$p <- p
      add_table$j <- j
      add_table$c <- 1:C
      add_table$g <- g[p]
      table_x <- rbind(table_x,add_table)
    }
  }
  return(table_x)
} 

x[1,,] <- t(mvrnorm(n=N,mu=-((C+1)/2)+(1:C),Sigma=diag(0.01,nrow=C)))
table_x <- add_rows(table_x,c(1),j="0",g=g)

for (t in 2:n) {
  for (i in 1:N) {
    
    mean <- (1-tp[t])/(1-tp[t-1])*x[t-1,,i] + (tp[t]-tp[t-1])/(1-tp[t-1]) * mean(x[t-1,,i])
    cov <- diag((tp[t]-tp[t-1])*(1-tp[t])/(1-tp[t-1]),nrow=C) + matrix((tp[t]-tp[t-1])^2/(C*(1-tp[t-1])),nrow = C,ncol=C)
    x[t,,i] <- mvrnorm(n=1,mu=mean,Sigma=cov)
    table_x <- add_rows(table_x,c(t),c(i),j="0",g)
  }
  
  if (runif(1) <p && t<n){
    from_p <- sample(1:N,1)
    to_p <- sample(1:N,1)
    g[from_p] <- g[from_p]+1
    for (c in 1:C){
      swaps <- rbind(swaps,c(t,from_p,x[t,c,from_p],x[t,c,to_p],c,g[from_p]))
    }
    x[t,,from_p] = x[t,,to_p]
    table_x <- add_rows(table_x,c(t),c(from_p),j="1",g)
    
  }
}

colnames(swaps) <- c("t","p","x_from","x_to","c","g")


plot <- ggplot()+geom_line(data=table_x,aes(x=t,y=value,col=interaction(p,c),group=interaction(g,p,c)))+geom_segment(data=swaps,aes(x=t,y=x_from,yend=x_to,color=interaction(p,c)),linetype = 4)
print(plot)

