library(ggplot2)

ir <- function(t) {
  x <- rnorm(1,0,sqrt(0.5))
  s <- rexp(1,dnorm(x,0,sqrt(0.5))/dnorm(x))
  out <- c(x,min(s,t))
  while(s<t) {
    x <- rnorm(1,0,sqrt(0.5))
    new_s <- rexp(1,dnorm(x,0,sqrt(0.5))/dnorm(x))
    out <- rbind(out,c(x,min(new_s,t-s)))
    s <- s+new_s
  }
  return(out)
}

ir_single <- function(id){
  x <- rnorm(1,2,sqrt(4))
  c(x,rexp(1,dnorm(x,2,sqrt(4))/dnorm(x)))
}

out <- do.call(rbind,pbmcapply::pbmclapply(c(1:1e6),ir_single))

plot <- ggplot()+
  xlim(-4,4)+
  geom_function(fun=dnorm,args=list(mean=2,sd=sqrt(4)))+
  geom_function(fun=dnorm,col='red')+
  geom_density(aes(x=out[,1],weight=out[,2]),col='blue')
  

print(plot)