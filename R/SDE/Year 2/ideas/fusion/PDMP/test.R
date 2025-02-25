library(ggplot2)


rejection <- function(){
  while(TRUE){
    x <- rnorm(1)
    if (runif(1)<exp(-0.5*x^2)) {
      return(x)
    }
  }
}

out <- replicate(1e4,rejection())

plot <- ggplot()+
  xlim(-3,3)+
  geom_function(fun=dnorm,args=list(sd=sqrt(0.5)),col='green')+
  geom_function(fun=dnorm)+
  geom_density(aes(x=out),col='red')

print(plot)