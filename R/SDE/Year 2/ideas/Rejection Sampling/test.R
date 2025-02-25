library(ggplot2)


rej <- function(){
  x <- c()
  w <- c()
  while(TRUE){
    new_x <- rnorm(1,0,sqrt(2))
    new_w <- dnorm(new_x)/dnorm(new_x,0,sqrt(2))
    if (new_w<1 & length(x)>100){
      if (runif(1)<new_w) {
        if (length(w)==0){
          return(new_x)
        } else {
          x <- append(x,new_x)
          w <- append(w,new_w)
          return(sample(x,1,prob=w))
        }
      }
    } else {
      x <- append(x,new_x)
      w <- append(w,new_w)
    }
  }
}

out <- pbapply::pbreplicate(1e4,rej())

plot <- ggplot()+
  geom_density(aes(x=out))+
  geom_function(fun=dnorm,col='red')+
  geom_function(fun=dnorm,args=list(sd=sqrt(2)),col='blue')


print(plot)