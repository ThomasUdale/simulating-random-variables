library(ggplot2)

f <- function(x){
  dnorm(x)
}

log_f <- function(x){
  -0.5*x^2
}

log_rej <- function(){
  while(T){
    y <- rnorm(1,2,sqrt(2))
    if (rexp(1)>(-(y-2)^2/4-log_f(y))){
      return(y)
    }
  }
}

s <- replicate(1e4,log_rej())

plot <- ggplot()+
  xlim(-3,3)+
  geom_function(fun=f)+
  geom_density(aes(x=s))

print(plot)


