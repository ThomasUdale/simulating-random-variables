library(ggplot2)

f <- function(x){
  exp(-x^4/2)
}


h <- 0.0001
target <- 1
z <- 0
time <- 0
integral <- 0

while(time<target){
  time <- time + h
  z <- z + f(z)*h
  integral <- integral + h*2/(f(z)^2)
}

est_q <- function(h,x){
  z <- 0
  time <- 0
  while(time<abs(x)){
    if (x>0){
      z <- z + f(z)*h
    } else {
      z <- z-f(z)*h
    }
    time <- time+h
  }
  z
}



sim <- function(){
  B_g <- rnorm(1,0,integral^2)
  X <- est_q(h,B_g)
  return(B_g)
}

simulation <- pbapply::pbreplicate(1e3,sim())
hist(simulation,freq=F)

# plot <- ggplot()+
#   geom_line(aes(x=times,y=zs))+
#   geom_line(aes(x=times,y=ps))
# print(plot)