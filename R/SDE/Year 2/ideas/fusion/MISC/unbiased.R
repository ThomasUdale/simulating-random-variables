library(ggplot2)

f <- function(x) {
  dnorm(x,4,sqrt(2/3))
}


mu = c(10,1)
sigma = c(sqrt(2),1)

fc <- function(x,c) {
  dnorm(x,mu[c],sigma[c])
}

samplefc <- function(c) {
  rnorm(1,mu[c],sigma[c])
}

s1 <- rnorm(1e4,mu[1],sigma[1])
s2 <- rnorm(1e4,mu[2],sigma[2])

d1 <- dnorm(s1,mu[1],sigma[1])
d2 <- dnorm(s2,mu[2],sigma[2])

te <- function() {
  N <- rgeom(1,2/5)
  x <- sapply(c(1:2),samplefc)
  cur_d <- sd(x)
  cur_f <- mean(x)
  z <- cur_f
  if (N==0){
    return(z)
  } else {
    for (i in 1:N) {
      new_x <- x
      for (j in 1:100){
        id <- sample(2,1)
        new_x[id] <- samplefc(id)
        if (sd(new_x)<cur_d) {
          x <- new_x
          cur_f <- mean(x)
          print(cur_f)
          
        }
      }
      z <- z + (mean(x)-cur_f)*exp(i)
      cur_d <- sd(x)
    }
    return(z)
  }
}

plot <- ggplot()+
  xlim(-5,15)+
  geom_function(fun=f)+
  geom_function(fun=fc,args=list(c=1),col='blue',alpha=0.5)+
  geom_function(fun=fc,args=list(c=2),col='blue',alpha=0.5)+
  geom_vline(xintercept=4,alpha=0.5)

print(plot)