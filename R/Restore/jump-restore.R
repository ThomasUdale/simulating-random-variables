C = 4
library(ggplot2)
library(progress)


f <- function(x) {
  exp(-x^4/2)/2.1558
}

samplef <- function(){
  while (TRUE) {
    x <- rnorm(1)
    if (runif(1)<(fc(x)/exp(-x^2/2+1/8))) {
      return(x)
    }
  }
}


fc <- function(x) {
  exp(-x^4/(2*C))
}

samplefc <- function(c) {
  while (TRUE) {
    x <- rnorm(1)
    if (runif(1)<(fc(x)/exp(-x^2/2+C/8))) {
      return(x)
    }
  }
}

mu <- function(x){
  dnorm(x,4,4)
}

sample_mu <- function() {
  rnorm(1,4,4)
}

sample_p <- function(x) {
  y <- rnorm(1,x)
  a <- dnorm(y)/dnorm(x)
  if (runif(1)<a){
    return(y)
  } else {
    return(x)
  }
}

lambda <- function(x) {
  dnorm(x)
}

k_lower <- 0.1

jump_sample <- function(index){
  x = sample_mu()
  i = 0
  while(TRUE){
    i <- i + 1
    p_step <- lambda(x)
    p_regen <- mu(x)-k_lower*f(x)
    p_exit <- k_lower*f(x)
    
    tau <- sample(c(1,2,3),1,prob=c(p_step,p_regen,p_exit))
    
    if (tau==1){
      x <- sample_p(x)
    } else if (tau==2){
      x <- sample_mu()
    } else {
        return(x)
    }
  }
}

output <- jump_sample()
print(output)
# # plot <- ggplot(output, aes(x)) + geom_histogram(aes(y=..density..,weight=state_time),binwidth=0.1)+geom_function(fun = f, colour = "red")+xlim(-4, 4)
# # print(plot)
# 
samples <- pbmcapply::pbmclapply(seq(100000),jump_sample,mc.cores=12)

hist(unlist(samples),freq=F,nclass=50,xlim=c(-3.5,3.5),ylim=c(0,0.6))
curve(f,add=T,col='red')
lines(density(unlist(samples)),col='blue')

