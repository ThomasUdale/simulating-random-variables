library(ggplot2)
library(progress)


f <- function(x) {
  exp(-x^4/2)/2.1558
}

f <- function(x) {
  dnorm(x,0,sqrt(0.5))
}

mu_0 = 6
sub_means = c(0.5*mu_0,-0.5*mu_0)
C = length(sub_means)

fc <- function(x,c) {
  dnorm(x,sub_means[c],1)
}

samplefc <- function(c) {
  rnorm(1,sub_means[c],1)
}

mu <- function(x) {
  1/C*sum(fc(x,c(1:C)))
}

sample_mu <- function() {
  samplefc(sample(c(1:C),1))
}

sample_p_lambda <- function(x) {
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

jump_sample <- function(total_time){
  samples <- data.frame()
  x = sample_mu()
  i = 0
  t = 0
  
  pb <- progress_bar$new(
    format = "sampling [:bar] :current/:total (:percent) eta: :eta, elapsed: :elapsed",
    clear = FALSE, total = total_time, width = 80)
  
  while(t<total_time){
    pb$update(t/total_time)
    i <- i + 1
    f <- f(x)
    tau_step <- rexp(1,lambda(x)/f)
    tau_regen <- rexp(1,mu(x)/f)
    tau <- min(tau_step,tau_regen)
    samples <- rbind(samples,c(i,tau,x,tau_step<tau_regen))
    t <- t + tau
    if (tau_step<tau_regen){
      x <- sample_p_lambda(x)
    } else {
      x <- sample_mu()
    }
    
    colnames(samples) <- c('index','sample_time','x','step')
  }
  pb$terminate()
  return(samples)
}

output <- jump_sample(10000)
plot <- ggplot(output, aes(x)) + geom_histogram(aes(y=..density..,weight=sample_time),binwidth=0.1)+geom_function(fun = f, colour = "red")+xlim(-4, 4)+geom_density(aes(weight=sample_time),color='blue')
print(plot)
