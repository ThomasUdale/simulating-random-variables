library(ggplot2)

f <- function(x) {
  dnorm(x,4,sqrt(2/3))
}

log_f <- function(x) {
  log(dnorm(x,4,sqrt(2/3)))
}

mu = c(10,1)
sigma = c(sqrt(2),1)

fc <- function(x,c) {
  dnorm(x,mu[c],sigma[c])
}

samplefc <- function(c) {
  rnorm(1,mu[c],sigma[c])
}

grad_log_fc <- function(x,c) {
  -(x-mu[c])/sigma[c]^2
}

iad <- function(samples,weights){
  d <- density(samples,weights = weights/sum(weights),from=-10,to=10,n=2000)
  sum(abs(d$y - f(d$x)))/2000*20
}


N <- 1e5

mean_is <- function(N) {
  s1 <- rnorm(N,mu[1],sigma[1])
  s2 <- rnorm(N,mu[2],sigma[2])
  
  grad_weights1 <- exp(sapply(s1,grad_log_fc,c=1))
  grad_weights2 <- exp(sapply(s2,grad_log_fc,c=2))
  
  
  
  
  means <- (s1*grad_weights1+s2*grad_weights2)/(grad_weights1+grad_weights2)
  means <- rowMeans(cbind(s1,s2))
  
  s <- rnorm(N,means)
  
  weights <- (dnorm(s,mu[1],sigma[1])*dnorm(s,mu[2],sigma[2]))/dnorm(s,means)
  
  return(iad(s,weights))
}

sum_grad <- function(x) {grad_log_fc(x,1)+grad_log_fc(x,2)}

plot <- ggplot()+
  xlim(-5,15)+
  ylim(-5,1)+
  geom_function(fun=f)+
  geom_function(fun=fc,args=list(c=1),col='blue')+
  geom_function(fun=fc,args=list(c=2),col='blue')+
  geom_function(fun=grad_log_fc,args=list(c=1),col='red',alpha=0.1)+
  geom_function(fun=grad_log_fc,args=list(c=2),col='red',alpha=0.1)+
  geom_function(fun=sum_grad)+
  geom_function(fun=log_f)

print(plot)

# plot <- ggplot()+
#   xlim(-3,15)+
#   ylim(0,0.5)+
#   geom_density(aes(x=s1),col='red',alpha=0.7)+
#   geom_density(aes(x=s2),col='green',alpha=0.7)+
#   geom_density(aes(x=s,weight=weights),col='orange',alpha=0.2)+
#   geom_function(fun=f,col='blue')
# 
# print(plot)





################################ EXP 
# 
# C=4
# 
# f <- function(x) {
#   exp(-x^4/2)/2.155
# }
# 
# fc <- function(x) {
#   exp(-x^4/(2*C))
# }
# 
# fc_norm <- function(x) {
#   exp(-x^4/(2*4))/3.04876
# }
# 
# samplefc <- function() {
#   while (TRUE) {
#     x <- rnorm(1)
#     if (runif(1)<(fc(x)/exp(-x^2/2+C/8))) {
#       return(x)
#     }
#   }
# }
# 
# sample_1 <- function() {
#   x <- replicate(C,samplefc())
#   y <- rnorm(1,mean(x))
#   w <- fc(y)^C/dnorm(y,mean(x))
#   return(c(y,w))
# }
# 
# sample <- replicate(1e5,sample_1())
# 
# plot <- ggplot()+
#   xlim(-5,5)+
#   ylim(0,1)+
#   geom_density(aes(x=sample[1,],weight=sample[2,]),col='orange',alpha=0.5)+
#   geom_function(fun=fc_norm,col='blue',alpha=0.2)+
#   geom_function(fun=f)
# 
# print(plot)