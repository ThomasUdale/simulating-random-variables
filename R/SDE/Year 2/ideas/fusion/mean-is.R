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

N <- 1e5

s1 <- rnorm(N,mu[1],sigma[1])
s2 <- rnorm(N,mu[2],sigma[2])

means <- rowMeans(cbind(s1,s2))

s <- rnorm(N,means)

weights <- (dnorm(s,mu[1],sigma[1])*dnorm(s,mu[2],sigma[2]))/dnorm(s,means)

weighted_means <- sample(s,N,replace=T,prob=weights)

plot <- ggplot()+
  xlim(-3,15)+
  geom_density(aes(x=s1),col='red',alpha=0.7)+
  geom_density(aes(x=s2),col='green',alpha=0.7)+
  geom_density(aes(x=s,weight=weights),col='orange',alpha=0.2)+
  geom_density(aes(x=s),alpha=0.2)+
  geom_function(fun=f,col='blue')

print(plot)


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