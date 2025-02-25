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

d <- d1*d2

means <- rowMeans(cbind(s1,s2))

split_mag <- exp(-abs(s1-s2))

plot <- ggplot()+
  geom_function(fun=f)+
  geom_density(aes(x=means),alpha=0.5,col='red')+
  geom_density(aes(x=means,weight=split_mag),alpha=0.5,col='green')+
  # geom_histogram(aes(x=means,y=after_stat(density)),alpha=0.5)+
  # geom_histogram(aes(x=means,weight=split_mag,y=after_stat(density)),col='red',alpha=0.5)+
  geom_vline(xintercept=4,alpha=0.5)+
  geom_vline(xintercept=5.5,col="red",alpha=0.5)

print(plot)