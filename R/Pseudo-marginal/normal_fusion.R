
num_dist = 2
num_is = 100
mu = 1:num_dist
sd = (num_dist:1)*2

prod_of_norms <- function(x) {
  prod(dnorm(x,mu,sd))
}



h <- Vectorize(prod_of_norms)

q = integrate(h,lower=-Inf,upper=Inf)$value

normh <- function(x) {
  h(x)/q
}


curve(normh,-10,30,ylim=c(0,0.5))

y = 0
r = 1e-6
samples = c()
n = 0
accepts = 0
burn = 10000
target= 10000

pb = txtProgressBar(min = 0, max = (burn+target), initial = 0,style=3) 

while (n<(burn+target)) {
  n = n+1
  setTxtProgressBar(pb,n)
  
  ynew = rnorm(1,y)
  
  new_r = 0
  for (i in 1:num_is) {
    r_star = 1
    for (c in 1:num_dist) {
      x = rnorm(1,mu[c],sd[c])
      r_star = r_star*(dnorm(ynew,mu[c],sd[c])/dnorm(x,mu[c],sd[c]))*dnorm(ynew,x)
    }
    new_r = new_r + r_star/num_is
  }
  
  a = new_r/r
  
  if (runif(1)<a){
    accepts = accepts+1
    y = ynew
    r = new_r
  }
  
  if (n>burn) {samples = append(samples,y)}
  
  
}

close(pb)
hist(samples,freq=F,add=T,nclass=50)