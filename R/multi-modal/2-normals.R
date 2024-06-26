# Simulate from mixture of two normals:
p = 0.3
mu1 = 1
sigma1 = 1
mu2 = 20
sigma2 = 4


density <- function(x){
  p*dnorm(x,mu1,sigma1) + (1-p)*dnorm(x,mu2,sigma2)
}

curve(density,-5,30,ylim=c(0,0.5))


# Metropolis-Hastings
x = 0
sample = c(x)
burn_in = 10000
target_ss = 10000

proposals = 0

while (length(sample)<burn_in+target_ss){
  proposals = proposals + 1
  candidate = rnorm(1,x)
  acceptance = min(1,(density(candidate))/(density(x)))
  u = runif(1)
  if (u<acceptance){
    x <- candidate
  }
  sample <- append(sample,x)
}


hist(sample[burn_in+1:length(sample)],add=T,freq=F,nclass=100)

# Mixing Disaster!


