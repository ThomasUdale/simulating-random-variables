# Simulate from mixture of two normals:
sample_size = 1000000
p = 0.3
mu1 = 1
sigma1 = 1
mu2 = 20
sigma2 = 4
k <- rbinom(1,sample_size,p)
sample = c(rnorm(k,mu1,sigma1),rnorm(sample_size-k,mu2,sigma2))
hist(sample,freq = F)

density <- function(x,temp=1){
  (p*dnorm(x,mu1,sigma1) + (1-p)*dnorm(x,mu2,sigma2))^temp
}

curve(density,-5,30,ylim=c(0,0.5))


# Metropoolis-Hastings

sample = c(0)
burn_in = 10000
target_ss = 10000
temps = 0.5^(0:3)
nchains = length(temps)
x = rep(0,nchains)

proposals = 0

while (length(sample)<burn_in+target_ss){
  proposals = proposals + 1
  
  swaps = sample(1:nchains,2)
  accept_swap = runif(1) < (density(x[swaps[1]],temps[swaps[2]])*density(x[swaps[2]],temps[swaps[1]]))/(density(x[swaps[1]],temps[swaps[1]])*density(x[swaps[2]],temps[swaps[2]]))
  if (accept_swap){x[swaps] = rev(x[swaps])}
  
  
  for (i in 1:nchains){
    candidate = rnorm(1,x[i])
    acceptance = min(1,(density(candidate,temp=temps[i]))/(density(x[i],temp=temps[i])))
    u = runif(1)
    if (u<acceptance){
      x[i] <- candidate
      if (i == 1){sample <- append(sample,candidate)}
    }
  }
  
}


hist(sample[burn_in+1:length(sample)],add=T,freq=F,nclass=100)
