# Simulate from mixture of two normals:
p = 0.3
mu1 = 1
sigma1 = 1
mu2 = 20
sigma2 = 4

density <- function(x,temp=1){
  (p*dnorm(x,mu1,sigma1) + (1-p)*dnorm(x,mu2,sigma2))^temp
}

curve(density,-5,30,ylim=c(0,0.5))


# Metropoolis-Hastings

sample = c()
burn_in = 20000
target_ss = 20000
temps = 0.5^(0:5)
nchains = length(temps)
x = rep(0,nchains)

proposals = 0
scount = 0

pb = txtProgressBar(min = 0, max = (burnin+num_samples), initial = 0,style=3) 
while (length(sample)<burn_in+target_ss){
  proposals = proposals + 1
  setTxtProgressBar(pb,proposals)
  
  swaps = sample(1:(nchains-1),1)
  accept_swap = runif(1) < (density(x[swaps],temp=temps[swaps+1])*density(x[swaps+1],temp=temps[swaps]))/(density(x[swaps],temp=temps[swaps])*density(x[swaps+1],temp=temps[swaps+1]))
  if (accept_swap){
    x[swaps:(swaps+1)] = rev(x[swaps:(swaps+1)])
    scount = scount + 1
  }
  
  
  for (i in 1:nchains){
    candidate = rnorm(1,x[i])
    if (runif(1)<min(1,(density(candidate,temp=temps[i]))/(density(x[i],temp=temps[i])))){
      x[i] <- candidate
    }
    if (i == 1){sample <- append(sample,x[i])}
  }
  
}

close(pb)
hist(sample[(burn_in+1):length(sample)],add=T,freq=F,nclass=200)

