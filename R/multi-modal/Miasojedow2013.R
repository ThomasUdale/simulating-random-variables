# Simulate from mixture of two normals:
p = 0.3
mu1 = 1
sigma1 = 1
mu2 = 20
sigma2 = 4

density <- function(x,inv_temp=1){
  (p*dnorm(x,mu1,sigma1) + (1-p)*dnorm(x,mu2,sigma2))^inv_temp
}

curve(density,-30,30,ylim=c(0,0.2))


# Adaptive Parallel Tempering Metropoolis-Hastings

burnin = 10000
num_samples = 10000
n=0
L = 5
alpha_star = 0.234

rho = rep(1,L-1)
beta = 0.1^(0:(L-1))
x = rep(0,L)
samples = c()

s_count = 0


pb = txtProgressBar(min = 0, max = (burnin+num_samples), initial = 0,style=3) 
while (n+1<burnin+num_samples) {
  n = n+1
  
  setTxtProgressBar(pb,n)
  
  
  # Swap
  swap_ps = c()
  for (i in 1:(L-1)){
    swap_p = min(1,(density(x[i],beta[i+1])*density(x[i+1],beta[i])/(density(x[i],beta[i])*density(x[i+1],beta[i+1]))))
    swap_ps = append(swap_ps,swap_p)
  }
  
  swap_candidate = sample(1:(L-1),1)
  if (runif(1)<swap_ps[swap_candidate]){
    x[swap_candidate:swap_candidate+1] = rev(x[swap_candidate:swap_candidate+1])
    if(swap_candidate==1){s_count = s_count +1}
  }
  
  # MH
  new_x = c()
  for (i in 1:L){
    candidate = rnorm(1,x[i])
    if (runif(1)<min(1,(density(candidate,inv_temp=beta[i]))/(density(x[i],inv_temp=beta[i])))) {
      x[i] = candidate
    }
  }
  
  samples = append(samples,x[1])
  
  
  
  # Update temps
  gamma = (n+1)^(-0.6)
  
  T = 1
  for (l in 1:L) {
    if (l<L){
      rho[l] = rho[l]+gamma*(swap_ps[l]-alpha_star)
    }
    if (l>1){
      T = T + exp(rho[l-1])
      beta[l] = 1/T
    }
  }
}

close(pb)

sample = samples[(burnin+1):length(samples)]

hist(sample,add=T,freq=F,nclass=50)