# Simulate from mixture of two normals:

# devtools::load_all("~/Documents/simulating-random-variables/R-packages/apt.rust")
use("apt.rust")

p = 0.3
mu1 = 1
sigma1 = 1
mu2 = 20
sigma2 = 4

density <- function(x,inv_temp=1){
  (p*dnorm(x,mu1,sigma1) + (1-p)*dnorm(x,mu2,sigma2))^inv_temp
}

curve(density,-5,30,ylim=c(0,0.5))



# Adaptive Parallel Tempering Metropoolis-Hastings




burnin = 1e6
num_samples = 1e6

# rust implementation
sample2 = apt(burnin,num_samples,FALSE)
hist(sample2[1:1e6],add=T,freq=F,nclass=200,col='red')

burnin = 100000
num_samples = 100000

n=0
L = 5
alpha_star = 0.234

rho = rep(1,L-1)
beta = 0.5^(0:(L-1))
x = rep(0,L)
samples = c()

s_count = 0

var = rep(1,L)
Rho = rep(1,L)
mu_var = 0

T_var = rep(1,L)

pb = txtProgressBar(min = 0, max = (burnin+num_samples), initial = 0,style=3) 
while (n<burnin+num_samples) {
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
    x[swap_candidate:(swap_candidate+1)] = rev(x[swap_candidate:(swap_candidate+1)])
    if(swap_candidate==1){s_count = s_count +1}
  }
  
  # MH
  T_gamma_var = n^(-0.6)
  gamma_var = 0.5*n^(-0.6)
  
  for (i in 1:L){
    candidate = rnorm(1,x[i],var[i])
    
    # Update variance
    
    mu_var = (1-gamma_var)*mu_var + gamma_var*mean(x)
    Rho[i] = (1-gamma_var)*Rho[i] + gamma_var*1/L*sum((x-mu_var)^2)
    
    acceptance_p = min(1,(density(candidate,inv_temp=beta[i]))/(density(x[i],inv_temp=beta[i])))
    
    
    T_var[i] = T_var[i] + T_gamma_var*(acceptance_p-alpha_star)
    var[i] = exp(T_var[i])*Rho[i]
    
    
    
    if (runif(1)<acceptance_p) {
      x[i] = candidate
    }
  }
  
  samples = append(samples,x[1])
  
  
  
  # Update temps
  gamma = (n)^(-0.6)

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
hist(sample,add=T,freq=F,nclass=200)
