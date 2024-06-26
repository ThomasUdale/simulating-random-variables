mu = 0
sd = 1


density <- function(x,noise=0) {
  dnorm(x,mu,sd)*(rnorm(1,1,noise))
  # # dnorm(x)*rexp(1,0.1+10*x*x)
  # dnorm(x)*rgamma(1,0.1+10*x*x,0.1+10*x*x)
}

op=par(mfrow=c(1,2))
curve(dnorm,-5,5,ylim=c(0,0.5))

noise = 0.1

x = 0
w=density(x,noise)
sample = c()
burn_in = 10000
target_ss = 10000
n = 0
accept_count = 0


while (n<burn_in+target_ss) {
  n=n+1
  
  new_x = rnorm(1,x)
  new_w = density(new_x,noise)
  
  r = new_w/w
  a = min(1,r)
  
  if (runif(1) < a) {
    x = new_x
    w = new_w
    accept_count = accept_count + 1 
  }
  
  sample = append(sample,x)
}


hist(sample[(burn_in+1):(burn_in+target_ss)],freq=F,add=T)
acf(sample)