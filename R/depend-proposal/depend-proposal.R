p = 0.3
mu1 = 1
sigma1 = 1
mu2 = 20
sigma2 = 4


density <- function(x){
  p*dnorm(x,mu1,sigma1) + (1-p)*dnorm(x,mu2,sigma2)
}

curve(density,-5,30,ylim=c(0,0.5))

step_p <- function(x,new_x) {
  
  if ((new_x-x)>0) {
    if ((new_x-x)<1) {
      u = density(new_x)
      l = density(2*x-new_x)
      
      u/(u+l)
    } else {
      0
    }
    
  } else if ((new_x-x)<0) {
    if ((new_x-x)>-1) {
      l = density(new_x)
      u = density(2*x-new_x)
      
      l/(u+l)
    } else {
      0
    }
  }
  
  
}

x = 0
sample = c()
burn_in = 10000
target_ss = 10000
n = 0
accept_count = 0

while (n<burn_in+target_ss) {
  n=n+1
  
  h = runif(1,-1,1)
  
  p = step_p(x,x+h)
  
  if (runif(1)<p) {
    new_x = x+h
  } else {
    new_x = x-h
  }
  
  r = (density(new_x)*step_p(new_x,x))/(density(x)*step_p(x,new_x))
  a = min(1,r)
  
  if (runif(1) < a) {
    x = new_x
    accept_count = accept_count + 1 
  }
  
  sample = append(sample,x)
}

hist(sample,freq=F,add=T,nclass=100)