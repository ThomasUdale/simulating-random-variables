library(layeredBB)

omega <- function(x,t) {
  (cos(x+t)+sin(x+t)^2+2*cos(x+t))/2
}

L <- -3/2
U <- 2

x0 <- 0

euler_approx <- function(t,n) {
  x <- x0
  elapsed <- 0 
  h <- t/n
  for (i in 1:n){
    elapsed <- elapsed + h
    x <- x + sin(x+elapsed)*h + rnorm(1,0,sqrt(h))
  }
  return(x)
}

ea1 <- function(t) {
  while(TRUE){
    k <- rpois(1,(U-L)*t)
    times <- runif(k,0,t)
    v <- runif(k,0,1)
    bb <- Brownian_motion_path_sampler(x0,c(0,times,t))
    phis <- omega(bb[1,2:(k+1)],times)
    if (all(v > (phis-L)/(U-L))) {
      return(bb[1,k+2])
    }
  }
} 


target <- pi
n <- 1e4

eu_out <- pbapply::pbreplicate(n,euler_approx(target,1000))
print(c(mean(eu_out),sd(eu_out)))
ea_out <- pbapply::pbreplicate(n,ea1(target))
print(c(mean(ea_out),sd(ea_out)))

hist(ea_out,freq=F)
lines(density(ea_out),col='red')
lines(density(eu_out))