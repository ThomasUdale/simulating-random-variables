mu <- mean(replicate(1e5,exp(-rexp(1))))

f <- function(x,n) {
  if (n/2<x){
    n/2
  } else {
    x
  }
}

f_prime <- function(x,n) {
  if (n/2<x) {
    0.5
  } else {
    0
  }
}

x_prime <- function(x,n) {
  (f(x,n)/n-f_prime(x,n))*(1-f(x,n)/n)^(n-1)
}

sim <- function() {
  N <- rexp(1,1)
  x <- rexp(1)
  Z <- x_prime(x,N)/dexp(N)
  if (Z>10){
    print(c(N,x,Z))
  }
  return(Z)
}

out <- pbapply::pbreplicate(1e5,sim())
print(summary(out))
hist(out[out<=1],freq=F,xlim=c(0,1),nclass =100)