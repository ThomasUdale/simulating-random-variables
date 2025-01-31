mu = 0.05
sigma = 0.2
x_0 = 1

f <- function(x) {
  exp(-mu)*max((x-1),0)
}

exact <- function() {
  x_0 * exp((mu-sigma^2/2)*1 + sigma*rnorm(1))
}

o_e <- replicate(1e4,f(exact()))

biased <- function(n) {
  x <- x_0
  for (i in 1:2^n) {
    x <- x + mu * x * 1/2^n + sigma * x * rnorm(1,0,(1/2^n)^2)
  }
  return(x)
}

o_b <- replicate(1e4,f(biased(5)))