f.phi <- function(x) {
  (sin(x)^2 + cos(x))/2 + 0.5
}

e.y <- exp(4)
e.inv <- exp(-4)

exprej <- function() {
  N <- rgeom(1,1-exp(-1))
  if(N==0){
    return(1)
  }
  z <- 1
  for (i in 1:N){
    z <- z+ exp(i)*prod(rexp(i,1/4))/factorial(i)
  }
  e <- rexp(1,z)
  return(e)
}

