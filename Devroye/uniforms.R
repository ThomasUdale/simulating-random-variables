
produnif <- function() {
  while(TRUE){
    u <- runif(1,0,1)
    v <- runif(1,0,1.5/exp(1))
    y <- v/(u^2)
    if(u^3<dexp(y)) {
      return(y)
    }
  }
}