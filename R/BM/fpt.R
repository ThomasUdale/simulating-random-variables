# devtools::load_all("~/Documents/simulating-random-variables/R-packages/cts.smc")



g <- function(k,t) {
  (-1)^k*(1+2*k)/sqrt(2*pi*t^3) * exp(-(1+2*k)^2/(2*t))
}

sn <- function(n,t) {
  sum(sapply(c(-n:n),g,t=t))
}



sim <- function() {
  while(TRUE){
    x <- rgamma(1,1.088870,1.233701)
    w <- runif(1)*1.243707*dgamma(x,1.088870,1.233701)
    n <- max(ceiling(0.275*x),3)
    if (n%%2==0){
      n = n+1
    }
    s <- sn(n,x)
    resolved=FALSE
    while(!resolved) {
      n <- n+2
      s_step <- g(-n,x) + g(n,x) + g(-(n-1),x)+g(n-1,x)
      s <- s + s_step
      if (abs(w-s)>s_step) {
        resolved = TRUE
        if (w<s) {
          return(x)
        } 
      }
    }
  }
}



g2 <- function(k,x) {
  (-1)^k*2*(2*k+1)*exp(-((2*k+1)^2)/(2*x))
}

sn2 <- function(x) {
  sapply(x,function(x) {sum(sapply(c(0:2),g2,x=x))}*(1/sqrt(2*pi*x^3)))
}

sn3 <- function(x) {
  2*(dinvgamma(x,1/2,rate=1/2)+dinvgamma(x,1/2,rate=25/2)-dinvgamma(x,1/2,rate=9/2))
}


ss = 1e4

e <- brownian_motion_fpt(ss)
hist(e,freq=F,ylim=c(0,1))
lines(density(e),col='green')

e <- replicate(ss,sim())
lines(density(e),col='red')

# curve(sn2,add=T,xlim=c(1e-4,6),col='orange')
curve(sn3,add=T,xlim=c(1e-4,6),col='blue')


