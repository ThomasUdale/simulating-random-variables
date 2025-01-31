library(ggplot2)
library(layeredBB)

target <- 0.33


phi <- function(x) {
  0.5*(4*x^6 - 6*x^2) + sqrt(2)
}

alpha <- function(x) {
  -2*x^3
}

euler_sim <- function(x,y,t,n){
  steps <- seq(0,t,length.out=n+1)
  l <- 1
  path <- c(l)
  bb_path <- Brownian_bridge_path_sampler(x,y,0,t,steps)$full_path
  for (i in 1:n){
    l <- l+l*alpha(bb_path[1,i+1])*rnorm(1,0,(t/n)^2)
    path <- append(path,l)
  }
  print(path)
}

euler_sim(0,1,1,10)