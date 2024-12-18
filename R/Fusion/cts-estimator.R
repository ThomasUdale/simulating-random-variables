library(ggplot2)
library(layeredBB)
# devtools::load_all("~/Documents/simulating-random-variables/R-packages/bm.rust")

phi <- function(x,c) {
  0.5*(4*x^6 - 6*x^2) + sqrt(2)
}

t = 1
s = 0
x = 0
y = 1


N <- rexp(1)
times <- c(seq(s,t,by=1/N),t)
bb <- brownian_bridge(x,y,s,t,times)
phi <- 

