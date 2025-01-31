# devtools::load_all("~/Documents/simulating-random-variables/R-packages/bm.rust")

temp = 1
x = 0
y = 1

phi <- function(x) {
  0.5*(4*x^6 - 6*x^2) + sqrt(2)
}

find_bound <- function(bes_layer) {
  m <- max(c(abs(bes_layer$L),bes_layer$U))
  if (m < (3/2)^(1/4)) {
    phi(0)
  } else {
    phi(m)
  }
}


approx_integral <- function(delta,path){
  exp(-delta*sum(phi(path)))
}

particle_filter <- function(x,y,t,N,SS){
  step <- 1
  paths <- matrix(NaN,nrow=N+1,ncol=2)
  paths[,1]=x
  paths[,2]=y
  paths[1,]=c(0,t)
  print(paths)
  for (i in 1:SS){
    step <- step+1
    new_t <- runif(1,0,t)
    u_t <- min(which(paths[1,]>new_t))
    l_t <- u_t-1
    paths <- cbind(paths[,1:l_t],NaN,paths[,u_t:step])
    paths[1,u_t] <- new_t
    for (j in 1:N){
     paths[j+1,u_t] <- brownian_bridge(paths[j+1,l_t],paths[j+1,u_t+1],paths[1,l_t],paths[1,u_t+1],new_t)
    }
  }
}

particle_filter(0,1,2,5,3)

