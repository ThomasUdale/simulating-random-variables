library(ggplot2)
library(layeredBB)
# devtools::load_all("~/Documents/simulating-random-variables/R-packages/bm.rust")

# f.alpha <- function(x) {
#   -2*x^3
# }
# 
# f.phi <- function(x,c) {
#   0.5*(4*x^6 - 6*x^2) + sqrt(2)
# }
# 
# f.find_bound <- function(bes_layer) {
#   m <- max(c(abs(bes_layer$L),bes_layer$U))
#   if (m < (3/2)^(1/4)) {
#     f.phi(0)
#   } else {
#     f.phi(m)
#   }
# }

f.alpha <- function(x){
  sin(x)
}
f.phi <- function(x) {
  (sin(x)^2 + cos(x))/2 + 0.5
}

f.find_bound <- function(bes_layer){
  9/8
}

euler.skel.forward <- function(s,t,x,y,n){
  while(TRUE){
    cur_z_f <- x
    cur_z_b <- y
    for (i in 1:n) {
      cur_z_f <- append(cur_z_f,rnorm(1,mean=((cur_z_f[length(cur_z_f)]+f.alpha(cur_z_f[length(cur_z_f)])*((t-s)/n))),sd=sqrt((t-s)/n)))
      cur_z_b <- append(cur_z_b,rnorm(1,mean=((cur_z_b[length(cur_z_b)]+f.alpha(cur_z_b[length(cur_z_b)])*((t-s)/n))),sd=sqrt((t-s)/n)))
    }
    
    diff <- cur_z_f-rev(cur_z_b)
    
    for (i in 1:n) {
      if (diff[i]*diff[i+1] <0) {
        # print(i)
        # print(cur_z_f[1:i])
        # print(rev(cur_z_b)[i:n+1])
        return(c(cur_z_f[1:i],rev(cur_z_b)[i:n+1]))
      }
    }
  }
}

euler.likelihood <- function(s,t,x,y,n,skel) {
  l <- -(skel[1]-(x + f.alpha(x)*(t-s)/n)*(1-1/n) + y/n)^2/(2*(1-1/n)^2*(t-s)/n)
  l <- l + sum(sapply(2:(n-1),function(i) {-((skel[i]-(skel[i-1] + f.alpha(skel[i-1])*(t-s)/n)*(1-1/n) + y/n)^2)/(2*(1-1/n)^2*(t-s)/n)}))
  return(exp(l))
}

bb.var.off <- function(i,j,n) {
  if (i==j){
    i*(1-i/n)^2 + (n-i)*(i/n)^2
  } else if(i<j){
    i*(1-j/n)
  } else {
    j*(1-i/n)
  }
}

bb.var <- function(n) {
  pp <- 1/n*outer(1:(n-1),1:(n-1),FUN=Vectorize(bb.var.off),n=n)
  return(pp)
}

bb.likelihood <- function(s,t,x,y,n,skel) {
  pp <- bb.var(n)
  mu <- sapply(1:(n-1), function(i){(1-i/n)*x + (i/n)*y})
  l <- exp(-0.5*t(skel-mu) %*% solve(pp) %*% (skel-mu))
  return(l)
}



ea3 <- function(x,y,t) {
  while(TRUE){
    bes_layer <- bessel_layer_simulation(x = x, y = y, s = 0, t = t, mult = 0.2)
    rbound <- f.find_bound(bes_layer)
    k <- rpois(1,t*rbound)
    skeleton_u <- runif(k,0,rbound)
    skeleton_v <- runif(k,0,t)
    layered_bb <- layered_brownian_bridge(x = x,
                                          y = y,
                                          s = 0,
                                          t = t,
                                          bessel_layer = bes_layer,
                                          times = skeleton_v)
    if (all(sapply(layered_bb$simulated_path[1,],f.phi) < skeleton_u)) {
      return(t(layered_bb$full_path))
    }
  }
}

ea3.fill_in <- function(full_path,times) {
  cur_path_id <- 1
  cur_time_id <- 1
  bb <- c()
  cur_times <- c()
  while(length(bb)<length(times)) {
    if (times[cur_time_id]>=full_path[cur_path_id,2]) {
      if (length(cur_times)>0){
        bb <- append(bb,brownian_bridge(
          full_path[cur_path_id-1,1],
          full_path[cur_path_id,1],
          full_path[cur_path_id-1,2],
          full_path[cur_path_id,2],
          cur_times
        ))
        cur_times <- c()
      }
      if (times[cur_time_id]==full_path[cur_path_id,2]) {
        bb <- append(bb,full_path[cur_path_id,1])
        cur_time_id=cur_time_id+1
      }
      cur_path_id = cur_path_id+1
    } else {
      cur_times <- append(cur_times,times[cur_time_id])
      cur_time_id = cur_time_id+1
    }
  }
  return(bb)
}

ea3.approx <- function(x,y,t,n) {
  skel <- euler_skel(0,t,x,y,n)
  p <- bb.likelihood(0,t,x,y,n,skel)/euler.likelihood(0,t,x,y,n,skel)
  return(p)
  for (i in 1:n) {
    if (i==1) {
      lower = x
      upper = skel[i]
    } else if (i==n) {
      lower = skel[i-1]
      upper=y
    } else {
      lower = skel[i-1]
      upper = skel[i]
    }
    bes_layer <- bessel_layer_simulation(x = lower, y = upper, s = (i-1)*t/n, t = i*t/n, mult = 0.2)
    rbound <- f.find_bound(bes_layer)
    k <- rpois(1,t*rbound/n)
    skeleton_u <- runif(k,0,rbound)
    skeleton_v <- runif(k,(i-1)*t/n,i*t/n)
    layered_bb <- layered_brownian_bridge(x = lower, 
                                          y = upper, 
                                          s = (i-1)*t/n, 
                                          t = i*t/n,
                                          bessel_layer = bes_layer,
                                          times = skeleton_v)
    p <- p * prod((rbound - f.phi(layered_bb$simulated_path[1,]))/rbound)
  }
  p
}

s=0
t=1
x=0
y=0
n=100

r_n <- 1e4

# e.t <- seq(s,t,length.out=n+1)
# ea.y <- pbapply::pbreplicate(r_n,ea3.fill_in(ea3(x,y,t),e.t))
# e.y <- pbapply::pbreplicate(r_n,euler.skel.forward(s,t,x,y,n))
# b.y <- pbapply::pbreplicate(r_n,brownian_bridge(x,y,s,t,e.t))

qu = 0.95
ql = 0.05
e.u=apply(e.y,1,quantile,probs=qu)
e.l=apply(e.y,1,quantile,probs=ql)
b.u=apply(b.y,1,quantile,probs=qu)
b.l=apply(b.y,1,quantile,probs=ql)
ea.u=apply(ea.y,1,quantile,probs=qu)
ea.l=apply(ea.y,1,quantile,probs=ql)
df <- data.frame(
  x=e.t,
  e.y=rowMeans(e.y),
  b.y=rowMeans(b.y),
  ea.y=rowMeans(ea.y),
  e.u=e.u,
  e.l=e.l,
  b.u=b.u,
  b.l=b.l,
  ea.u=ea.u,
  ea.l=ea.l
)


plot <-
  ggplot(data=df,aes(x=x))+
  geom_line(aes(y=e.y))+
  geom_line(aes(y=b.y),col='blue')+
  geom_line(aes(y=ea.y),col='red')+
  geom_ribbon(aes(ymin=e.l,ymax=e.u),col='black',alpha=0.2)+
  geom_ribbon(aes(ymin=b.l,ymax=b.u),col='blue',alpha=0.2)+
  geom_ribbon(aes(ymin=ea.l,ymax=ea.u),col='red',alpha=0.2)
print(plot)

# bb <- replicate(1e4, brownian_bridge(0,1,0,1,0.5))
# eout <- replicate(1e4,euler.skel(0,1,0,1,2))
# 
# 
# plot <-
#   ggplot()+
#   geom_density(data=data.frame(x=bb),aes(x=x),col='red')+
#   geom_density(data=data.frame(x=eout),aes(x=x),col='blue')+
#   xlim(-2,2)
# print(plot)
