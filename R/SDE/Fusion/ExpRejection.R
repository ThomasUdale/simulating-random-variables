f.phi <- function(x) {
  (sin(x)^2 + cos(x))/2 + 0.5
}

c <- 2

e.y <- exp(c)
e.inv <- exp(-c)

exp.unbiased <- function() {
  N <- rgeom(1,1-exp(-1))
  if (N==0){
    return(1)
  }
  z <- 1
  cur_x <- 1
  for (i in 1:N) {
    cur_x <- cur_x * rexp(1,1/c)/i
    z <- z + cur_x * exp(i)
  }
  return(z)
}

inv.unbiased <- function() {
  u <- runif(1)
  cur_l <- 0
  cur_u <- 1
  cur_est <- 0
  i <- 0
  while((u-cur_l)*(cur_u-u)>0){
    i <- i+1
    cur_est <- cur_est*(i-1)/i + 1/i*1/exp.unbiased()
    cur_l <- cur_est - 1/i * c^2
    cur_u <- cur_est
    print(c(i,u,cur_l,cur_u))
  }
  print('---')
  print(c(i,u,cur_l,cur_u))
  
  if (u<cur_l){
    return(1)
  } else {
    return(0)
  }
}

