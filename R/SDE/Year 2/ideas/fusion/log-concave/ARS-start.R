library(ggplot2)

f <- function(x){dnorm(x)}

log_f <- function(x){-0.5*x^2}

x2 <- sort(rnorm(10))
log_fx <- log_f(x2)

x <- cbind(x2,log_fx)

joining_lines <- function(x){
  lines <- c()
  for (i in 2:nrow(x)){
    m <- (x[i,2]-x[i-1,2])/(x[i,1]-x[i-1,1])
    a <- -m*x[i,1]+x[i,2]
    lines <- rbind(lines,c(x[i-1],x[i],a,m))
  }
  return(lines)
}

envelope_lines <- function(lines){
  lu <- c(-Inf,lines[1,1],lines[1,3],lines[1,4])
  for (l in 1:nrow(lines)) {
    if (l==1){
      lu <- rbind(lu,c(lines[l,1:2],lines[l+1,3:4]))
    } else if (l==nrow(lines)) {
      lu <- rbind(lu,c(lines[l,1:2],lines[l-1,3:4]))
      lu <- rbind(lu,c(lines[l,2],Inf,lines[l,3:4]))
    } else {
      mid <- -(lines[l+1,3]-lines[l-1,3])/(lines[l+1,4]-lines[l-1,4])
      lu <- rbind(lu,c(lines[l,1],mid,lines[l-1,3:4]))
      lu <- rbind(lu,c(mid,lines[l,2],lines[l+1,3:4]))
    }
  }
  return(lu)
}

envelope_u <- function(x,e_lines) {
  for (l in 1:nrow(e_lines)) {
    if (x<e_lines[l,2]) {
      return(e_lines[l,3]+e_lines[l,4]*x)
    }
  }
}

envelope_l <- function(x,lines) {
  if (x<lines[1,1]){
    return(-Inf)
  }
  for (l in 1:nrow(lines)) {
    if (x<lines[l,2]) {
        return(lines[l,3]+lines[l,4]*x)
    }
  }
  return(-Inf)
}

sim_g <- function(e_lines){
  probs <- apply(e_lines,1,function(x){exp(x[3])/x[4]*(exp(x[4]*x[2])-exp(x[4]*x[1]))})
  it <- sample(1:nrow(e_lines),1,prob=probs)
  u <- runif(1)
  x <- 1/e_lines[it,4]*log(
    exp(
      e_lines[it,4]*e_lines[it,1])+
      u*(exp(e_lines[it,4]*e_lines[it,2])-exp(e_lines[it,4]*e_lines[it,1]))
  )
  return(x)
}

a.r <- function(x){
  cur_x <- x
  evals <- 0
  while(TRUE){
    jl <- joining_lines(cur_x)
    el <- envelope_lines(jl)
    prop <- sim_g(el)
    u <- runif(1)
    if (u<exp(envelope_l(prop,jl))/exp(envelope_u(prop,el))) {
      return(list(sample=prop,x=cur_x,evals=evals))
    } else {
      evals <- evals+1
      f_prop <- f(prop)
      cur_x <- rbind(cur_x,c(prop,log_f(prop)))
      cur_x <- cur_x[order(cur_x[,1]),]
      if (u<f_prop/exp(envelope_u(prop,el))) {
        return(list(sample=prop,x=cur_x,evals=evals))
      }
    }
  }
}

a.r.s <- function(n,x) {
  cur_x <- x
  out <- c()
  evals <- 0
  for (i in 1:n) {
    s <- a.r(cur_x)
    cur_x <- s$x
    out <- append(out,s$sample)
    if (s$evals>0){
      evals = evals+s$evals
      print(c(i,evals))
    }
  }
  return(list(sample=out,cur_x=cur_x))
}

out <- a.r.s(1e4,x)


plot <- ggplot()+
  xlim(-4,4)+
  ylim(-2,1)+
  geom_function(fun=log_f)+
  geom_point(aes(x=x2,y=log_fx))+
  geom_function(fun=Vectorize(envelope_u,vectorize.args=c("x")),args=list(e_lines=envelope_lines(joining_lines(x))),col='blue')+
  geom_function(fun=Vectorize(envelope_l,vectorize.args=c("x")),args=list(lines=joining_lines(x)),col='green')+
  geom_function(fun=f)+
  geom_density(aes(x=out$sample),col='red')

print(plot)