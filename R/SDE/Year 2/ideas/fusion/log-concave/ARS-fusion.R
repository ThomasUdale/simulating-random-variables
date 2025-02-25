library(ggplot2)

f <- function(x) {
  dnorm(x,0,sqrt(0.5))
}

mu = c(-5,5)
sigma = c(1,1)

fc <- function(x,c) {
  exp(log_fc(x,c))
}

log_fc <- function(x,c){-0.5*(x-mu[c])^2/sigma[c]^2}

samplefc <- function(c) {
  rnorm(1,mu[c],sigma[c])
}



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

combine_lines <- function(lines1,lines2) {
  cur_x <- -Inf 
  cur_id1 <- 1
  cur_id2 <- 1
  out_lines <- c()
  while(cur_id1<=nrow(lines1) | cur_id2<=nrow(lines2)) {
    if (lines1[cur_id1,2]<lines2[cur_id2,2]) {
      out_lines <- rbind(
        out_lines, 
        c(
          cur_x,
          lines1[cur_id1,2],
          lines1[cur_id1,3]+lines2[cur_id2,3],
          lines1[cur_id1,4]+lines2[cur_id2,4]
        )
      )
      cur_x = lines1[cur_id1,2]
      cur_id1 = cur_id1+1
      
    } else if (lines1[cur_id1,2]>lines2[cur_id2,2]){
      out_lines <- rbind(
        out_lines, 
        c(
          cur_x,
          lines2[cur_id2,2],
          lines1[cur_id1,3]+lines2[cur_id2,3],
          lines1[cur_id1,4]+lines2[cur_id2,4]
        )
      )
      cur_x = lines2[cur_id2,2]
      cur_id2 = cur_id2+1
    } else {
      out_lines <- rbind(
        out_lines, 
        c(
          cur_x,
          lines2[cur_id2,2],
          lines1[cur_id1,3]+lines2[cur_id2,3],
          lines1[cur_id1,4]+lines2[cur_id2,4]
        )
      )
      cur_x = lines2[cur_id2,2]
      if (lines2[cur_id2,2]==Inf){
        return(out_lines)
      }
      cur_id2=cur_id2+1
      cur_id1=cur_id1+1
    }
  }
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

a.r.fusion <- function(x1,x2){
  cur_x1 <- x1
  cur_x2 <- x2
  evals <- 0
  while(TRUE){
    jl1 <- joining_lines(cur_x1)
    jl2 <- joining_lines(cur_x2) 
    
    el1 <- envelope_lines(jl1)
    el2 <- envelope_lines(jl2)
    el <- combine_lines(el1,el2)
    prop <- sim_g(el)
    u <- runif(1)
    if (u<exp(envelope_l(prop,jl1)+envelope_l(prop,jl2))/exp(envelope_u(prop,el))) {
      return(list(sample=prop,x1=cur_x1,x2=cur_x2,evals=evals))
    } else {
      evals <- evals+1
      f_prop <- prod(sapply(c(1,2),fc,x=prop))
      cur_x1 <- rbind(cur_x1,c(prop,log_fc(prop,1)))
      cur_x1 <- cur_x1[order(cur_x1[,1]),]
      cur_x2 <- rbind(cur_x2,c(prop,log_fc(prop,2)))
      cur_x2 <- cur_x2[order(cur_x2[,1]),]
      
      if (u<f_prop/exp(envelope_u(prop,el))) {
        return(list(sample=prop,x1=cur_x1,x2=cur_x2,evals=evals))
      }
    }
  }
}

a.fusion <- function(n) {
  s1 <- sort(rnorm(1e1,mu[1],sigma[1]))
  x1 <- cbind(s1,sapply(s1,log_fc,c=1))
  s2 <- sort(rnorm(1e1,mu[2],sigma[2]))
  x2 <- cbind(s2,sapply(s2,log_fc,c=2))
  
  cur_x1 <- x1
  cur_x2 <- x2
  out <- c()
  evals <- 0
  
  for (i in 1:n){
    s <- a.r.fusion(cur_x1,cur_x2)
    cur_x1 <- s$x1
    cur_x2 <- s$x2
    out <- append(out,s$sample)
    if (s$evals>0){
      evals <- evals+s$evals
      print(c(i,evals))
    }
  }
  
  return(list(sample=out))
}


for (n in c(1e4)){
  out <- a.fusion(n)
}


plot <- ggplot()+
  xlim(mu[1]-3,mu[2]+3)+
  ylim(-2,1)+
  geom_function(fun=log_fc,args=list(c=1),col='blue',alpha=0.2)+
  geom_function(fun=fc,args=list(c=1),col='blue',alpha=0.2)+
  geom_function(fun=log_fc,args=list(c=2),col='blue',alpha=0.2)+
  geom_function(fun=fc,args=list(c=2),col='blue',alpha=0.2)+
  geom_function(fun=f)+
  geom_density(aes(x=out$sample),color="red",alpha=0.5)

print(plot)