library(ggplot2)
library(progress)

mu <- 0
sigma <- 1

f <- function(x) {
  dnorm(x,mu,sqrt(sigma))
}



grad_log_f <- function(x) {
  -(x-mu)/sigma
}

pdmp <- function(np) {
  v <- rnorm(1)
  x <- rnorm(1,mu,sqrt(sigma))
  s <- 0
  event_count <- 1
  out <- matrix(nrow=np,ncol=4)
  out[1,] <- c(s,x,v,0)
  pb <- progress_bar$new(
    format = "sampling [:bar] :current/:total (:percent) eta: :eta, elapsed: :elapsed",
    clear = FALSE, total = np, width = 80)
  while (event_count<np) {
    pb$tick()
    event_count=event_count+1
    e <- rexp(1)
    v_event <- 1e9
    if (v*(x-mu)<0){
      event_time <- abs((x-mu)/v) + sqrt(2*e*sigma)
    } else {
      event_time <- -(x-mu)/v + sqrt(((x-mu)/v)^2 + 2*e*sigma/v^2)
    }
    # if (s+min(event_time,v_event)>t) {
    #   out <- rbind(out,c(t,x+(t-s)*v,v,0))
    #   pb$terminate()
    #   return(out)
    # }
    
    if (event_time<v_event){
      x <- x + event_time*v
      v_prop <- rnorm(1,-v,0.01)
      if (runif(1)<dnorm(v_prop)/dnorm(v)){
        v <- v_prop
      }
      s <- s + event_time
      out[event_count,] <-c(s,x,v,1)
    } else {
      v_event_count = v_event_count+1
      x <- x + v_event*v
      v <- rnorm(1)
      s <- s+v_event
      out[event_count,] <-c(s,x,v,2)
      # if (v_event_count==10){
      #   pb$terminate()
      #   return(out)
      # }
    }
  }
  pb$terminate()
  return(out)
}

interpolate <- function(out,t,m) {
  h <- t/m
  s <- h
  id <- 1
  idx <- 1
  int_out <- matrix(nrow=m,ncol=2)
  while (idx<m) {
    if (out[id,1]>=s) {
      int_out[idx,] <- c(s,out[id-1,2]+(s-out[id-1,1])*out[id-1,3])
      s <- s+h
      idx <- idx+1
    } else {
      id <- id+1
    }
  }
  return(int_out)
}

out <- pdmp(5000)
int_out <- interpolate(out,out[nrow(out),1],10000)

plot <- ggplot()+
  xlim(mu-(sqrt(sigma)*5),mu+(sqrt(sigma)*5))+
  geom_function(fun=f)+
  geom_density(aes(x=int_out[,2]),alpha=0.5,col='red')

print(plot)

plot <- ggplot()+
  geom_point(aes(x=out[,1],y=out[,2],colour=as.factor(out[,4])),alpha=0.4)+
  geom_line(aes(x=out[,1],y=out[,2]),alpha=0.1)

print(plot)