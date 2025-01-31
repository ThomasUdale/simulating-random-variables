library(ggplot2)

pois_thinning <- function(){
  set.seed(2)
  t=2*pi
  
  k <- rpois(1,2*t)
  u <- runif(k,0,t)
  v <- runif(k,0,2)
  a <- 2-as.integer(v < sin(u)+1)
  times <- seq(0,t,length.out=100)
  lambda <- sin(times)+1
  
  plot <- ggplot()+
    geom_line(aes(x=times,y=lambda))+
    geom_ribbon(aes(x=times,ymin=0,ymax=lambda),fill='green',alpha=0.2)+
    geom_ribbon(aes(x=times,ymin=lambda,ymax=2),fill='red',alpha=0.2)+
    geom_point(aes(x=u,y=v),color=a)+
    labs(x="Time",y="Rate Function")+
    theme_classic()+
    annotate(geom = "text", x = 1.3, y = 1.72, label = expression((psi[1]*','*mu[1]*lambda)))+
    annotate(geom = "text", x = 2.6, y = 1.1, label = "Accept")+
    annotate(geom = "text", x = 3.36, y = 1.2, label = "Reject")+
    scale_y_continuous(breaks=c(0,2))+
    scale_x_continuous(breaks=c(0,pi, 2*pi),labels=c(0,expression(pi),expression(2*pi)))
  
  print(plot)
}

rejection_sampling <- function(){
  set.seed(1)
  times <- seq(0,1,length.out=100)
  target_density <- function(x) {dbeta(x,2,5)}
  M <- 2.5
  proposal_density <- function(x) {dunif(x)*M}
  
  plot <- ggplot()+
    geom_function(fun=target_density)+
    geom_function(fun=proposal_density)+
    geom_ribbon(aes(x=times,ymin=0,ymax=target_density(times)),fill='green',alpha=0.2)+
    geom_ribbon(aes(x=times,ymin=target_density(times),ymax=proposal_density(times)),fill='red',alpha=0.2)+
    annotate(geom = "text", x = 0.39, y = 1.1, label = "Accept")+
    annotate(geom = "text", x = 0.53, y = 1.2, label = "Reject")+
    labs(x="X",y="UMq(X)")+
    scale_y_continuous(breaks=c(0,2.5))+
    scale_x_continuous(breaks=c(0,1))+
    theme_classic()
  print(plot)
}

rejection_sampling()

