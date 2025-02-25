library(ggplot2)
library(ggpmisc)
# devtools::load_all("~/Documents/simulating-random-variables/R-packages/bm.rust")


x <- 0
y <- 1
s <- 0
t <- 1

phi <- function(x) {
  0.5 * ((x)^2-1) + 0.5
}

cur_length <- 1
phis = c(phi(x),phi(y))
cur_I = (t-s)/(cur_length+1)*sum(phis)
times <- (0:cur_length)/cur_length*(t-s) + s
path = rbind(brownian_bridge(x,y,s,t,times),times,phis,rep(0,2))
graph_path <- path

integrals <- c(cur_length,round(cur_I,4))

for (j in 1:6) {
  new_points <- matrix(nrow=4,ncol=cur_length)
  for (i in 1:cur_length) {
    # sim new point
    new_x <- brownian_bridge(
      path[1,i],
      path[1,i+1],
      path[2,i],
      path[2,i+1],
      (path[2,i]+path[2,i+1])/2
    ) 
    new_points[,i] <- c(new_x,(path[2,i]+path[2,i+1])/2,phi(new_x),j)
  }
  phis <- append(phis,phi(new_points[1,]))
  new_path <- matrix(nrow=4,ncol=(cur_length*2+1))
  new_path[,(1:(cur_length*2+1))[(1:(cur_length*2+1))%%2==1]] = path
  new_path[,(1:(cur_length*2+1))[(1:(cur_length*2+1))%%2==0]] = new_points
  path <- new_path
  new_path[4,] <- j
  graph_path <- cbind(graph_path,new_path)
  cur_length = cur_length*2
  new_I <- 0.5*(t-s)/cur_length*(2*sum(phis)-(phis[1]+phis[2]))
  integrals <- rbind(integrals,c(j,round(new_I,4)))
}

integrals <- data.frame(integrals)
colnames(integrals) <- c("level","integral estimate")

graph_path <- data.frame(t(graph_path))
colnames(graph_path) <- c("x","t","phi",'level')

plot <- ggplot(graph_path)+
  geom_line(aes(x=t,y=phi,color=as.factor(level)))+
  annotate(
    geom='table',
    x=0.05,
    y=1,
    label=list(integrals)
  )

print(plot)

