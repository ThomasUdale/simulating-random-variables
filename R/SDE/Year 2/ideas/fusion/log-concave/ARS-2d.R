library(ggplot2)
library(mvnfast)
library(geometry)

# f <- function(x){dnorm(x)}
# 
# log_f <- function(x){-0.5*x^2}
# 
# x2 <- sort(rnorm(10))
# log_fx <- log_f(x2)
# 
# x <- cbind(x2,log_fx)
# 
# joining_lines <- function(x){
#   lines <- c()
#   for (i in 2:nrow(x)){
#     m <- (x[i,2]-x[i-1,2])/(x[i,1]-x[i-1,1])
#     a <- -m*x[i,1]+x[i,2]
#     lines <- rbind(lines,c(x[i-1],x[i],a,m))
#   }
#   return(lines)
# }
# 
# envelope_lines <- function(lines){
#   lu <- c(-Inf,lines[1,1],lines[1,3],lines[1,4])
#   for (l in 1:nrow(lines)) {
#     if (l==1){
#       lu <- rbind(lu,c(lines[l,1:2],lines[l+1,3:4]))
#     } else if (l==nrow(lines)) {
#       lu <- rbind(lu,c(lines[l,1:2],lines[l-1,3:4]))
#       lu <- rbind(lu,c(lines[l,2],Inf,lines[l,3:4]))
#     } else {
#       mid <- -(lines[l+1,3]-lines[l-1,3])/(lines[l+1,4]-lines[l-1,4])
#       lu <- rbind(lu,c(lines[l,1],mid,lines[l-1,3:4]))
#       lu <- rbind(lu,c(mid,lines[l,2],lines[l+1,3:4]))
#     }
#   }
#   return(lu)
# }
# 
# envelope_u <- function(x,e_lines) {
#   for (l in 1:nrow(e_lines)) {
#     if (x<e_lines[l,2]) {
#       return(e_lines[l,3]+e_lines[l,4]*x)
#     }
#   }
# }
# 
# envelope_l <- function(x,lines) {
#   if (x<lines[1,1]){
#     return(-Inf)
#   }
#   for (l in 1:nrow(lines)) {
#     if (x<lines[l,2]) {
#       return(lines[l,3]+lines[l,4]*x)
#     }
#   }
#   return(-Inf)
# }
# 
# sim_g <- function(e_lines){
#   probs <- apply(e_lines,1,function(x){exp(x[3])/x[4]*(exp(x[4]*x[2])-exp(x[4]*x[1]))})
#   it <- sample(1:nrow(e_lines),1,prob=probs)
#   u <- runif(1)
#   x <- 1/e_lines[it,4]*log(
#     exp(
#       e_lines[it,4]*e_lines[it,1])+
#       u*(exp(e_lines[it,4]*e_lines[it,2])-exp(e_lines[it,4]*e_lines[it,1]))
#   )
#   return(x)
# }
# 
# a.r <- function(x){
#   cur_x <- x
#   evals <- 0
#   while(TRUE){
#     jl <- joining_lines(cur_x)
#     el <- envelope_lines(jl)
#     prop <- sim_g(el)
#     u <- runif(1)
#     if (u<exp(envelope_l(prop,jl))/exp(envelope_u(prop,el))) {
#       return(list(sample=prop,x=cur_x,evals=evals))
#     } else {
#       f_prop <- f(prop)
#       if (u<f_prop/exp(envelope_u(prop,el))) {
#         return(list(sample=prop,x=cur_x,evals=evals))
#       } else {
#         evals <- evals+1
#         cur_x <- rbind(cur_x,c(prop,log_f(prop)))
#         cur_x <- cur_x[order(cur_x[,1]),]
#       }
#     }
#   }
# }
# 
# a.r.s <- function(n,x) {
#   cur_x <- x
#   out <- c()
#   evals <- 0
#   for (i in 1:n) {
#     s <- a.r(cur_x)
#     cur_x <- s$x
#     out <- append(out,s$sample)
#     if (s$evals>0){
#       evals <- evals+s$evals
#       print(c(i,evals))
#     }
#   }
#   return(list(sample=out,cur_x=cur_x))
# }
# 
# out <- a.r.s(1e4,x)

f <- function(x){
  dmvn(x,c(0,0),sigma=matrix(c(1,3/5,3/5,2),nrow=2))
}

log_f <- function(x){
  log(f(x))
}

x1 <- seq(-4,4,0.1)
x2 <- seq(-4,4,0.1)

f_grid <- expand.grid(X1=x1,X2=x2)
f_grid$f <- apply(f_grid,1,log_f)

x <- rmvn(10,c(0,0),sigma=matrix(c(1,3/5,3/5,2),nrow=2))
fx <- log_f(x)
hull <- melt(convhulln(cbind(x,fx)))
hull$x <- x[hull$value,1]
hull$y <- x[hull$value,2]


plot <- ggplot()+
  geom_contour_filled(data=f_grid,aes(X1,X2,z=f))+
  geom_polygon(data=hull,aes(x=x,y=y,group=Var1,fill=as.factor(Var1)))+
  geom_point(aes(x=x[,1],y=x[,2]))
 

print(plot)
