library(ggplot2)

gen_prop <- function(n){
  rexp(n,1)
}

gen_d <- function(x){
  dexp(x,1)
}

tar_d <- function(x){
  dchisq(x,4)
}

weight <- function(x) {
  tar_d(x)/gen_d(x)
}

est_n <- function(n) {
  x <- gen_prop(n)
  w <- weight(x)
  mean(x*w)/mean(w)
}

# look at reduction in variance
# for (i in 0:4){
#   es <- pbapply::pbreplicate(1e4,est_n(10^i))
#   print(summary(es))
#   print(sd(es))
# }

unbiased <- function() {
  N <- rgeom(1,1-3^(-1))
  # if (N>6){
  #   print(N)
  # }
  x <- gen_prop(1)
  w <- weight(x)
  I <- mean(w*x)/mean(w)
  Z <- I
  if (N==0){
    return(Z)
  }
  for (i in 1:N){
    # print(i)
    new_x <- gen_prop(floor(10^log(i+1)-10^log(i)))
    x <- append(x,new_x)
    w <- append(w,weight(new_x))
    new_I <- mean(w*x)/mean(w)
    Z <- Z + (new_I-I)*exp(i)
    I <- new_I
    # print(c(100*i,new_I,Z))
  }
  return(Z)
}
# 
est <- pbapply::pbreplicate(1e5,unbiased())
print(summary(est))
print(sd(est))
biased_est <- pbapply::pbreplicate(1e5,est_n(1000))
print(summary(biased_est))
print(sd(biased_est))

# 
# plot <- ggplot()+
#   geom_histogram(aes(x=w,y=..density..))
# 
# print(plot)