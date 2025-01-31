use("pseudo.fusion")

target_draw = list(
  d1 = function() rnorm(1,2,5),
  d2 = function() rt(1,df=1),
  d3 = function() rcauchy(1,-2,1)
)


target_dists = list(
  d1 = function(x) dnorm(x,2,5),
  d2 = function(x) dt(x,df=1),
  d3 = function(x) dcauchy(x,-2,1)
)
num_dist = length(target_dists)

call_fun <- function(f, ...) f(...)
prod_of_norms <- function(x) {
  Reduce(lapply(target_dists, call_fun, x),f="*")
}

q = integrate(prod_of_norms,lower=-Inf,upper=Inf)$value

f <- function(x) {
  prod_of_norms(x)/q
}


curve(f,-10,10,ylim=c(0,0.6))

y = 0
r = 1e-10
samples = c()
n = 0
accepts = 0
# burn = 20000
# target= 20000
# num_is = 10000
# 
# sample_rust = mixture_fusion_example(burn,target,num_is)
# hist(sample_rust,col='red',add=T,freq=F,nclass=100)

burn = 10000
target= 10000
num_is = 100

pb = txtProgressBar(min = 0, max = (burn+target), initial = 0,style=3)

while (n<(burn+target)) {
  n = n+1
  setTxtProgressBar(pb,n)

  ynew = rnorm(1,y)

  new_r = 0
  for (i in 1:num_is) {
    r_star = 1
    for (c in 1:num_dist) {
      x = target_draw[[c]]()
      r_star = r_star*min((target_dists[[c]](ynew)/target_dists[[c]](x))*dnorm(x,ynew),dnorm(ynew,x))
    }
    new_r = new_r + r_star/num_is
  }

  a = new_r/r

  if (runif(1)<a){
    accepts = accepts+1
    y = ynew
    r = new_r
  }

  if (n>burn) {samples = append(samples,y)}


}

close(pb)
hist(samples,freq=F,add=T,nclass=50)