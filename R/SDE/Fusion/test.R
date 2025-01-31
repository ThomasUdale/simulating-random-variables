N <- 1000

q <- 0.5



sample_geom <- function(){
  sample <- 0
  
  for (i in 1:N){
    x <- rgeom(1,prob=rnorm(1,mean=q,sd=0.05))+1
    if (sum(rexp(x))>1){
      sample <- sample + 1
    }
  }
  
  sample/N
}

sample_geom_bounded <- function() {
  prod(replicate(8,sample_geom()))
}


print(sample_geom())
print(exp(-q))

print(sample_geom_bounded())
print(exp(-q)^8)
