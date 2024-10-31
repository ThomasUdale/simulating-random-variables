f <- function(x) {
  dnorm(x,0,sqrt(0.5))
}

mu_0 = 3
sub_means = c(0.5*mu_0,-0.5*mu_0)
C = length(sub_means)

fc <- function(x,c) {
  dnorm(x,sub_means[c],1)
}

samplefc <- function(n,c) {
  rnorm(n,sub_means[c],1)
}

sample_p_fc <- function(x,c) {
  y <- rnorm(1,x)
  if (runif(1)<fc(y,c)/fc(x,c)){
    return(y)
  } else {
    return(x)
  }
}

total_time <- 10000
n_particles <- 1000

x <- array(rep(NaN,(total_time+1)*n_particles*C),dim=c((total_time+1),n_particles,C))
weights <- array(rep(NaN,(total_time+1)*n_particles),dim=c((total_time+1),n_particles))

x[1,,] <- array(unlist(lapply(c(1:C),samplefc,n=n_particles)),dim=c(n_particles,C))
weights[1,] <- 1

pb <- progress_bar$new(
      format = "sampling [:bar] :current/:total (:percent) eta: :eta, elapsed: :elapsed",
      clear = FALSE, total = total_time, width = 80)

for (i in 1:(total_time)){
  pb$tick()
  id_new <- sample(c(1:n_particles),n_particles,prob=weights[i,],replace=TRUE)
  for (n in 1:n_particles){
    x[i+1,n,] <- mapply(sample_p_fc,x[i,id_new[n],],c(1:C))
    weights[i+1,n] <- exp(-var(x[i+1,n,]))
    print(min(weights[i+1,n]))
  }
  
}
pb$terminate()

out <- rowMeans(x[total_time+1,,])

hist(out,freq=F,xlim=c(-3,3),ylim=c(0,0.8),nclass=50)
curve(f,add=T,col='red')



