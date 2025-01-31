library(layeredBB)
library(ggplot2)
library(pbmcapply)
library(progress)

devtools::load_all("~/Documents/simulating-random-variables/R-packages/cts.smc")
# set.seed(1)



l = -0.5
r = 9/8
phi <- function(x) {
  (sin(x-pi)^2 + cos(x-pi))/2
}

sim_end <- function(t) {
  while (T) {
    u <- rnorm(1,0,sqrt(t))
    if (runif(1) < exp(-cos(u-pi)-1)) {
      return(u)
    }
  }
}

po_est_bound <- function(path,tt) {
  k <- rpois(1,r*tt)
  
  if(k==0){
    return(list(rho=1,path=path))
  }
  
  skeleton = sort(runif(k,0,tt))
  end_inds <- c()
  for (s in skeleton){
    end_inds <- append(end_inds, which(path$t>=s,arr.ind=T)[1])
  }
  
  skel_points <- c()
  for (id in unique(end_inds)){
    bb <- Brownian_bridge_path_sampler(
      x=path$w[id-1],
      y=path$w[id],
      s=path$t[id-1],
      t=path$t[id],
      times=skeleton[end_inds==id]
    )
    skel_points <- append(skel_points,bb$simulated_path[1,])
  }

  return(list(rho=prod((r - phi(skel_points)+l)/r),path=rbind(path,data.frame(t=skeleton,w=skel_points))))
}


p_smc <- function(N,n_steps,tt) {
  initial_paths <- replicate(N,sim_end(tt))
  
  paths <- list()
  
  for (n in 1:N) {
    paths[[length(paths)+1]] <- data.frame(t=c(0,tt),w=c(0,initial_paths[n]))
  }
  
  for (t_step in 1:n_steps) {
    print(t_step)
    # mutation
    for (n in 1:N){
      times <- diff(paths[[n]]$t)
      update_noise <- rnorm(length(times)-1,0,0.01)
      new_w <- paths[[n]]$w + c(0,update_noise)
      if (runif(1) < prod(dnorm(new_w,c(0,new_w[1:(length(times)-1)]),sqrt(times)))/prod(dnorm(paths[[n]]$w,c(0,paths[[n]]$w[1:(length(times)-1)]),sqrt(times)))){
        paths[[n]]$w <- new_w
      }
    }
    
    # selection
    weights <- c()
    for (n in 1:N) {
      out <-po_est_bound(paths[[n]],tt)
      paths[[n]] <- out$path
      paths[[n]] <- paths[[n]][with(paths[[n]], order(t)),] 
      row.names(paths[[n]]) <- NULL
      weights = append(weights,out$rho)
    }
    
    old_paths <- paths
    for (n in 1:N){
        paths[n] <- old_paths[sample(N,1,prob=weights)]
    }
    print(c(t_step,mean(weights),min(weights)))
  }
  return(list(paths=paths,weights=weights))
}



sample_at_t <- function(outp,sample_t){
  pb <- progress_bar$new(
    format = "sampling [:bar] :current/:total (:percent) eta: :eta, elapsed: :elapsed",
    clear = FALSE, total = length(outp$weights), width = 80)
  samples <- c()
  for (path in outp$paths) {
    pb$tick()
    id <- which(path$t>=sample_t,arr.ind=T)[1]
    bb <- Brownian_bridge_path_sampler(
      x=path$w[id-1],
      y=path$w[id],
      s=path$t[id-1],
      t=path$t[id],
      times=sample_t
    )
    samples <- append(samples,bb$simulated_path[1])
  }
  pb$terminate()
  samples
}






ea1 <- function(id,t,return_path=F){
  while (T){
    y <- sim_end(t)
    k <- rpois(1,r*t)
    if (k==0){
      if (return_path){
        return(data.frame(t=c(0,t),w=c(0,y)))
      } else {
        return(y)
      }
    }
    skel_t <- sort(runif(k,0,t))
    skel_u <- runif(k,0,r)
    bb <- Brownian_bridge_path_sampler(
      x=0,
      y=y,
      s=0,
      t=t,
      times=skel_t
    )
    if (all((phi(bb$simulated_path[1,])-l)<skel_u)){
      if (return_path){
        return(data.frame(t=bb$full_path[2,],w=bb$full_path[1,]))
      } else {
        return(y)
      }
    }
  }
}

ss= 100000
interval_time = 3
target_time = 1

# discrete-path-smc
# start_time <- proc.time()
# out_paths <- p_smc(N=ss,n_steps=1,tt=interval_time)
# print(proc.time()-start_time)
# out_sample <- data.frame(w=sample_at_t(out_paths,target_time),weights=out_paths$weights)

# rust_smc
start_time <- proc.time()
out_sample_r <- path_smc_rust(n=ss,n_steps=5,tt=interval_time,target_time)
print(proc.time()-start_time)
out_sample_r <- data.frame(out_sample_r)

# exact algorithm

start_time <- proc.time()
exact = pbmclapply(c(1:100000),ea1,t=target_time,mc.cores=8)
print(proc.time()-start_time)
exact <- data.frame(w=unlist(exact))



plot <- ggplot()+
  geom_density(data=out_sample_r,aes(x=w),col='red')+
  # geom_density(data=out_sample,aes(x=w,weight=weights),col='blue')+
  geom_density(data=exact,aes(x=w),col='green')
print(plot)

