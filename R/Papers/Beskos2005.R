# Simulate start point
starting_sample_size = 1000000
end_h = rnorm(starting_sample_size)
end_accept = ifelse(runif(starting_sample_size) < exp(-cos(end_h)-1),TRUE, FALSE)
end_h = end_h[end_accept]
sample_size = sum(end_accept)

# Simulate poisson process
fun <- function(x,high,ordered=FALSE) {
  if (ordered) {
    sort(runif(x,0,high))
  } else {
    runif(x,0,high)
  }
  
}


po_rate = 9.0/8.0
t = 1
pp = data.frame(k = rpois(sample_size,po_rate*t),end_h=end_h)
pp$points = apply(pp['k'],1,fun,high=t,ordered=TRUE)
pp$marks = apply(pp['k'],1,fun,high=1)

# Proposal Distribution
bb <- function(points,end) {
  x <- c(0,unlist(points,use.names=FALSE),t)
  vars <- x[2:length(x)] - x[1:length(x)-1]
  b <- cumsum(rnorm(length(x)-1,0,vars))
  b <- b - x[2:length(x)]/t * tail(b,n=1)
  b + x[2:length(x)]/t*end
}

pp$bb = mapply(bb,points=pp$points,end=pp$end_h)

# Acceptance
accept <- function(bb,marks,po_rate){
  num_points <- length(bb)-1
  all(marks < 1/po_rate * (sin(bb[1:num_points]-pi)**2 + cos(bb[1:num_points]-pi))/2)
}

pp$accept = mapply(accept,bb=pp$bb,marks=pp$marks,po_rate=po_rate)

end_sample_size <- sum(pp$accept)

