sample <- data.frame(proposal = runif(100000))
sample$targetDensity <- dbeta(sample$proposal, 3, 6)
maxDens = max(sample$targetDensity, na.rm = T)
sample$accepted = ifelse(runif(100000) < sample$targetDensity / maxDens, TRUE, FALSE)
hist(sample$proposal[sample$accepted], freq=F, col = "grey", breaks = 100)
curve(dbeta(x,3,6),0,1,add=TRUE,col='red')


