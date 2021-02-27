library(mvtnorm)
library(coda)

setwd("~/simplicialSampler/")

source("R/simplicialSampler.R")

maxIt <- 100000
N     <- 5

# spherical target / no adapt scales
output <- simplicialSampler(N=N, x0=rep(0,N), maxIt = maxIt, adaptStepSize=TRUE,
                     adaptScales = FALSE, target = "sphericalGaussian")
effectiveSize(as.mcmc(output[[1]]))
for(i in 1:N){
  qqplot(output[[1]][,i],rnorm(maxIt,sd=1))
  qqline(output[[1]][,i], distribution = function(p) qnorm(p,sd=1 ))
}


# ill conditioned target / no adapt scales
output <- simplicialSampler(N=N, x0=rep(0,N), maxIt = maxIt, adaptStepSize=TRUE,
                            adaptScales = FALSE, target = "diagGaussian")
effectiveSize(as.mcmc(output[[1]]))
for(i in 1:N){
  qqplot(output[[1]][,i],rnorm(maxIt,sd=sqrt(i)))
  qqline(output[[1]][,i], distribution = function(p) qnorm(p,sd=sqrt(i)))
}


# spherical target / adapt scales
output <- simplicialSampler(N=N, x0=rep(0,N), maxIt = maxIt, adaptStepSize=TRUE,
                            adaptScales = TRUE, target = "sphericalGaussian")
effectiveSize(as.mcmc(output[[1]]))
for(i in 1:N){
  qqplot(output[[1]][,i],rnorm(maxIt,sd=1))
  qqline(output[[1]][,i], distribution = function(p) qnorm(p,sd=1 ))
}


# ill conditioned target / adapt scales
output <- simplicialSampler(N=N, x0=rep(0,N), maxIt = maxIt, adaptStepSize=TRUE,
                            adaptScales = TRUE, target = "diagGaussian")
effectiveSize(as.mcmc(output[[1]]))
for(i in 1:N){
  qqplot(output[[1]][,i],rnorm(maxIt,sd=sqrt(i)))
  qqline(output[[1]][,i], distribution = function(p) qnorm(p,sd=sqrt(i)))
}

# banana target / no adapt scales
output <- simplicialSampler(N=N, x0=c(0,-20,rep(0,N-2)), maxIt = maxIt, adaptStepSize=TRUE,
                            adaptScales = FALSE, target = "banana")
effectiveSize(as.mcmc(output[[1]]))
plot(output[[1]][,1],output[[1]][,2])
plot(output[[1]][,1],output[[1]][,3])
plot(output[[1]][,2],output[[1]][,3])
output[[3]]

# banana target / adapt scales
output <- simplicialSampler(N=N, x0=c(0,-20,rep(0,N-2)), maxIt = maxIt, adaptStepSize=TRUE,
                            adaptScales = TRUE, target = "banana")
effectiveSize(as.mcmc(output[[1]]))
plot(output[[1]][,1],output[[1]][,2])
plot(output[[1]][,1],output[[1]][,3])
plot(output[[1]][,2],output[[1]][,3])
output[[3]]

