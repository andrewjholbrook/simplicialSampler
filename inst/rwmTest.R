library(mvtnorm)
library(coda)

setwd("~/simplicialSampler/")

source("R/simplicialSampler.R")

maxIt <- 10000
N     <- 5

# soherical target / adapt cov
output <- randomWalk(N=N, x0=rep(0,N), maxIt = maxIt,
                  adaptCov = TRUE, target = "sphericalGaussian")

for(i in 1:N){
  qqplot(output[[1]][,i],rnorm(maxIt,sd=1))
  qqline(output[[1]][,i], distribution = function(p) qnorm(p,sd=1 ))
}


# ill conditioned target / adapt cov
output <- randomWalk(N=N, x0=rep(0,N), maxIt = maxIt,
                     adaptCov = TRUE, target = "diagGaussian")

for(i in 1:N){
  qqplot(output[[1]][,i],rnorm(maxIt,sd=sqrt(i)))
  qqline(output[[1]][,i], distribution = function(p) qnorm(p,sd=sqrt(i)))
}

# spherical target / no adapt cov
output <- randomWalk(N=N, x0=rep(0,N), maxIt = maxIt,
                     adaptCov = FALSE, target = "sphericalGaussian")

for(i in 1:N){
  qqplot(output[[1]][,i],rnorm(maxIt,sd=1))
  qqline(output[[1]][,i], distribution = function(p) qnorm(p,sd=1 ))
}


# ill conditioned target / no adapt cov
output <- randomWalk(N=N, x0=rep(0,N), maxIt = maxIt,
                     adaptCov = FALSE, target = "diagGaussian")

for(i in 1:N){
  qqplot(output[[1]][,i],rnorm(maxIt,sd=sqrt(i)))
  qqline(output[[1]][,i], distribution = function(p) qnorm(p,sd=sqrt(i)))
}

# banana target / no adapt cov
output <- randomWalk(N=N, x0=c(0,-20,rep(0,N-2)), maxIt = maxIt,
                     adaptCov = FALSE, target = "banana")
effectiveSize(as.mcmc(output[[1]]))
plot(output[[1]][,1],output[[1]][,2])
plot(output[[1]][,1],output[[1]][,3])
plot(output[[1]][,2],output[[1]][,3])
output[[3]]

# banana target / adapt cov
output <- randomWalk(N=N, x0=c(0,-20,rep(0,N-2)), maxIt = maxIt,
                     adaptCov = TRUE, target = "banana")
effectiveSize(as.mcmc(output[[1]]))
plot(output[[1]][,1],output[[1]][,2])
plot(output[[1]][,1],output[[1]][,3])
plot(output[[1]][,2],output[[1]][,3])
output[[3]]
