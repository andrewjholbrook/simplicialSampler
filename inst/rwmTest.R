library(mvtnorm)
library(coda)

setwd("~/simplicialSampler/")

source("R/simplicialSampler.R")

maxIt <- 100000
N     <- 5

# soherical target
output <- randomWalk(N=N, x0=rep(0,N), maxIt = maxIt,
                  adaptCov = TRUE, target = "sphericalGaussian")

for(i in 1:N){
  qqplot(output[[1]][,i],rnorm(maxIt,sd=1))
  qqline(output[[1]][,i], distribution = function(p) qnorm(p,sd=1 ))
}


# ill conditioned target
output <- randomWalk(N=N, x0=rep(0,N), maxIt = maxIt,
                     adaptCov = TRUE, target = "diagGaussian")

for(i in 1:N){
  qqplot(output[[1]][,i],rnorm(maxIt,sd=sqrt(i)))
  qqline(output[[1]][,i], distribution = function(p) qnorm(p,sd=sqrt(i)))
}

