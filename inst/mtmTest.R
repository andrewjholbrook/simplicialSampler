library(mvtnorm)
library(coda)

setwd("~/simplicialSampler/")

source("R/simplicialSampler.R")

maxIt <- 20000
N     <- 5

# soherical target / adapt cov
output <- MTM(N=N, x0=rep(0,N), maxIt = maxIt,
                     adaptCov = TRUE, targetName = "sphericalGaussian")

for(i in 1:N){
  qqplot(output[[1]][,i],rnorm(maxIt,sd=1))
  qqline(output[[1]][,i], distribution = function(p) qnorm(p,sd=1 ))
}


# ill conditioned target / adapt cov
output <- MTM(N=N, x0=rep(0,N), maxIt = maxIt,
                     adaptCov = TRUE, targetName = "diagGaussian")

for(i in 1:N){
  qqplot(output[[1]][,i],rnorm(maxIt,sd=sqrt(i)))
  qqline(output[[1]][,i], distribution = function(p) qnorm(p,sd=sqrt(i)))
}


# ill conditioned target / adapt cov
output <- MTM(N=N, x0=rep(0,N), maxIt = maxIt,
              adaptCov = TRUE, targetName = "diagGaussian")

for(i in 1:N){
  qqplot(output[[1]][,i],rnorm(maxIt,sd=sqrt(i)))
  qqline(output[[1]][,i], distribution = function(p) qnorm(p,sd=sqrt(i)))
}






