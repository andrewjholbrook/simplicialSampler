library(mvtnorm)
library(coda)

setwd("~/simplicialSampler/")

source("R/simplicialSampler.R")

maxIt <- 10000
N     <- 2

# spherical target / no adapt scales
output <- tjelP1(N=N, x0=rep(0,N), maxIt = maxIt, adaptStepSize=TRUE,
                 adaptScales = FALSE, target = "sphericalGaussian",
                 nProps=10)
effectiveSize(as.mcmc(output[[1]]))
for(i in 1:N){
  qqplot(output[[1]][,i],rnorm(maxIt,sd=1))
  qqline(output[[1]][,i], distribution = function(p) qnorm(p,sd=1 ))
}