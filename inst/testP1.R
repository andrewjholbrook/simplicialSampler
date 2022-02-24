library(mvtnorm)
library(coda)

setwd("~/simplicialSampler/")

source("R/simplicialSampler.R")

maxIt <- 10000
N     <- 10

# spherical target / no adapt scales
ptm <- proc.time()[3]
output <- tjelP1(N=N, x0=rep(0,N), maxIt = maxIt, adaptStepSize=TRUE,
                 adaptScales = FALSE, target = "sphericalGaussian",
                 nProps=N)
effectiveSize(as.mcmc(output[[1]]))
time1 <- proc.time()[3] - ptm
for(i in 1:N){
  qqplot(output[[1]][,i],rnorm(maxIt,sd=1))
  qqline(output[[1]][,i], distribution = function(p) qnorm(p,sd=1 ))
}


# slow
ptm <- proc.time()[3]
output <- tjelP1(N=N, x0=rep(0,N), maxIt = maxIt, adaptStepSize=TRUE,
                 adaptScales = FALSE, target = "sphericalGaussian",
                 nProps=N, slow=TRUE)
effectiveSize(as.mcmc(output[[1]]))
time2 <- proc.time()[3] - ptm
for(i in 1:N){
  qqplot(output[[1]][,i],rnorm(maxIt,sd=1))
  qqline(output[[1]][,i], distribution = function(p) qnorm(p,sd=1 ))
}