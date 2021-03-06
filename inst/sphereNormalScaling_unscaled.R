library(mvtnorm)
library(coda)
library(pracma)

setwd("~/simplicialSampler/")

source("R/simplicialSampler.R")

#
######
############ experiment to find good lambda
#####
#
maxIt <- 100000
for(N in seq(from=60,to=100,by=5)) {
  for(targetAccept in seq(from=0.2,to=0.95,length.out = 20)) {
    output <- simplicialSampler(N=N, x0=rep(0,N), maxIt = maxIt,
                                adaptStepSize=TRUE,
                                adaptScales = FALSE,
                                targetAccept = targetAccept,
                                target = "sphericalGaussian")
    out.mcmc <- as.mcmc(output[[1]])
    eff <- effectiveSize(out.mcmc)
    cat(N," ", output[[3]], " ",median(eff),
        " ", min(eff), " ", output[[2]]," ",targetAccept,
        " ", "sphericalGaussian",  " ", "unscaled", "\n",
        file="inst/output/lambdaSearch.txt",
        append=TRUE)  }
}