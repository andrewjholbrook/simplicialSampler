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
maxIt <- 10000
for(N in seq(from=5,to=500,by=5)) {
  for(targetAccept in seq(from=0.2,to=0.95,length.out = 20)) {
    output <- MTM(N=N, x0=rep(0,N), maxIt = maxIt, adaptStepSize=TRUE,
                  targetAccept = targetAccept,
                  adaptCov = FALSE, targetName = "sphericalGaussian")
    out.mcmc <- as.mcmc(output[[1]])
    eff <- effectiveSize(out.mcmc)
    cat(N," ", output[[3]], " ",mean(eff), " ", output[[2]]," ",targetAccept, "\n",
        file="inst/output/mtmSigmaSearch.txt",
        append=TRUE)  }
}