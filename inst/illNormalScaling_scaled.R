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
for(N in seq(from=100,to=500,by=10)) {
  for(targetAccept in seq(from=0.2,to=0.95,length.out = 20)) {
    output <- simplicialSampler(N=N, x0=rep(0,N), maxIt = maxIt,
                                adaptStepSize=TRUE,
                                adaptScales = TRUE,
                                targetAccept = targetAccept,
                                target = "diagGaussian")
    out.mcmc <- as.mcmc(output[[1]])
    eff <- effectiveSize(out.mcmc)*10 # hack to compare to old ESSs (out of 100,000)
    cat(N," ", output[[3]], " ",median(eff),
        " ", min(eff), " ", output[[2]]," ",targetAccept,
        " ", "diagGaussian",  " ", "scaled", "\n",
        file="inst/output/lambdaSearch.txt",
        append=TRUE)  }
}