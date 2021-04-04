library(mvtnorm)
library(coda)
library(pracma)

setwd("~/simplicialSampler/")

source("R/simplicialSampler.R")

#
######
############ number of proposals experiment
#####
#
maxIt <- 10000
for(k in 1:100) {
  for(nProps in seq(from=5,to=100,length.out = 20)) {
    output <- simplicialSampler(N=100, x0=rep(0,100), maxIt = maxIt,
                                adaptStepSize=TRUE,
                                adaptScales = FALSE,
                                targetAccept = 0.675,
                                target = "sphericalGaussian",
                                nProps=nProps)
    out.mcmc <- as.mcmc(output[[1]])
    eff <- effectiveSize(out.mcmc)
    cat(nProps, " ",mean(eff),
        " ", min(eff) , "\n",
        file="inst/output/nPropsResults.txt",
        append=TRUE)  
  }
}