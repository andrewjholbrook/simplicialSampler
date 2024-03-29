library(mvtnorm)
library(coda)
library(pracma)

setwd("~/simplicialSampler/")

source("R/simplicialSampler.R")

#
######
############ comparison to RWM
#####
#
maxIt <- 10000
targets <- c(0.33,0.4,0.47,0.52,0.57,0.675,0.675,0.675)
dimensions <- c(4,8,16,32,64,128,256,512)
for(i in 1:8) {
    for(k in 1:10) {
        N <- dimensions[i]
        ptm <- proc.time()
        output1 <- simplicialSampler(N=N,x0=rep(0,N), maxIt = maxIt, lambda = 6,
                                     adaptStepSize = TRUE, targetAccept = targets[i],
                                     target = "sphericalGaussian")
        time1 <- proc.time() - ptm
        out.mcmc1 <- as.mcmc(output1[[1]])
        eff1 <- effectiveSize(out.mcmc1)
        ptm <- proc.time()
        output2 <- randomWalk(N=N,x0=rep(0,N), maxIt = maxIt,adaptCov = FALSE,
                              target = "sphericalGaussian")
        time2 <- proc.time() - ptm
        out.mcmc2 <- as.mcmc(output2[[1]])
        eff2 <- effectiveSize(out.mcmc2)
        
        ptm <- proc.time()
        output3 <- MTM(N=N,x0=rep(0,N), maxIt = maxIt,adaptCov = FALSE,
                              target = "sphericalGaussian")
        time3 <- proc.time() - ptm
        out.mcmc3 <- as.mcmc(output3[[1]])
        eff3 <- effectiveSize(out.mcmc3)
        
        cat(N, " ",median(eff1)," ", min(eff1), " ", time1[3],
            " ",median(eff2)," ", min(eff2), " ", time2[3],
            " ",median(eff3)," ", min(eff3), " ", time3[3],"\n",
            file="inst/output/rwmComparison.txt",
            append=TRUE)   
    }
}

