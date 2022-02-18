library(mvtnorm)
library(coda)
library(pracma)

setwd("~/simplicialSampler/")

source("R/simplicialSampler.R")

#
######
############ comparison to P1
#####
#
maxIt <- 20000
dimensions <- 4:64
for(i in 1:length(4:64)) {
  for(k in 1:10) {
    N <- dimensions[i]
    ptm <- proc.time()
    output1 <- simplicialSampler(N=N,x0=rep(0,N), maxIt = maxIt, lambda = 0.5/sqrt(N),
                                 adaptStepSize = TRUE, targetAccept = 0.5,
                                 target = "sphericalGaussian")
    time1 <- proc.time() - ptm
    out.mcmc1 <- as.mcmc(output1[[1]])
    eff1 <- effectiveSize(out.mcmc1)
    
    ptm <- proc.time()
    output2 <- tjelP1(N=N, x0=rep(0,N), maxIt = maxIt, adaptStepSize=TRUE,
                     adaptScales = FALSE, target = "sphericalGaussian",
                     lambda = 0.5/sqrt(N),
                     nProps=N)
    time2 <- proc.time() - ptm
    out.mcmc2 <- as.mcmc(output2[[1]])
    eff2 <- effectiveSize(out.mcmc2)
    
    ptm <- proc.time()
    output3 <- tjelP1(N=N, x0=rep(0,N), maxIt = maxIt, adaptStepSize=TRUE,
                         adaptScales = FALSE, target = "sphericalGaussian",
                         lambda = 0.5/sqrt(N),
                         nProps=2*N)
    time3 <- proc.time() - ptm
    out.mcmc3 <- as.mcmc(output3[[1]])
    eff3 <- effectiveSize(out.mcmc3)
    
    ptm <- proc.time()
    output4 <- tjelP1(N=N, x0=rep(0,N), maxIt = maxIt, adaptStepSize=TRUE,
                      adaptScales = FALSE, target = "sphericalGaussian",
                      lambda = 0.5/sqrt(N),
                      nProps=4*N)
    time4 <- proc.time() - ptm
    out.mcmc4 <- as.mcmc(output4[[1]])
    eff4 <- effectiveSize(out.mcmc4)
    
    cat(N, " ",mean(eff1), " ", time1[3],
        " ",mean(eff2), " ", time2[3],
        " ",mean(eff3)," ", time3[3],
        " ",mean(eff4), " ", time4[3],"\n",
        file="inst/output/p1P2Comparison.txt",
        append=TRUE)   
  }
}