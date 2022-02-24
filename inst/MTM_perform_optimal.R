library(mvtnorm)
library(coda)
library(pracma)

setwd("~/simplicialSampler/")

source("R/simplicialSampler.R")

df <- readRDS("inst/output/optimalAcceptanceMTM.rds")

set.seed(1)
#
######
############ 
#####
#
maxIt <- 10000
dimensions <- c(128,256,512) #c(4,8,16,32,64,128,256,512)
for(l in 1:3) {
  diffs <- 1/abs(dimensions[i]-df$Dimension) / sum(1/abs(dimensions[i]-df$Dimension))
  acceptParam <- sum(df$Acceptance*diffs)
  
  for(k in 1:10) {
    N <- dimensions[l]

    ptm <- proc.time()
    output1 <- MTM(N=N,x0=rep(0,N), maxIt = maxIt,adaptCov = TRUE,
                   targetName = "fullyIllGaussian",
                   adaptStepSize=TRUE,
                   targetAccept=acceptParam)
    time1 <- proc.time() - ptm
    out.mcmc1 <- as.mcmc(output1[[1]])
    eff1 <- effectiveSize(out.mcmc1)
    cat(N," ","fullyIll",
        " ",median(eff1)," ", min(eff1), " ", time1[3],"\n",
        file="inst/output/mtmUpdatedRuns.txt",
        append=TRUE) 
    
    ptm <- proc.time()
    output2 <- MTM(N=N,x0=rep(0,N), maxIt = maxIt,adaptCov = TRUE,
                   targetName = "diagGaussian",
                   adaptStepSize=TRUE,
                   targetAccept=acceptParam)
    time2 <- proc.time() - ptm
    out.mcmc2 <- as.mcmc(output2[[1]])
    eff2 <- effectiveSize(out.mcmc2)
    cat(N," ","diag",
        " ",median(eff2)," ", min(eff2), " ", time2[3],"\n",
        file="inst/output/mtmUpdatedRuns.txt",
        append=TRUE) 
    
    ptm <- proc.time()
    output3 <- MTM(N=N,x0=rep(0,N), maxIt = maxIt,adaptCov = FALSE,
                   targetName = "sphericalGaussian",
                   adaptStepSize=TRUE,
                   targetAccept=acceptParam)
    time3 <- proc.time() - ptm
    out.mcmc3 <- as.mcmc(output3[[1]])
    eff3 <- effectiveSize(out.mcmc3)
    
    cat(N," ","spherical",
        " ",median(eff3)," ", min(eff3), " ", time3[3],"\n",
        file="inst/output/mtmUpdatedRuns.txt",
        append=TRUE)   
  }
}
