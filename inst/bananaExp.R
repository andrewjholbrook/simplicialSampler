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
set.seed(666)
N <- 20
maxIt <- 10000
for(k in 1:100) {
  ptm <- proc.time()
  output1 <- simplicialSampler(N=N,x0=rep(0,N), maxIt = maxIt, lambda = 3,
                               adaptStepSize = TRUE, targetAccept = 0.25,
                               target = "banana")
  time1 <- proc.time() - ptm
  out.mcmc1 <- as.mcmc(output1[[1]])
  eff1 <- effectiveSize(out.mcmc1)

  ptm <- proc.time()
  output2 <- simplicialSampler(N=N,x0=rep(0,N), maxIt = maxIt, lambda = 3,
                               adaptStepSize = TRUE, targetAccept = 0.25,
                               target = "banana", adaptScales = TRUE)
  time2 <- proc.time() - ptm
  out.mcmc2 <- as.mcmc(output2[[1]])
  eff2 <- effectiveSize(out.mcmc2)

  ptm <- proc.time()
  output3 <- randomWalk(N=N,x0=rep(0,N), maxIt = maxIt,adaptCov = FALSE,
                        target = "banana")
  time3 <- proc.time() - ptm
  out.mcmc3 <- as.mcmc(output3[[1]])
  eff3 <- effectiveSize(out.mcmc3)
  
  ptm <- proc.time()
  output4 <- randomWalk(N=N,x0=rep(0,N), maxIt = maxIt,adaptCov = TRUE,
                        target = "banana")
  time4 <- proc.time() - ptm
  out.mcmc4 <- as.mcmc(output4[[1]])
  eff4 <- effectiveSize(out.mcmc4)
  
  cat(N, " ",mean(eff1)," ", min(eff1), " ", time1[3],
      " ",mean(eff2)," ", min(eff2), " ", time2[3],
      " ",mean(eff3)," ", min(eff3), " ", time3[3],
      " ",mean(eff4)," ", min(eff4), " ", time4[3], "\n",
      file="inst/output/bananaComparisons.txt",
      append=TRUE)   
}

