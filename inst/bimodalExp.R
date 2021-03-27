library(mvtnorm)
library(coda)
library(pracma)

setwd("~/simplicialSampler/")

source("R/simplicialSampler.R")

dimensions <- 2:10
maxIt     <- 10000

numberJumps <- function(chain) {
  N <- dim(chain)[1]
  count <- 0
  
  for (i in 2:N) {
    if ( (sum(chain[i-1,]) < 5 & sum(chain[i,]) > 5) |
         (sum(chain[i-1,]) > 5 & sum(chain[i,]) < 5) ) {
      count <- count + 1
    }
  }
  return(count)
}

set.seed(666)

for(i in 1:9) {
  for(k in 1:100){
    N <- dimensions[i]
    
    ptm <- proc.time()
    output2 <- randomWalk(N=N,x0=rep(0,N), maxIt = maxIt,
                           target = "bimodalGaussian", adaptCov = TRUE)
    time2 <- proc.time() - ptm
    nJumps2 <- numberJumps(output2[[1]])
    
    ptm <- proc.time()
    output1 <- simplicialSampler(N=N,x0=rep(0,N), maxIt = maxIt,lambda=2.4/sqrt(N),
                                 adaptStepSize = TRUE, targetAccept = output2[[2]],
                                 target = "bimodalGaussian", adaptScales = TRUE,
                                 Gaussians = TRUE)
    time1 <- proc.time() - ptm
    nJumps1 <- numberJumps(output1[[1]])
    

    ptm <- proc.time()
    output3 <- randomWalk(N=N,x0=rep(0,N), maxIt = maxIt,
                          target = "bimodalGaussian", adaptCov = FALSE)
    time3 <- proc.time() - ptm
    nJumps3 <- numberJumps(output3[[1]])
    
    
    ptm <- proc.time()
    output4 <- simplicialSampler(N=N,x0=rep(0,N), maxIt = maxIt,lambda=2.4/sqrt(N),
                                 adaptStepSize = TRUE, targetAccept = output3[[2]],
                                 target = "bimodalGaussian", adaptScales = FALSE,
                                 Gaussians = TRUE)
    time4 <- proc.time() - ptm
    nJumps4 <- numberJumps(output4[[1]])
    
    cat(N, " ",nJumps1, " ", time1[3],
        " ", nJumps2, " ", time2[3],
        " ", nJumps3, " ", time3[3],
        " ", nJumps4, " ", time4[3], "\n",
        file="inst/output/bimodalComparison.txt",
        append=TRUE)
  }
}







