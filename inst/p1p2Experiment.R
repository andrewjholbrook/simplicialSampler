library(mvtnorm)
library(coda)
library(pracma)

setwd("~/simplicialSampler/")

source("R/simplicialSampler.R")

set.seed(1)

#
######
############ comparison to P1
#####
#
maxIt <- 11000
dimensions <- 4:64
numReps <- 30
for(i in 1:length(4:64)) {
    N <- dimensions[i]
    
    firstMomentEstimators  <- matrix(0,numReps,N)
    secondMomentEstimators <- matrix(0,numReps,N)
    time1 <- 0
    meanEff1 <- 0
    for (k in 1:numReps) {
      ptm <- proc.time()
      output1 <-
        simplicialSampler(
          N = N,
          x0 = rep(0, N),
          maxIt = maxIt,
          lambda = 0.5 / sqrt(N),
          adaptStepSize = TRUE,
          targetAccept = 0.5,
          target = "sphericalGaussian",
          burnin = 1000
        )
      time1 <- time1 + proc.time() - ptm
      out.mcmc1 <- as.mcmc(output1[[1]])
      eff1 <- effectiveSize(out.mcmc1)
      meanEff1 <- meanEff1 + mean(eff1)
      
      firstMomentEstimators[k,]  <- colMeans(output1[[1]])
      secondMomentEstimators[k,] <- colMeans(output1[[1]]^2)
    }
    
    cat("simpl", " " ,N, " ",meanEff1/numReps,
        " ", time1[3]/numReps," ", mean(colMeans(firstMomentEstimators^2))," ",
        mean(colMeans((secondMomentEstimators-1)^2))," ",
        "\n",
        file="inst/output/p1P2Comparison.txt",
        append=TRUE)
    
    firstMomentEstimators  <- matrix(0,numReps,N)
    secondMomentEstimators <- matrix(0,numReps,N)
    time1 <- 0
    meanEff1 <- 0
    for (k in 1:numReps) {
      ptm <- proc.time()
      output1 <-
        tjelP1(
          N = N,
          x0 = rep(0, N),
          maxIt = maxIt,
          lambda = 0.5 / sqrt(N),
          adaptStepSize = TRUE,
          targetAccept = 0.5,
          target = "sphericalGaussian",
          nProps = N,
          burnin = 1000
        )
      time1 <- time1 + proc.time() - ptm
      out.mcmc1 <- as.mcmc(output1[[1]])
      eff1 <- effectiveSize(out.mcmc1)
      meanEff1 <- meanEff1 + mean(eff1)
      
      firstMomentEstimators[k,]  <- colMeans(output1[[1]])
      secondMomentEstimators[k,] <- colMeans(output1[[1]]^2)
    }
    
    cat("pNProps", " " ,N, " ",meanEff1/numReps,
        " ", time1[3]/numReps," ", mean(colMeans(firstMomentEstimators^2))," ",
        mean(colMeans((secondMomentEstimators-1)^2))," ",
        "\n",
        file="inst/output/p1P2Comparison.txt",
        append=TRUE)
    
    firstMomentEstimators  <- matrix(0,numReps,N)
    secondMomentEstimators <- matrix(0,numReps,N)
    time1 <- 0
    meanEff1 <- 0
    for (k in 1:numReps) {
      ptm <- proc.time()
      output1 <-
        tjelP1(
          N = N,
          x0 = rep(0, N),
          maxIt = maxIt,
          lambda = 0.5 / sqrt(N),
          adaptStepSize = TRUE,
          targetAccept = 0.5,
          target = "sphericalGaussian",
          nProps = 2*N,
          burnin = 1000
        )
      time1 <- time1 + proc.time() - ptm
      out.mcmc1 <- as.mcmc(output1[[1]])
      eff1 <- effectiveSize(out.mcmc1)
      meanEff1 <- meanEff1 + mean(eff1)
      
      firstMomentEstimators[k,]  <- colMeans(output1[[1]])
      secondMomentEstimators[k,] <- colMeans(output1[[1]]^2)
    }
    
    cat("p2NProps", " " ,N, " ",meanEff1/numReps,
        " ", time1[3]/numReps," ", mean(colMeans(firstMomentEstimators^2))," ",
        mean(colMeans((secondMomentEstimators-1)^2))," ",
        "\n",
        file="inst/output/p1P2Comparison.txt",
        append=TRUE)
    
    firstMomentEstimators  <- matrix(0,numReps,N)
    secondMomentEstimators <- matrix(0,numReps,N)
    time1 <- 0
    meanEff1 <- 0
    for (k in 1:numReps) {
      ptm <- proc.time()
      output1 <-
        tjelP1(
          N = N,
          x0 = rep(0, N),
          maxIt = maxIt,
          lambda = 0.5 / sqrt(N),
          adaptStepSize = TRUE,
          targetAccept = 0.5,
          target = "sphericalGaussian",
          nProps = 4*N,
          burnin = 1000
        )
      time1 <- time1 + proc.time() - ptm
      out.mcmc1 <- as.mcmc(output1[[1]])
      eff1 <- effectiveSize(out.mcmc1)
      meanEff1 <- meanEff1 + mean(eff1)
      
      firstMomentEstimators[k,]  <- colMeans(output1[[1]])
      secondMomentEstimators[k,] <- colMeans(output1[[1]]^2)
    }
    
    cat("p14nProps", " " ,N, " ",meanEff1/numReps,
        " ", time1[3]/numReps," ", mean(colMeans(firstMomentEstimators^2))," ",
        mean(colMeans((secondMomentEstimators-1)^2))," ",
        "\n",
        file="inst/output/p1P2Comparison.txt",
        append=TRUE)
}