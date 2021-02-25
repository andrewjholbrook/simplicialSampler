setwd("~/simplicialSampler/")

source("R/simplicialSampler.R")

library(coda)

n <- 6
results <- simplicialSampler(N=n,x0=rep(0,n), maxIt = 100000,lambda=2,
                     adaptStepSize = TRUE,targetAccept = 0.25,
                     target = "banana")
out.mcmc <- as.mcmc(results[[1]])
effectiveSize(out.mcmc)
plot(out.mcmc)
plot(results[[1]][,1],results[[1]][,2])
plot(results[[1]][,1],results[[1]][,3])
plot(results[[1]][,2],results[[1]][,3])
results[[3]] # lambda

# MH
results2 <- randomWalk(N=n,x0=rep(0,n), maxIt = 100000,
                     adaptStepSize = TRUE,targetAccept = 0.25,
                     target = "banana")
out.mcmc2 <- as.mcmc(results2[[1]])
effectiveSize(out.mcmc2)
plot(out.mcmc2)
plot(results2[[1]][,1],results2[[1]][,2])
plot(results2[[1]][,1],results2[[1]][,3])
plot(results2[[1]][,2],results2[[1]][,3])
results2[[3]]



