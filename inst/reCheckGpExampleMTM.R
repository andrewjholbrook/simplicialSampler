library(mvtnorm)
library(coda)
library(pracma)

setwd("~/simplicialSampler/")
source("R/simplicialSampler.R")

df <- readRDS("inst/data/electionResults.rds")
y <- as.numeric(factor(df$lead)) - 1 # 1 hillary, 0 trump
x <- as.matrix(df[,2:4])
x <- as.matrix(scale(x))

set.seed(1)

nIter = 100000 # 0.3 best for vanilla, for precond

# output <- gpReg(x=x,
#                 y=y,
#                 nIter=nIter,
#                 sampler="MTM",
#                 targetAccept = 0.234,
#                 precond=TRUE)
# 
# mcmc.obj <- as.mcmc(t(output$Z[,1000:nIter]))
# effs <- mean(effectiveSize(mcmc.obj))
# effs
# 
# output <- gpReg(x=x,
#                 y=y,
#                 nIter=nIter,
#                 sampler="MTM",
#                 targetAccept = 0.3,
#                 precond=TRUE)
# 
# mcmc.obj <- as.mcmc(t(output$Z[,1000:nIter]))
# effs2 <- mean(effectiveSize(mcmc.obj))
# effs2
# 
# output <- gpReg(x=x,
#                 y=y,
#                 nIter=nIter,
#                 sampler="MTM",
#                 targetAccept = 0.4,
#                 precond=TRUE)
# 
# mcmc.obj <- as.mcmc(t(output$Z[,1000:nIter]))
# effs3 <- mean(effectiveSize(mcmc.obj))
# effs3

nIter = 100000 # simpl: 290, 1444; mtm: 124, 674; rwm: 66, 383

# No PC: 0.234 MH/MTM, 0.5 Simpl

for (k in 1:100) {
  
  ptm <- proc.time()[3]
  output <- gpReg(x=x,
                  y=y,
                  nIter=nIter,
                  sampler="MTM",
                  targetAccept = 0.3)
  time <- proc.time()[3] - ptm
  
  mcmc.obj <- as.mcmc(t(output$Z[,1000:nIter]))
  effs <- effectiveSize(mcmc.obj)
  
  preds <- round(1/(1+exp(-t(output$Z))))
  misses <- apply(X=preds,MARGIN = 1, FUN = function(x) sum(x!=y))
  i <- 0
  found <- FALSE
  while (found==FALSE & i < nIter+1) {
    i <- i + 1
    if(misses[i] == 10) {
      firstTen <- i
      found <- TRUE
    }
  }
  
  cat("vanilla ", mean(effs)," ", min(effs), " ",
      effectiveSize(as.mcmc(output$lambda[1000:nIter]))," ",
      effectiveSize(as.mcmc(output$eta[1000:nIter]))," ",
      effectiveSize(as.mcmc(output$rho[1000:nIter]))," ",
      effectiveSize(as.mcmc(output$sigma[1000:nIter]))," ",
      firstTen, " ",
      time,"\n",
      file="inst/output/gpClassificationMTM.txt",
      append=TRUE)
  
  ptm <- proc.time()[3]
  output <- gpReg(x=x,
                  y=y,
                  nIter=nIter,
                  sampler="MTM",
                  targetAccept = 0.4,
                  precond=TRUE)
  time <- proc.time()[3] - ptm
  
  mcmc.obj <- as.mcmc(t(output$Z[,1000:nIter]))
  effs <- effectiveSize(mcmc.obj)
  
  preds <- round(1/(1+exp(-t(output$Z))))
  misses <- apply(X=preds,MARGIN = 1, FUN = function(x) sum(x!=y))
  i <- 0
  found <- FALSE
  while (found==FALSE & i < nIter+1) {
    i <- i + 1
    if(misses[i] == 10) {
      firstTen <- i
      found <- TRUE
    }
  }
  
  cat("precond ", mean(effs)," ", min(effs), " ",
      effectiveSize(as.mcmc(output$lambda[1000:nIter]))," ",
      effectiveSize(as.mcmc(output$eta[1000:nIter]))," ",
      effectiveSize(as.mcmc(output$rho[1000:nIter]))," ",
      effectiveSize(as.mcmc(output$sigma[1000:nIter]))," ",
      firstTen, " ",
      time,"\n",
      file="inst/output/gpClassificationMTM.txt",
      append=TRUE)
}
