library(mvtnorm)
library(coda)
library(pracma)

#
######
############ experiment to find good lambda
#####
#
maxIt <- 100000
for(N in c(2,4,8,16,32,64,128)) {
  for(target in seq(from=0.2,to=0.6,length.out = 10)) {
    output <- inference(N=N,x0=rep(0,N), maxIt = maxIt, lambda = 3,
                        adaptStepSize = TRUE, targetAccept = target)
    out.mcmc <- as.mcmc(output[[1]])
    eff <- effectiveSize(out.mcmc)
    cat(N," ", output[[3]], " ",mean(eff),
        " ", min(eff), " ", output[[2]]," ",target,"\n",
        file="~/hedgingMCMC/output/lambdaSearch.txt",
        append=TRUE)
    traceplot(out.mcmc[,1])
  }
}


#
######
############ comparison to RWM
#####
#
maxIt <- 100000
targets <- c(0.27,0.33,0.42,0.47,0.56)
dimensions <- c(2,4,8,16,32)
for(i in 1:5) {
    N <- dimensions[i]
    ptm <- proc.time()
    output1 <- inference(N=N,x0=rep(0,N), maxIt = maxIt, lambda = 3,
                        adaptStepSize = TRUE, targetAccept = targets[i])
    time1 <- proc.time() - ptm
    out.mcmc1 <- as.mcmc(output1[[1]])
    eff1 <- effectiveSize(out.mcmc1)
    ptm <- proc.time()
    output2 <- randomWalk(N=N,x0=rep(0,N), maxIt = maxIt,
                         adaptStepSize = TRUE, targetAccept = 0.25)
    time2 <- proc.time() - ptm
    out.mcmc2 <- as.mcmc(output2[[1]])
    eff2 <- effectiveSize(out.mcmc2)
    
    cat(N, " ",mean(eff1)," ", min(eff1), " ", time1[3]," ",mean(eff2)," ", min(eff2), " ", time2[3],"\n",
        file="~/hedgingMCMC/output/rwmComparison.txt",
        append=TRUE)
    #traceplot(out.mcmc[,1])
}

#
######
############ time stewart vs rando
#####
#
N <- 100
U <- stewart(N)
D <- U%*%t(U) - diag(N)
norm(D)

system.time(stewart(N))
system.time(randortho(n=N))

