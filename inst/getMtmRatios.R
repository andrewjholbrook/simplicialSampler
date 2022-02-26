library(mvtnorm)
library(coda)

setwd("~/simplicialSampler/")

source("R/simplicialSampler.R")

maxIt <- 20000

for(N in seq())
N     <- 5

# soherical target / adapt cov
maxIt <- 10000
for(N in seq(from=250,to=300,by=5)) {
    output <- MTM(N=N, x0=rep(0,N), maxIt = maxIt,
                  adaptCov = FALSE, targetName = "sphericalGaussian")
    cat(N," ", output[[2]], "\n",
        file="inst/output/mtmAcceptanceRates.txt",
        append=TRUE)  
    
}
