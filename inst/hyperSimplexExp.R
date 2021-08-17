library(mvtnorm)
library(coda)

setwd("~/simplicialSampler/")

source("R/simplicialSampler.R")

maxIt <- 10000
N     <- 3

proposal <- function(N,x,lambda,distrib,gaussians=FALSE) {
    
    Dim <- length(x)
    x <- c(x,rep(0,N-Dim))
    
    v_Nplus1 <- rep((1+sqrt(N+1))/N,N)
    M <- rbind(diag(N), v_Nplus1)
    M <- M - matrix(v_Nplus1,N+1,N)
    M4 <- M * lambda /sqrt(2) 
    if(gaussians) M4 <- M4 * sqrt(rchisq(n=1,df=N))
    U <- pracma::randortho(n=N)
    M4 <- M4 %*% U
    M4 <- M4 + matrix(x,N+1,N,byrow = TRUE)
    M4 <- M4[,1:Dim]
    output <- M4[sample(1:(N+1),1,prob = target(M4,distrib)),]
    
    return(output)
}

simplicialSampler <- function(nProps,N,x0,lambda=1,maxIt=10000,adaptStepSize=TRUE,
                              targetAccept=0.5,target=NULL,
                              Gaussians=FALSE) {
    #if(N!=length(x0)) stop("Dimension mismatch.")
    
    chain <- matrix(0,maxIt,N)
    
    Acceptances = 0 # total acceptances within adaptation run (<= SampBound)
    SampBound = 5   # current total samples before adapting radius
    SampCount = 0   # number of samples collected (adapt when = SampBound)
    Proposed = 0
    
    accept <- rep(0,maxIt)
    chain[1,] <- x0
    for (i in 2:maxIt){
        Proposed = Proposed + 1
        
        
        prpsl <- proposal(nProps,chain[i-1,],lambda,target,
                          gaussians = Gaussians)
        
        chain[i,] <- prpsl
        SampCount <- SampCount + 1
        if(any(chain[i,] != chain[i-1,])){
            accept[i] <- 1
            Acceptances = Acceptances + 1
        }  
        
        if (SampCount == SampBound) { # adjust lambda at increasing intervals
            AcceptRatio <- Acceptances / SampBound
            AdaptRatio  <- AcceptRatio / targetAccept
            if (AdaptRatio>2) AdaptRatio <- 2
            if (AdaptRatio<0.5) AdaptRatio <- 0.5
            lambda <- lambda * AdaptRatio
            
            SampCount <- 0
            SampBound <- ceiling(SampBound^1.01)
            Acceptances <- 0
        }
        
        
        if(i %% 1000 == 0) cat(i,"\n")
    }
    ratio <- sum(accept)/(maxIt-1)
    cat("Acceptance ratio: ", ratio,"\n")
    return(list(chain,ratio,lambda))
}

# spherical target / no adapt scales
output <- simplicialSampler(nProps=100,N=N, x0=rep(0,N), maxIt = 100000, adaptStepSize=TRUE,
                            target = "sphericalGaussian",targetAccept = 0.5)
effectiveSize(as.mcmc(output[[1]]))
for(i in 1:N){
    qqplot(output[[1]][,i],rnorm(maxIt,sd=1))
    qqline(output[[1]][,i], distribution = function(p) qnorm(p,sd=1 ))
}


#
### make figure
#

prpsl <- proposal(nProps,chain[i-1,],lambda,target,
                  gaussians = Gaussians)

