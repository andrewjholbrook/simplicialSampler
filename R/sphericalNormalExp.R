library(rstiefel)
library(mvtnorm)
library(coda)
library(pracma)
library(Matrix)

f <- function(X) {
  densities <- dmvnorm(X)
  return(densities)
}

householder <- function(x,r=NULL) {
  # takes vector x and returns householder matrix
  N <- length(x)
  if(is.null(r)) r <- - sqrt(sum(x^2)) * sign(x[N])
  u <- x - r*c(1,rep(0,N-1))
  v <- u / sqrt(sum(u^2))
  H <- diag(N) - 2*v%*%t(v)
}

stewart <- function(N) {
  # the stewart 1980 algorithm for generating a haar distributed ortho matrix
  xs <- list()
  rs <- rep(0,N)
  Hs <- list()
  U  <- diag(N)
  for(i in (N-1):1) {
    x       <- rnorm(i+1)
    xs[[i]] <- x
    rs[i]   <- - sqrt(sum(x^2)) * sign(x[i+1])
    U <- U %*% bdiag(diag(N-1-i), householder(x,rs[i]))
  }
  rs[N] <- 1
  U <- diag(sign(rs)) %*% U
  return(U)
}

proposal <- function(N,x,lambda) {
  # N is dimensionality
  # x is current state
  # lambda is scale
  v_Nplus1 <- rep((1+sqrt(N+1))/N,N)
  M <- rbind(diag(N), v_Nplus1)
  M <- M - matrix(v_Nplus1,N+1,N)
  M4 <- M * lambda /sqrt(2)
  U <- randortho(n=N)
  M4 <- M4 %*% U
  M4 <- M4 + matrix(x,N+1,N,byrow = TRUE)
  
  output <- M4[sample(1:(N+1),1,prob = f(M4)),]
  
  return(output)
}

inference <- function(N,x0,lambda=1,maxIt=10000,adaptStepSize=TRUE,
                      targetAccept=0.5) {
  if(N!=length(x0)) stop("Dimension mismatch.")
  
  chain <- matrix(0,maxIt,N)
  
  Acceptances = 0 # total acceptances within adaptation run (<= SampBound)
  SampBound = 5   # current total samples before adapting radius
  SampCount = 0   # number of samples collected (adapt when = SampBound)
  Proposed = 0;
  
  accept <- rep(0,maxIt)
  chain[1,] <- x0
  for (i in 2:maxIt){
    Proposed = Proposed + 1
    chain[i,] <- proposal(N,chain[i-1,],lambda)
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

randomWalk <- function(N,x0,maxIt=10000,adaptStepSize=TRUE,
                      targetAccept=0.25) {
  if(N!=length(x0)) stop("Dimension mismatch.")
  
  chain <- matrix(0,maxIt,N)
  
  sigma <- 1/sqrt(N)
  
  Acceptances = 0 # total acceptances within adaptation run (<= SampBound)
  SampBound = 5   # current total samples before adapting radius
  SampCount = 0   # number of samples collected (adapt when = SampBound)
  Proposed = 0;
  
  accept <- rep(0,maxIt)
  chain[1,] <- x0
  for (i in 2:maxIt){
    Proposed = Proposed + 1
    xStar <- rnorm(N,sd=sigma) + chain[i-1,]
    if(log(runif(1)) < sum(dnorm(xStar,log=TRUE)) -
       sum(dnorm(chain[i-1,],log=TRUE))){
      accept[i] <- 1
      Acceptances = Acceptances + 1
      chain[i,] <- xStar
    } else {
      chain[i,] <- chain[i-1,]
    }
    SampCount <- SampCount + 1
    
    
    if (SampCount == SampBound) { # adjust lambda at increasing intervals
      
      AcceptRatio <- Acceptances / SampBound
      AdaptRatio  <- AcceptRatio / targetAccept
      if (AdaptRatio>2) AdaptRatio <- 2
      if (AdaptRatio<0.5) AdaptRatio <- 0.5
      sigma <- sigma * AdaptRatio
      
      SampCount <- 0
      SampBound <- ceiling(SampBound^1.01)
      Acceptances <- 0
    }
    
    if(i %% 1000 == 0) cat(i,"\n")
  }
  ratio <- sum(accept)/(maxIt-1)
  cat("Acceptance ratio: ", ratio,"\n")
  return(list(chain,ratio,sigma))
}



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

