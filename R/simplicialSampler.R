
target <- function(X,distrib=NULL) {
  if(is.null(distrib)) stop("Target distribution must be specified.")
  if (distrib=="sphericalGaussian") {
    densities <- mvtnorm::dmvnorm(X)
  } else if (distrib=="diagGaussian") {
    densities <- mvtnorm::dmvnorm(X,sigma = diag(1:dim(X)[2]))
  } else if (distrib=="banana") {
    densities <- banana(X,B=0.1)
  } else {
    stop("Specified distribution not supported.")
  }
  return(densities)
}

banana <- function(X,B) {
  # B is bananicity constant
  N <- dim(X)[2]
  eponent <- - X[,1]^2/200 - 0.5*(X[,2]+B*X[,1]^2-100*B)^2 -
    0.5*rowSums(X[,3:N]^2)
  return(exp(exponent))
}

proposal <- function(N,x,lambda) {
  # N is dimensionality
  # x is current state
  # lambda is scale
  v_Nplus1 <- rep((1+sqrt(N+1))/N,N)
  M <- rbind(diag(N), v_Nplus1)
  M <- M - matrix(v_Nplus1,N+1,N)
  M4 <- M * lambda /sqrt(2)
  U <- pracma::randortho(n=N)
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