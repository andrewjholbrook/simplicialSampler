
recursion <- function(Ct,XbarMinus,Xt,epsilon,sd,t,warmup=100) {
  if(t>warmup) {
    XbarT  <- Xt/t + (t-1)/t*XbarMinus
    CtPlus <- Ct*(t-1)/t + sd/t*( t*XbarMinus%*%t(XbarMinus) -
                                    (t+1)*XbarT%*%t(XbarT) +
                                    Xt%*%t(Xt) + epsilon*diag(length(Xt)))
  } else {
    XbarT  <- Xt/t + (t-1)/t*XbarMinus
    CtPlus <- Ct
  }
  return(list(CtPlus,XbarT))
}

logsumexp <- function(x) {
  c <- max(x)
  return(c+log(sum(exp(x-c))))
}

target <- function(X,distrib=NULL,parallel=FALSE) {
  if(is.null(distrib)) stop("Target distribution must be specified.")
  if (distrib=="sphericalGaussian") {
    if (is.vector(X)) { # RWM
      if (parallel) {
        densities <- mvnfast::dmvn(X,sigma = diag(rep(1,length(X))),
                                   mu= rep(0,length(X)),
                                   log = TRUE, ncores=1)
      } else {
        densities <- mvtnorm::dmvnorm(X,sigma = diag(rep(1,length(X))),
                                      log = TRUE)
      }
    } else if (is.matrix(X)) { # SS
      if (parallel) {
        densities <- mvnfast::dmvn(X,sigma = diag(rep(1,dim(X)[2])),
                                   mu=rep(0,dim(X)[2]),
                                   log = TRUE, ncores=10)
      } else {
        densities <- mvtnorm::dmvnorm(X,sigma = diag(rep(1,dim(X)[2])),log = TRUE)
      }
      densities <- exp(densities-logsumexp(densities)) # helps with underflow
    } else {
      stop("States must be vectors or matrices.")
    }
  } else if (distrib=="diagGaussian") {
    if (is.vector(X)) { # RWM
      if (parallel) {
        densities <- mvnfast::dmvn(X,sigma = diag(1:length(X)),
                                   mu= rep(0,length(X)),
                                   log = TRUE, ncores=1)
      } else {
        densities <- mvtnorm::dmvnorm(X,sigma = diag(1:length(X)), log = TRUE)
      }
    } else if (is.matrix(X)) { # SS
      if (parallel) {
        densities <- mvnfast::dmvn(X,sigma = diag(1:dim(X)[2]),
                                   mu=rep(0,dim(X)[2]),
                                   log = TRUE, ncores=10)
      } else {
        densities <- mvtnorm::dmvnorm(X,sigma = diag(1:dim(X)[2]),log = TRUE)
      }
      densities <- exp(densities-logsumexp(densities)) # helps with underflow
    } else {
      stop("States must be vectors or matrices.")
    }
  } else if (distrib=="banana") {
    densities <- banana(X,B=0.1) 
  } else if (distrib=="bimodalGaussian") {
    if (is.vector(X)) { 
      densities <- log( 0.5*(mvtnorm::dmvnorm(X,mean=rep(5,length(X))) +
                               mvtnorm::dmvnorm(X)) )
    } else if (is.matrix(X)) { 
      densities <- log( 0.5*(mvtnorm::dmvnorm(X,mean=rep(5,ncol(X))) +
                               mvtnorm::dmvnorm(X)) )
      densities <- exp(densities-logsumexp(densities)) # helps with underflow
    } else {
      stop("States must be vectors or matrices.")
    }
  } else if (distrib == "fullyIllGaussian") {
    if (is.vector(X)) { # RWM
      C <- readRDS(file = paste0("inst/savedCovs/cov",length(X),".rds"))
      if (parallel) {
        densities <- mvnfast::dmvn(X,sigma = C,
                                   mu= rep(0,length(X)),
                                   log = TRUE, ncores=1)
      } else {
        densities <- mvtnorm::dmvnorm(X,sigma = C, log = TRUE)
      }
    } else if (is.matrix(X)) { # SS
      C <- readRDS(file = paste0("inst/savedCovs/cov",dim(X)[2],".rds"))
      if (parallel) {
        densities <- mvnfast::dmvn(X,sigma = C,
                                   mu=rep(0,dim(X)[2]),
                                   log = TRUE, ncores=10)
      } else {
        densities <- mvtnorm::dmvnorm(X,sigma = C,log = TRUE)
      }
      densities <- exp(densities-logsumexp(densities)) # helps with underflow
    } else {
      stop("States must be vectors or matrices.")
    }
  } else {
    stop("Specified distribution not supported.")
  }
  return(densities)
}

banana <- function(X,B) { 
  # B is bananicity constant
  if (is.vector(X)) {
    N <- length(X)
    output <- - X[1]^2/200 - 0.5*(X[2]+B*X[1]^2-100*B)^2 -
      0.5*sum(X[3:N]^2)
  } else if (is.matrix(X)) {
    N <- dim(X)[2]
    output <- exp( - X[,1]^2/200 - 0.5*(X[,2]+B*X[,1]^2-100*B)^2 -
                     0.5*rowSums(X[,3:N]^2) )
  } else {
    stop("States must be vectors or matrices.")
  }
  return(output)
}

hawkes_posterior <- function(X,engine) {
  if(is.vector(X)) {
    output <- - hpHawkes::Potential(engine,exp(X))
  } else {
    stop("Input must be vector.")
  }
  return(output)
}


proposal <- function(N,x,lambda,distrib,adaptScales=FALSE,Ct=NULL) {
  # N is dimensionality
  # x is current state
  # lambda is scale
  # distrib is target distribution
  v_Nplus1 <- rep((1+sqrt(N+1))/N,N)
  M <- rbind(diag(N), v_Nplus1)
  M <- M - matrix(v_Nplus1,N+1,N)
  M4 <- M * lambda /sqrt(2)
  U <- pracma::randortho(n=N)
  if (adaptScales) {
    M4 <- M4 %*% U %*% chol(Ct)
  } else {
    M4 <- M4 %*% U
  }
  M4 <- M4 + matrix(x,N+1,N,byrow = TRUE)
  output <- M4[sample(1:(N+1),1,prob = target(M4,distrib)),]
  return(output)
}

simplicialSampler <- function(N,x0,lambda=1,maxIt=10000,adaptStepSize=TRUE,
                              targetAccept=0.5,target=NULL,adaptScales=FALSE) {
  if(N!=length(x0)) stop("Dimension mismatch.")
  
  chain <- matrix(0,maxIt,N)
  
  Acceptances = 0 # total acceptances within adaptation run (<= SampBound)
  SampBound = 5   # current total samples before adapting radius
  SampCount = 0   # number of samples collected (adapt when = SampBound)
  Proposed = 0
  
  if (adaptScales) {
    Ct <- diag(N)
    xbar <- x0 
  } else {
    Ct <- NULL
  }
  
  accept <- rep(0,maxIt)
  chain[1,] <- x0
  for (i in 2:maxIt){
    Proposed = Proposed + 1
    
    prpsl <- NULL
    attempt <- 0
    while( is.null(prpsl) && attempt <= 100 ) {
      attempt <- attempt + 1
      try(
        prpsl <- proposal(N,chain[i-1,],lambda,target,
                          adaptScales = adaptScales, Ct=Ct)
      )
    } 
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
    
    if (adaptScales) {
      updt <- recursion(Ct=Ct,     # multivariate scales update
                        XbarMinus=xbar,
                        Xt=chain[i,],
                        epsilon = 0.000001,
                        sd=1,
                        t=i,
                        warmup=maxIt/100)
      Ct <- updt[[1]]
      xbar <- updt[[2]] 
    }
    
    if(i %% 1000 == 0) cat(i,"\n")
  }
  ratio <- sum(accept)/(maxIt-1)
  cat("Acceptance ratio: ", ratio,"\n")
  return(list(chain,ratio,lambda))
}

randomWalk <- function(N, x0, maxIt=10000,
                       adaptCov=FALSE, target=NULL) {
  if(N!=length(x0)) stop("Dimension mismatch.")
  
  chain <- matrix(0,maxIt,N)
  
  sigma <- 2.4/sqrt(N)
  Ct <- sigma^2*diag(N) 
  xbar <- x0
  
  accept <- rep(0,maxIt)
  chain[1,] <- x0
  for (i in 2:maxIt){
    if (adaptCov==FALSE) {
      xStar <- rnorm(N,sd=sigma) + chain[i-1,]
      if(log(runif(1)) < sum(target(xStar,distrib = target)) -
         sum(target(chain[i-1,],distrib = target))){
        accept[i] <- 1
        chain[i,] <- xStar
      } else {
        chain[i,] <- chain[i-1,]
      }
    } else { # with covariance
      xStar <- as.vector(t(chol(Ct))%*%rnorm(N) + chain[i-1,])
      if(log(runif(1)) < sum(target(xStar,distrib = target)) -
         sum(target(as.vector(chain[i-1,]),distrib = target))){
        accept[i] <- 1
        chain[i,] <- xStar
      } else {
        chain[i,] <- chain[i-1,]
      }
      updt <- recursion(Ct=Ct,
                        XbarMinus=xbar,
                        Xt=chain[i,],
                        epsilon = 0.000001,
                        sd=sigma^2,
                        t=i,
                        warmup=maxIt/10)
      Ct <- updt[[1]]
      xbar <- updt[[2]]
    }
    
    if(i %% 1000 == 0) cat(i,"\n")
  }
  ratio <- sum(accept)/(maxIt-1)
  cat("Acceptance ratio: ", ratio,"\n")
  if (adaptCov) {
    return(list(chain,ratio,sigma,Ct))
  } else{
    return(list(chain,ratio,sigma,diag(N)))
  }
}


MTM <- function(N, x0, maxIt=10000,
                       adaptCov=FALSE, targetName=NULL) {
  if(N!=length(x0)) stop("Dimension mismatch.")
  
  chain <- matrix(0,maxIt,N)
  
  sigma <- 2.4/sqrt(N)
  Ct <- sigma^2*diag(N) 
  xbar <- x0
  
  accept <- rep(0,maxIt)
  chain[1,] <- x0
  for (i in 2:maxIt){
    if (adaptCov==FALSE) {
      ys   <- matrix(rnorm(N*N,sd=sigma),N,N) + chain[i-1,]
      targStars <- target(t(ys),distrib = targetName)
      yStar     <- as.vector(ys[,sample(1:N,1,prob = targStars)])
      xs   <- cbind(matrix(rnorm(N*(N-1),sd=sigma),N,N-1) + yStar, chain[i-1,])
      allTargets <- target(t(cbind(xs,ys)),distrib = targetName)
      xTargsSum <- sum( allTargets[1:N] )
      yTargsSum <- sum( allTargets[(N+1):(2*N)] )
      
      if(runif(1) < yTargsSum / xTargsSum){
        accept[i] <- 1
        chain[i,] <- yStar
      } else {
        chain[i,] <- chain[i-1,]
      }
    } else { # with covariance
      ys   <- t(chol(Ct))%*%matrix(rnorm(N*N,sd=sigma),N,N) + chain[i-1,]
      targStars <- target(t(ys),distrib = targetName)
      yStar     <- as.vector(ys[,sample(1:N,1,prob = targStars)])
      xs        <- cbind(t(chol(Ct))%*%matrix(rnorm(N*(N-1),sd=sigma),N,N-1) + yStar, chain[i-1,])
      allTargets <- target(t(cbind(xs,ys)),distrib = targetName)
      xTargsSum <- sum( allTargets[1:N] )
      yTargsSum <- sum( allTargets[(N+1):(2*N)] )
      
      if(runif(1) < yTargsSum / xTargsSum){
        accept[i] <- 1
        chain[i,] <- yStar
      } else {
        chain[i,] <- chain[i-1,]
      }
      updt <- recursion(Ct=Ct,
                        XbarMinus=xbar,
                        Xt=chain[i,],
                        epsilon = 0.000001,
                        sd=sigma^2,
                        t=i,
                        warmup=maxIt/10)
      Ct <- updt[[1]]
      xbar <- updt[[2]]
    }
    
    if(i %% 1000 == 0) cat(i,"\n")
  }
  ratio <- sum(accept)/(maxIt-1)
  cat("Acceptance ratio: ", ratio,"\n")
  if (adaptCov) {
    return(list(chain,ratio,sigma,Ct))
  } else{
    return(list(chain,ratio,sigma,diag(N)))
  }
}

# requires hpHawkes
ssHawkes <- function(locations=NULL,
                     times=NULL,
                     params=c(1, 1/1.6, 1/(14*24),1,1,1),
                     latentDimension=2,
                     threads=1,                     # number of CPU cores
                     simd=0,                        # simd = 0, 1, 2 for no simd, SSE, and AVX, respectively
                     gpu=0,
                     single=0,
                     lambda=1,
                     maxIt=10000,
                     adaptStepSize=TRUE,
                     targetAccept=0.35,
                     adaptScales=TRUE) {
  
  N <- 4
  x0 <- c(1, 1/1.6, 1/(14*24),1,1,1)
  
  if(is.null(locations)){
    stop("No locations found.")
  }
  if(is.null(times)){
    stop("No times found.")
  }
  numLocs <- dim(locations)[1]
  P <- latentDimension
  
  chain <- matrix(0,maxIt,N)
  
  Acceptances = 0 # total acceptances within adaptation run (<= SampBound)
  SampBound = 5   # current total samples before adapting radius
  SampCount = 0   # number of samples collected (adapt when = SampBound)
  Proposed = 0
  
  if (adaptScales) {
    Ct <- diag(N)
    xbar <- rep(0,N) 
  } else {
    Ct <- NULL
  }
  
  engine <- engineInitial(locations,numLocs,P,times,x0,threads,simd,gpu,single)
  engine <- hpHawkes::setParameters(engine,x0)
  
  accept <- rep(0,maxIt)
  chain[1,] <- log(x0[c(1,4:6)])
  for (i in 2:maxIt){
    Proposed = Proposed + 1
    
    v_Nplus1 <- rep((1+sqrt(N+1))/N,N)
    M <- rbind(diag(N), v_Nplus1)
    M <- M - matrix(v_Nplus1,N+1,N)
    M4 <- M * lambda /sqrt(2)
    U <- pracma::randortho(n=N)
    if (adaptScales) {
      M4 <- M4 %*% U %*% chol(Ct)
    } else {
      M4 <- M4 %*% U
    }
    M4 <- M4 + matrix(chain[i-1,],N+1,N,byrow = TRUE)
    hawkesTargets <- apply(M4,MARGIN=1,
                           FUN=function(x) {
                             hpHawkes::setParameters(engine,c(exp(x[1]),
                                                              x0[2],
                                                              x0[3],
                                                              exp(x[2:4])))
                             output <- -hpHawkes::Potential(engine,
                                                                c(exp(x[1]),
                                                                  x0[2],
                                                                  x0[3],
                                                                  exp(x[2:4])))
                             return(output)})
    hawkesTargets <- exp(hawkesTargets-logsumexp(hawkesTargets))
    
    prpsl <- M4[sample(1:(N+1),1,prob = hawkesTargets),]
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
    
    if (adaptScales) {
      updt <- recursion(Ct=Ct,     # multivariate scales update
                        XbarMinus=xbar,
                        Xt=chain[i,],
                        epsilon = 0.000001,
                        sd=1,
                        t=i,
                        warmup=20)
      Ct <- updt[[1]]
      xbar <- updt[[2]] 
    }
    
    if(i %% 100 == 0) cat(i,"\n")
  }
  ratio <- sum(accept)/(maxIt-1)
  cat("Acceptance ratio: ", ratio,"\n")
  return(list(chain,ratio,lambda))
}


