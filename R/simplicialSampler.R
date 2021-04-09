
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


proposal <- function(N,x,lambda,distrib,adaptScales=FALSE,Ct=NULL,
                     gaussians=FALSE,nProps=NULL) {
  # N is dimensionality
  # x is current state
  # lambda is scale
  # distrib is target distribution
  if(is.null(nProps)) {
    v_Nplus1 <- rep((1+sqrt(N+1))/N,N)
    M <- rbind(diag(N), v_Nplus1)
    M <- M - matrix(v_Nplus1,N+1,N)
    M4 <- M * lambda /sqrt(2) 
    if(gaussians) M4 <- M4 * sqrt(rchisq(n=1,df=N))
    U <- pracma::randortho(n=N)
    if (adaptScales) {
      M4 <- M4 %*% U %*% chol(Ct)
    } else {
      M4 <- M4 %*% U
    }
    M4 <- M4 + matrix(x,N+1,N,byrow = TRUE)
    output <- M4[sample(1:(N+1),1,prob = target(M4,distrib)),]
  } else { # fewer proposals (bad idea)
    if(nProps>N) stop("nProps greater than N.")
    v_Nplus1 <- rep((1+sqrt(N+1))/N,N)
    M <- rbind(diag(N)[1:nProps,], v_Nplus1)
    M <- M - matrix(v_Nplus1,nProps+1,N)
    M4 <- M * lambda /sqrt(2) 
    if(gaussians) M4 <- M4 * sqrt(rchisq(n=1,df=N))
    U <- pracma::randortho(n=N)
    if (adaptScales) {
      M4 <- M4 %*% U %*% chol(Ct)
    } else {
      M4 <- M4 %*% U
    }
    M4 <- M4 + matrix(x,nProps+1,N,byrow = TRUE)
    output <- M4[sample(1:(nProps+1),1,prob = target(M4,distrib)),]
  }
  return(output)
}

simplicialSampler <- function(N,x0,lambda=1,maxIt=10000,adaptStepSize=TRUE,
                              targetAccept=0.5,target=NULL,adaptScales=FALSE,
                              Gaussians=FALSE,nProps=NULL) {
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
        if(is.null(nProps)) {
          prpsl <- proposal(N,chain[i-1,],lambda,target,
                            adaptScales = adaptScales, Ct=Ct,
                            gaussians = Gaussians)
        } else {
          prpsl <- proposal(N,chain[i-1,],lambda,target,
                            adaptScales = adaptScales, Ct=Ct,
                            gaussians = Gaussians, nProps=nProps)
        }
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
                        warmup=maxIt/100)
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
      tCholCt <- t(chol(Ct))
      ys   <- tCholCt%*%matrix(rnorm(N*N,sd=1),N,N) + chain[i-1,]
      targStars <- target(t(ys),distrib = targetName)
      yStar     <- as.vector(ys[,sample(1:N,1,prob = targStars)])
      xs        <- cbind(tCholCt%*%matrix(rnorm(N*(N-1),sd=1),N,N-1) + yStar, chain[i-1,])
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
                        warmup=maxIt/100)
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

#
### GP functions
#

gpReg = function(x, y, nIter, sampler="simpl",
                 targetAccept=NULL, precond=FALSE){
  
  n <- dim(x)[1]
  d2 <- dim(x)[2]
  
  diffMatAll2 <- matrix(0,n,n)
  for(k in 1:d2) {
    diffMatAll = matrix(x[,k], nrow=n, ncol=n) - matrix(x[,k], nrow=n, ncol=n, byrow=TRUE)
    diffMatAll2 = diffMatAll^2 + diffMatAll2
  }
  
  post.lambda = rep(NA, nIter)
  post.eta = rep(NA, nIter)
  post.rho = rep(NA, nIter)
  post.sigma = rep(NA, nIter)
  post.Z = matrix(NA, n, nIter)
  
  lambda = 1
  eta = 1
  rho = 1
  sigma = 1
  Z     = -10*(y-0.5) #rnorm(n,sd=1/sqrt(n))
  
  if (precond) {
    Ct <- diag(n)
    xbar <- rep(0,n) 
  } else {
    Ct <- NULL
  }
  
  tuner = 1
  totalAccept <- rep(0,nIter)
  Acceptances = 0 # total acceptances within adaptation run (<= SampBound)
  SampBound = 5   # current total samples before adapting radius
  SampCount = 0   # number of samples collected (adapt when = SampBound)
  Proposed = 0
  
  for(i in 1:nIter){
    
    lambda <- get.lambda(x, y, Z, diffMatAll2, lambda, eta, rho, sigma)
    eta <- get.eta(x, y, Z, diffMatAll2, lambda, eta, rho, sigma)
    rho <- get.rho(x, y, Z, diffMatAll2, lambda, eta, rho, sigma)
    sigma <- get.sigma(x, y, Z, diffMatAll2, lambda, eta, rho, sigma)
    
    prpsl <- NULL
    attempt <- 0
    while( is.null(prpsl) && attempt <= 100 ) {
      attempt <- attempt + 1
      try(
        prpsl <- get.Z(x, y, Z, diffMatAll2, lambda,
                       eta, rho, sigma, sampler, tuner, Ct)
      )
    }
    Z <- prpsl
    
    post.lambda[i] <- lambda	
    post.eta[i] <- eta
    post.rho[i] <- rho
    post.sigma[i] <- sigma
    post.Z[,i] <- Z
    
    if (precond) {
      updt <- recursion(Ct=Ct,
                        XbarMinus=xbar,
                        Xt=Z,
                        epsilon = 0.000001,
                        sd=1,
                        t=i,
                        warmup=nIter/100)
      Ct <- updt[[1]]
      xbar <- updt[[2]] 
    }
    
    SampCount <- SampCount + 1
    if (i>1){
      if ( any(Z != post.Z[,i-1]) ) {
        totalAccept[i] <- 1
        Acceptances = Acceptances + 1
      }
    }
    
    if (SampCount == SampBound) { # adjust tuner at increasing intervals
      AcceptRatio <- Acceptances / SampBound
      AdaptRatio  <- AcceptRatio / targetAccept
      if (AdaptRatio>2) AdaptRatio <- 2
      if (AdaptRatio<0.5) AdaptRatio <- 0.5
      tuner <- tuner * AdaptRatio
      
      SampCount <- 0
      SampBound <- ceiling(SampBound^1.01)
      Acceptances <- 0
    }
    
    if(i %% 100 == 0){
      cat("Iteration ", i, "\n")
    }
  }
  
  
  return(list(lambda = post.lambda, eta = post.eta, rho = post.rho,
              sigma = post.sigma, Z = post.Z,
              accept.rate = sum(totalAccept)/(nIter-1)))
  
}

get.Z <- function(x, y, Z, diffMatAll2, lambda, eta, rho, sigma, sampler,
                  tuner=1, Ct=NULL){
  
  if(sampler=="simpl") { # simplicial sampler
    N <- length(Z)
    v_Nplus1 <- rep((1+sqrt(N+1))/N,N)
    M <- rbind(diag(N), v_Nplus1)
    M <- M - matrix(v_Nplus1,N+1,N)
    M4 <- M * tuner /sqrt(2) 
    U <- pracma::randortho(n=N)
    if (!is.null(Ct)) {
      M4 <- M4 %*% U %*% chol(Ct)
    } else {
      M4 <- M4 %*% U
    }
    M4 <- M4 + matrix(Z,N+1,N,byrow = TRUE)
    newParam <- M4[sample(1:(N+1),1,
                          prob = exp(getPost(x, y, t(M4), diffMatAll2,
                                             lambda, eta, rho, sigma))),]
    
  } else if(sampler=="RWM") { # mh
    N <- length(Z)
    
    if (! is.null(Ct)) {
      Zstar <- Z + t(chol(Ct))%*%rnorm(N,sd=1/sqrt(N)*tuner)
    } else {
      Zstar <- Z + rnorm(N,sd=1/sqrt(N)*tuner)
    }
    posts <- getPost(x, y, cbind(Z,Zstar), diffMatAll2, lambda, eta, rho, sigma)
    
    if(log(runif(1)) < posts[2]-posts[1]) {
      newParam <- Zstar
    } else {
      newParam <- Z
    }
    
  } else { #MTM
    N <- length(Z)
    if (! is.null(Ct)) {
      tCholCt <- t(chol(Ct))
      zs   <- tCholCt%*%matrix(rnorm(N*N,sd=1/sqrt(N)*tuner),N,N) + Z
    } else {
      zs   <- matrix(rnorm(N*N,sd=1/sqrt(N)*tuner),N,N) + Z
    }
    postAndChol <- getPost(x, y, zs, diffMatAll2, lambda, eta, rho, sigma, returnL = TRUE, Exp=TRUE) # so we don't need to decompose twice
    targStars <- postAndChol[[1]]
    L         <- postAndChol[[2]]
    ZStar     <- as.vector(zs[,sample(1:N,1,prob = targStars)])
    if (! is.null(Ct)) {
      zPrimes   <- cbind(tCholCt%*%matrix(rnorm(N*(N-1),sd=1/sqrt(N)*tuner),N,N-1) + ZStar, Z)
    } else {
      zPrimes   <- cbind(matrix(rnorm(N*(N-1),sd=1/sqrt(N)*tuner),N,N-1) + ZStar, Z)
    }
    allTargets <- getPost(x, y, zPrimes, diffMatAll2, lambda, eta, rho, sigma, L=L, Exp = TRUE)
    
    if(runif(1) < sum(targStars)/sum(allTargets)){
      newParam <- ZStar
    } else {
      newParam <- Z
    }
    
  }
  
  
  return(newParam)
}


get.lambda <- function(x, y, Z, diffMatAll2, lambda, eta, rho, sigma){
  
  w = 4
  m = 10
  
  z = getPost(x, y, Z, diffMatAll2, lambda, eta, rho, sigma) - rexp(1)
  
  # Stepping out to obtain the [L, R] range
  u = runif(1)
  L = lambda - w*u
  R = L + w
  v = runif(1)
  J = floor(m*v)
  K = (m-1) - J
  
  L = max(0, L)
  while (J>0 && L>0 && z < getPost(x, y, Z, diffMatAll2, L, eta, rho, sigma)) {
    L = L - w	
    L = max(0, L)		
    J = J - 1
  }
  
  while (K>0 && z < getPost(x, y, Z, diffMatAll2, R, eta, rho, sigma)) {
    R = R+w
    K = K-1
  }
  
  
  # Shrinkage to obtain a sample
  u = runif(1)
  newParam = L + u*(R-L)
  
  while (z > getPost(x, y, Z, diffMatAll2, newParam, eta, rho, sigma)) {
    if (newParam < lambda) {
      L = newParam
    }else{
      R = newParam
    }
    
    u = runif(1)
    newParam = L + u*(R-L)
  }
  
  
  return(newParam)
}




get.eta <- function(x, y, Z, diffMatAll2, lambda, eta, rho, sigma){
  
  w = 4
  m = 10
  
  z = getPost(x, y, Z, diffMatAll2, lambda, eta, rho, sigma) - rexp(1)
  
  # Stepping out to obtain the [L, R] range
  u = runif(1)
  L = eta - w*u
  R = L + w
  v = runif(1)
  J = floor(m*v)
  K = (m-1) - J
  
  L = max(0, L)
  while (J>0 && L>0 && z < getPost(x, y, Z, diffMatAll2, lambda, L, rho, sigma)) {
    L = L - w	
    L = max(0, L)		
    J = J - 1
  }
  
  while (K>0 && z < getPost(x, y, Z, diffMatAll2, lambda, R, rho, sigma)) {
    R = R+w
    K = K-1
  }
  
  
  # Shrinkage to obtain a sample
  u = runif(1)
  newParam = L + u*(R-L)
  
  while (z > getPost(x, y, Z, diffMatAll2, lambda, newParam, rho, sigma)) {
    if (newParam < eta) {
      L = newParam
    }else{
      R = newParam
    }
    
    u = runif(1)
    newParam = L + u*(R-L)
  }
  
  return(newParam)
}







get.rho <- function(x, y, Z, diffMatAll2, lambda, eta, rho, sigma){
  
  w = 4
  m = 10
  
  z = getPost(x, y, Z, diffMatAll2, lambda, eta, rho, sigma) - rexp(1)
  
  # Stepping out to obtain the [L, R] range
  u = runif(1)
  L = rho - w*u
  R = L + w
  v = runif(1)
  J = floor(m*v)
  K = (m-1) - J
  
  L = max(0, L)
  while (J>0 && L>0 && z < getPost(x, y, Z, diffMatAll2, lambda, eta, L, sigma)) {
    L = L - w	
    L = max(0, L)		
    J = J - 1
  }
  
  while (K>0 && z < getPost(x, y, Z, diffMatAll2, lambda, eta, R, sigma)) {
    R = R+w
    K = K-1
  }
  
  
  # Shrinkage to obtain a sample
  u = runif(1)
  newParam = L + u*(R-L)
  
  while (z > getPost(x, y, Z, diffMatAll2, lambda, eta, newParam, sigma)) {
    if (newParam < rho) {
      L = newParam
    }else{
      R = newParam
    }
    
    u = runif(1)
    newParam = L + u*(R-L)
  }
  
  return(newParam)
}





get.sigma <- function(x, y, Z, diffMatAll2, lambda, eta, rho, sigma){
  
  w = 4
  m = 10
  
  z = getPost(x, y, Z, diffMatAll2, lambda, eta, rho, sigma) - rexp(1)
  
  # Stepping out to obtain the [L, R] range
  u = runif(1)
  L = sigma - w*u
  R = L + w
  v = runif(1)
  J = floor(m*v)
  K = (m-1) - J
  
  L = max(0, L)
  while (J>0 && L>0 && z < getPost(x, y, Z, diffMatAll2, lambda, eta, rho, L)) {
    L = L - w	
    L = max(0, L)		
    J = J - 1
  }
  
  while (K>0 && z < getPost(x, y, Z, diffMatAll2, lambda, eta, rho, R)) {
    R = R+w
    K = K-1
  }
  
  
  # Shrinkage to obtain a sample
  u = runif(1)
  newParam = L + u*(R-L)
  
  while (z > getPost(x, y, Z, diffMatAll2, lambda, eta, rho, newParam)) {
    if (newParam < sigma) {
      L = newParam
    }else{
      R = newParam
    }
    
    u = runif(1)
    newParam = L + u*(R-L)
  }
  
  return(newParam)
}




getPost = function(x, y, Z, diffMatAll2, lambda, eta, rho, sigma,
                   L=NULL, returnL=FALSE,Exp=FALSE){
  
  if(is.null(L)) {
    C = lambda + eta*exp(-rho*(diffMatAll2)) + sigma*diag(1, length(y));
    L = chol(C); 
  }
  
  detC = prod(diag(L))^2
  
  ncols <- ncol(Z)
  if(is.null(ncols)){ 
    invL.Z = backsolve(L,Z)
    
    prbs <- 1/(1+exp(-Z)) #pnorm(Z)
    logLike = sum(y*log(prbs) + (1-y)*log(1-prbs))
    
    logPrior =  dgamma(lambda, 1, 1, log=TRUE) +
      dgamma(eta, 1, 1, log=TRUE) + dgamma(rho, 1, 1, log=TRUE) +
      dgamma(sigma, 1, 1, log=TRUE) -
      0.5*log(detC) - 0.5 * t(invL.Z)%*%invL.Z
    
    
    logPost = (logLike + logPrior)
    
  } else { # multiple evals 
    
    logPost <- rep(0,ncols)
    logPrior1 =  dgamma(lambda, 1, 1, log=TRUE) +
      dgamma(eta, 1, 1, log=TRUE) + dgamma(rho, 1, 1, log=TRUE) +
      dgamma(sigma, 1, 1, log=TRUE) -
      0.5*log(detC)
    
    invL.Z = backsolve(L,Z)
    
    for (k in 1:ncols) {
      prbs <- 1/(1+exp(-Z[,k])) #pnorm(Z[,k])
      logLike = sum(y*log(prbs) + (1-y)*log(1-prbs))
      if(any(is.na(logLike))) browser()
      
      logPrior2 = - 0.5 * t(invL.Z[,k])%*%invL.Z[,k]
      
      
      logPost[k] = (logLike + logPrior1 + logPrior2)
    }
  }
  
  if(Exp==TRUE) logPost <- exp(logPost) # for MTM
  
  if (returnL) {
    return(list(logPost,L))
  } else {
    return(logPost)
  }
}

