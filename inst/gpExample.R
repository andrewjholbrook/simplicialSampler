library(mvtnorm)
library(coda)

n = 50
set.seed(1)
x = runif(n, -3, 3)
#y = (1 + 2*x - 3*sin(x)) + rnorm(n, 0, .5)
y = .5*(x^3) - 3*x + 1 + rnorm(n, 0, .5)
y <- rbinom(n=50,size=1,prob=pnorm(y))
plot(x, y, xlim=c(-3, 3))



nIter = 100000

# No PC: 0.234 MH, 0.5 Simpl
output <- gpReg(x=x,
                y=y,
                nIter=nIter,
                sampler="simpl",
                targetAccept = 0.5,
                adaptCov = FALSE)
mcmc.obj <- as.mcmc(t(output$Z))
summary(effectiveSize(mcmc.obj))
traceplot(as.mcmc(output$lambda))
traceplot(as.mcmc(output$eta))
traceplot(as.mcmc(output$rho))
traceplot(as.mcmc(output$sigma[20:nIter]))
traceplot(as.mcmc(mcmc.obj[,1]))



gpReg = function(x, y, nIter, sampler="simpl", targetAccept=NULL, adaptCov=FALSE){
  
  
  diffMatAll = matrix(x, nrow=n, ncol=n) - matrix(x, nrow=n, ncol=n, byrow=TRUE)
  diffMatAll2 = diffMatAll^2
  
  post.lambda = rep(NA, nIter)
  post.eta = rep(NA, nIter)
  post.rho = rep(NA, nIter)
  post.sigma = rep(NA, nIter)
  post.Z = matrix(NA, 50, nIter)
  
  lambda = 1
  eta = 1
  rho = 1
  sigma = 1
  Z     = rnorm(50,sd=1/sqrt(50))
  
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
                       eta, rho, sigma, sampler, tuner)
      )
    }
    Z <- prpsl

    post.lambda[i] <- lambda	
    post.eta[i] <- eta
    post.rho[i] <- rho
    post.sigma[i] <- sigma
    post.Z[,i] <- Z
    
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
                  tuner=1){
  
  if(sampler=="simpl") { # simplicial sampler
    N <- length(Z)
    v_Nplus1 <- rep((1+sqrt(N+1))/N,N)
    M <- rbind(diag(N), v_Nplus1)
    M <- M - matrix(v_Nplus1,N+1,N)
    M4 <- M * tuner /sqrt(2) 
    U <- pracma::randortho(n=N)
    M4 <- M4 %*% U
    M4 <- M4 + matrix(Z,N+1,N,byrow = TRUE)
    newParam <- M4[sample(1:(N+1),1,
                          prob = exp(getPost(x, y, t(M4), diffMatAll2,
                                         lambda, eta, rho, sigma))),]
    
  } else if(sampler=="RWM") { # mh
    Zstar <- Z + rnorm(50,sd=1/sqrt(50)*tuner)
    posts <- getPost(x, y, cbind(Z,Zstar), diffMatAll2, lambda, eta, rho, sigma)
    
    if(log(runif(1)) < posts[2]-posts[1]) {
      newParam <- Zstar
    } else {
      newParam <- Z
    }
    
  } else { #MTM
    
    stop("MTM not implemented yet")
    
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




getPost = function(x, y, Z, diffMatAll2, lambda, eta, rho, sigma){
  
  C = lambda + eta*exp(-rho*(diffMatAll2)) + sigma*diag(1, nrow=n, ncol=n);
  
  L = chol(C);
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
      
      logPrior2 = - 0.5 * t(invL.Z[,k])%*%invL.Z[,k]
      
      
      logPost[k] = (logLike + logPrior1 + logPrior2)
    }
  }
  return(logPost)
}

