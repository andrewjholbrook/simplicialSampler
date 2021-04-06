library(mvtnorm)
library(coda)

n = 50
set.seed(2)
x = runif(n, -2.5, 2.5)
#y = (1 + 2*x - 3*sin(x)) + rnorm(n, 0, .5)
y = .5*(x^3) - 3*x + 1 + rnorm(n, 0, .5)
y <- rbinom(n=50,size=1,prob=pnorm(y))
plot(x, y, xlim=c(-3, 3))



nIter = 2000


output <- gpReg(x,y,nIter,FALSE)
mcmc.obj <- as.mcmc(t(output$Z))
effectiveSize(mcmc.obj)
traceplot(as.mcmc(output$lambda))
traceplot(as.mcmc(output$eta))
traceplot(as.mcmc(output$rho))
traceplot(as.mcmc(output$sigma[20:2000]))
traceplot(as.mcmc(mcmc.obj[,50]))



gpReg = function(x, y, nIter, simpl=TRUE){
  
  
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
  Z     = rep(0,50)
  
  for(i in 1:nIter){
    
    lambda <- get.lambda(x, y, Z, diffMatAll2, lambda, eta, rho, sigma)
    
    eta <- get.eta(x, y, Z, diffMatAll2, lambda, eta, rho, sigma)
    
    rho <- get.rho(x, y, Z, diffMatAll2, lambda, eta, rho, sigma)
    
    sigma <- get.sigma(x, y, Z, diffMatAll2, lambda, eta, rho, sigma)
    
    Z <- get.Z(x, y, Z, diffMatAll2, lambda, eta, rho, sigma, simpl)
    
    post.lambda[i] <- lambda	
    post.eta[i] <- eta
    post.rho[i] <- rho
    post.sigma[i] <- sigma
    post.Z[,i] <- Z
    
    if(i %% 100 == 0) cat("Iteration ", i, "\n")
  }
  
  
  return(list(lambda = post.lambda, eta = post.eta, rho = post.rho,
              sigma = post.sigma, Z = post.Z))
  
}

get.Z <- function(x, y, Z, diffMatAll2, lambda, eta, rho, sigma, simpl){
  
  if(simpl) { # simplicial sampler
    
  } else { # mh
    Zstar <- Z + rnorm(50,sd=1/50)
    posts <- getPost(x, y, cbind(Z,Zstar), diffMatAll2, lambda, eta, rho, sigma)
    
    if(log(runif(1)) < posts[2]-posts[1]) {
      newParam <- Zstar
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




getPost = function(x, y, Z, diffMatAll2, lambda, eta, rho, sigma){
  
  C = lambda + eta*exp(-rho*(diffMatAll2)) + sigma*diag(1, nrow=n, ncol=n);
  
  L = chol(C);
  invL = solve(L);
  invC = invL%*%t(invL)
  
  ncols <- ncol(Z)
  if(is.null(ncols)){ 
    invL.Z = invL%*%Z
    
    prbs <- pnorm(Z)
    logLike = sum(y*log(prbs) + (1-y)*log(1-prbs))
    
    logPrior =  dgamma(lambda, 1, 1, log=TRUE) +
      dgamma(eta, 1, 1, log=TRUE) + dgamma(rho, 1, 1, log=TRUE) +
      dgamma(sigma, 1, 1, log=TRUE) -
      0.5*log(det(C)) - 0.5 * t(invL.Z)%*%invL.Z
    
    
    logPost = (logLike + logPrior)
    
  } else { # multiple evals 
    
    logPost <- rep(0,ncols)
    logPrior1 =  dgamma(lambda, 1, 1, log=TRUE) +
      dgamma(eta, 1, 1, log=TRUE) + dgamma(rho, 1, 1, log=TRUE) +
      dgamma(sigma, 1, 1, log=TRUE) -
      0.5*log(det(C))
    
    for (k in 1:ncols) {
      invL.Z = invL%*%Z[,k]
      
      prbs <- pnorm(Z[,k])
      logLike = sum(y*log(prbs) + (1-y)*log(1-prbs))
      
      logPrior2 = - 0.5 * t(invL.Z)%*%invL.Z
      
      
      logPost[k] = (logLike + logPrior1 + logPrior2)
    }
  }
  return(logPost)
}

