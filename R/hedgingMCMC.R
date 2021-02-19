library(rstiefel)
library(mvtnorm)

f <- function(X) {
  densities <- dmvnorm(X,sigma = diag(1:dim(X)[2]))
  return(densities)
}

proposal <- function(N,x,lambda) {
  # N is dimensionality
  # x is current state
  # lambda is scale
  v_Nplus1 <- rep((1+sqrt(N+1))/N,N)
  M <- rbind(diag(N), v_Nplus1)
  M <- M - matrix(v_Nplus1,N+1,N)
  M4 <- M * lambda /sqrt(2)
  U <- rmf.matrix(matrix(0,N,N))
  M4 <- M4 %*% U
  M4 <- M4 + matrix(x,N+1,N,byrow = TRUE)
  
  output <- M4[sample(1:(N+1),1,prob = f(M4)),]
  
  return(output)
}

mcmc <- function(N,x0,lambda,maxIt=10000) {
  if(N!=length(x0)) stop("Dimension mismatch.")

  chain <- matrix(0,maxIt,N)
  accept <- rep(0,maxIt)
  chain[1,] <- x0
  for (i in 2:maxIt){
    chain[i,] <- proposal(N,chain[i-1,],lambda)
    if(any(chain[i,] != chain[i-1,])) accept[i] <- 1 
    if(i %% 100 == 0) cat(i,"\n")
  }
  cat("Acceptance ratio: ", sum(accept)/(maxIt-1),"\n")
  return(chain)
}

#
######
############ scratch
#####
#
N = 10
maxIt <- 10000

output <- mcmc(N=N,x0=rep(0,N),lambda = 3, maxIt = maxIt)

for(i in 1:N){
  qqplot(output[,i],rnorm(maxIt,sd=sqrt(i)))
  qqline(output[,i], distribution = function(p) qnorm(p, sd=sqrt(i)))
}


