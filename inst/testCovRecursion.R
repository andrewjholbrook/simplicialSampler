N <- 5
its <- 100000
sd <- 2.4^2/N
warmup <- 10
M <- matrix(rnorm(N^2),N,N)
S <- M%*%t(M)
Schol <- t(chol(S))
epsilon <- 0.000001
xbar <- 0
Ct <- sd*diag(N) # C0

for (t in 1:its) {
  x <- Schol %*% rnorm(N)
  updt <- recursion(Ct=Ct,
                    XbarMinus=xbar,
                    Xt=x,
                    epsilon = epsilon,
                    sd=sd,
                    t=t,
                    warmup=1000)
  xbar <- updt[[2]]
  Ct   <- updt[[1]]
}
