setwd("~/simplicialSampler/")

source("R/simplicialSampler.R")
library(hpHawkes)
library(coda)

load("inst/data/no_holidays.Rdata")

X <- as.matrix(df_no_holidays[,1:2])
times <- df_no_holidays$Time

# fixed spatial smoother (2nd parameter) corresponds to 1.6 km lengthscale
# fixed temporal smoother (3rd parameter) corresponds to 14 day lenghtscale
Max <- 10000
burn <- 0
set.seed(666)

ptm <- proc.time()
res <- sampler(n_iter=Max,burnIn=burn,locations = X,
                 params=c(1, 1/1.6, 1/(14*24),1,1,1),
                 times=times,gpu=2,radius=2)
time1 <- proc.time() - ptm

mcmc.obj1 <- as.mcmc(t(res$samples[c(1,4:6),200:Max]))
plot(mcmc.obj1)
eff1 <- effectiveSize(mcmc.obj1)
eff1/time1[3]

ptm <- proc.time()
res2 <- ssHawkes(maxIt = Max, locations = X, times=times, lambda = 1,
                 gpu=2)
time2 <- proc.time() - ptm
mcmc.obj2 <- as.mcmc(exp(res2[[1]][100:Max,]))
plot(mcmc.obj2)
eff2 <- effectiveSize(mcmc.obj2)
eff2/time2[3]

#save(samps,file = "data/no_holidays_samples_new_param.Rdata")


