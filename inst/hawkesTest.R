setwd("~/simplicialSampler/")

source("R/simplicialSampler.R")
library(hpHawkes)

load("inst/data/no_holidays.Rdata")

X <- as.matrix(df_no_holidays[,1:2])
times <- df_no_holidays$Time

# fixed spatial smoother (2nd parameter) corresponds to 1.6 km lengthscale
# fixed temporal smoother (3rd parameter) corresponds to 14 day lenghtscale
Max <- 10
burn <- 2
set.seed(666)

res <- sampler(n_iter=Max,burnIn=burn,locations = X,
               params = c(10/1.605,1/c(1.605,14*24),10/(14*24),0.5,0.5),
               times=times,gpu=2,radius=0.1)

res2 <- ssHawkes(maxIt = Max+burn, locations = X, times=times, lambda = 3,
                 simd = 2,threads = 2)

#save(samps,file = "data/no_holidays_samples_new_param.Rdata")


