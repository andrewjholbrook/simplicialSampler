library(pracma)
setwd("~/simplicialSampler/")

dims <- c(4,8,16,32,64,128,256,512)

for (i in dims) {
  O <- pracma::randortho(i)
  C <- O %*% diag(i:1) %*% t(O)
  saveRDS(C,file = paste0("inst/savedCovs/cov",i,".rds"))
}