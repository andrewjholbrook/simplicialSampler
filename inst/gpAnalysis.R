setwd("~/simplicialSampler/")

library(readr)

df <- read_table2("inst/output/gpClassification.txt", 
                            col_names = FALSE)
df <- df[,1:9]
colnames(df) <- c("Algorithm","MeanEssLatents","MinEssLatents","EssLambda",
                            "EssEta","EssRho", "EssSigma","ItsTo10", "Time")

df[,2:7] <- df[,2:7] / df$Time
df$TimeTo10 <- df$ItsTo10 * (df$ItsTo10/100000*df[,9])

df <- df[,c(1:7,10)]

df$Algorithm <- as.factor(as.numeric(factor(df$Algorithm)))
g <- df$Algorithm
splt <- split(df,g)
splt <- lapply(splt,FUN=function(x) apply(x, MARGIN = 2, as.numeric ))
sds <- lapply(splt,FUN=function(x) apply(x, MARGIN = 2,FUN=function(x) sd(x)/10 ))
sds[[1]][1] <- "MTM"
sds[[2]][1] <- "RWM"
sds[[3]][1] <- "simpl"


mns <- lapply(splt,FUN=function(x) apply(x, MARGIN = 2, FUN=mean))
mns[[1]][1] <- "MTM"
mns[[2]][1] <- "RWM"
mns[[3]][1] <- "simpl"

df2 <- rbind(sds[[1]],sds[[2]],sds[[3]])
df2 <- as.data.frame(df2)
Algorithm <- factor(df2$Algorithm)
df2 <- df2[,-1]
df2 <- as.matrix(df2) 
df2 <- data.frame(df2)
df2$Algorithm <- Algorithm

df3 <- rbind(mns[[1]],mns[[2]],mns[[3]])
df3 <- as.data.frame(df3)
Algorithm <- factor(df3$Algorithm)
df3 <- df3[,-1]
df3 <- as.matrix(df3)
df3 <- data.frame(df3)
df3$Algorithm <- Algorithm

df4 <- cbind(df3,df2)
df4 <- df4[,-8]
df4 <- df4[,order(colnames(df4))]

df4 <- df4[3:1,]
