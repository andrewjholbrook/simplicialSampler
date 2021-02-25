setwd("~/hedgingMCMC/")

library(readr)

lambdaSearch <- read_table2("output/lambdaSearch.txt", 
                            col_names = FALSE)
lambdaSearch <- lambdaSearch[,1:6]
colnames(lambdaSearch) <- c("Dimension","lambda","meanESS","minESS",
                            "Acceptance","TargetAcceptance")
keep <- vector()
for (i in unique(lambdaSearch$Dimension)) {
  keep <- c(keep,max(lambdaSearch$minESS[lambdaSearch$Dimension==i]))
}

df <- lambdaSearch[lambdaSearch$minESS %in% keep,]
View(df)

df <- df[order(df$Dimension),]
plot(df$Dimension,df$lambda)
lines(x=seq(from=0,to=64),y=seq(from=0,to=64)^(1/6))
mod <- lm(df$lambda ~ df$Dimension + I(df$Dimension^2))
