setwd("~/simplicialSampler/")

library(readr)
library(ggplot2)

lambdaSearch <- read_table2("inst/output/lambdaSearch.txt", 
                            col_names = FALSE)
lambdaSearch <- lambdaSearch[,1:8]
colnames(lambdaSearch) <- c("Dimension","lambda","meanESS","minESS",
                            "Acceptance","TargetAcceptance",
                            "TargetDistribution", "Scaled")
#
###
####### spherical gaussian unscaled
###
#
df1 <- lambdaSearch[lambdaSearch$Scaled=="unscaled" & 
                      lambdaSearch$TargetDistribution=="sphericalGaussian",]

keep <- vector() # get max min ESSs
for (i in unique(df1$Dimension)) {
  keep <- c(keep,max(df1$minESS[df1$Dimension==i]))
}
df1$isMaxMinESS <- df1$minESS %in% keep

keep <- vector() # get max mean ESSs
for (i in unique(df1$Dimension)) {
  keep <- c(keep,max(df1$meanESS[df1$Dimension==i]))
}
df1$isMaxMeanESS <- df1$meanESS %in% keep

df1 <- df1[df1$isMaxMeanESS | df1$isMaxMinESS,]
df1 <- df1[df1$Dimension > 2,]

# plot(df1$Dimension,df1$lambda)
# plot(df1$Dimension,df1$Acceptance)

gg <- ggplot(data = df1, aes(x=Dimension,y=Acceptance,color=isMaxMinESS)) +
  geom_point() +
  stat_smooth() +
  theme_bw() +
  scale_color_discrete(labels=c("Max-Min ESS","Max-Mean ESS")) +
  theme(legend.title = element_blank())

gg

gg2 <- ggplot(data = df1, aes(x=Dimension,y=lambda,color=isMaxMinESS)) +
  geom_point() +
  stat_smooth() +
  theme_bw() +
  scale_color_discrete(labels=c("Max-Min ESS","Max-Mean ESS")) +
  theme(legend.title = element_blank())

gg2
