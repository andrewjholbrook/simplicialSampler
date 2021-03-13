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
df1 <- lambdaSearch[(lambdaSearch$Scaled=="unscaled" & 
                      lambdaSearch$TargetDistribution=="sphericalGaussian"),]

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

for (i in 1:dim(df1)[1]) {
  if(df1$isMaxMinESS[i] == TRUE &
     df1$isMaxMeanESS[i] == TRUE) {
    dummyRow <- df1[i,]
    dummyRow$isMaxMinESS <- FALSE
    df1 <- rbind(df1, dummyRow)
  }
}
df2 <- df1


#
###
####### ill conditioned gaussian scaled
###
#

df1 <- lambdaSearch[(lambdaSearch$Scaled=="scaled" & 
                       lambdaSearch$TargetDistribution=="diagGaussian"),]

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

for (i in 1:dim(df1)[1]) {
  if(df1$isMaxMinESS[i] == TRUE &
     df1$isMaxMeanESS[i] == TRUE) {
    dummyRow <- df1[i,]
    dummyRow$isMaxMinESS <- FALSE
    df1 <- rbind(df1, dummyRow)
  }
}

df1 <- rbind(df1,df2)

colnames(df1)[9] <- "Criterion:"
df1$`Gaussian:` <- factor(df1$TargetDistribution,labels = c("Ill-conditioned","Spherical"))

gg <- ggplot(data = df1, aes(x=Dimension,y=Acceptance,color=`Criterion:`,
                             shape=`Gaussian:`)) +
  geom_hline(yintercept = 0.675) +
  annotate(label="0.675",x=20,y=0.69,color="black",geom="text") +
  geom_jitter() +
  stat_smooth() +
  theme_bw() +
  ylab("Target acceptance") +
  scale_color_discrete(labels=c("Max-Min ESS","Max-Mean ESS")) 

gg

gg2 <- ggplot(data = df1, aes(x=Dimension,y=lambda,color=`Criterion:`,
              shape=`Gaussian:`)) +
  geom_jitter() +
  stat_smooth() +
  theme_bw() +
  ylab("Lambda") +
  scale_color_discrete(labels=c("Max-Min ESS","Max-Mean ESS")) 

gg2

library(grid)
library(gridExtra)

source("inst/grid_arrange.R")
ggsave(grid_arrange_shared_legend(gg,gg2),file="inst/figures/acceptFigOrig.pdf",
       width = 9,height = 4)

system2(command = "pdftk", 
        args    = c("~/simplicialSampler/inst/figures/acceptFigOrig.pdf",
                    "cat 2-end",
                    "output ~/simplicialSampler/inst/figures/acceptFig.pdf") 
)
system2(command = "pdfcrop", 
        args    = c("~/simplicialSampler/inst/figures/acceptFig.pdf", 
                    "~/simplicialSampler/inst/figures/acceptFig.pdf") 
)
system2(command = "rm", 
        args    = c("~/simplicialSampler/inst/figures/acceptFigOrig.pdf") 
)


