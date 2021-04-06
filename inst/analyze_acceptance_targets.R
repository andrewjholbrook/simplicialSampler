setwd("~/simplicialSampler/")

library(readr)
library(ggplot2)
library(wesanderson)
pal <- wes_palette("FantasticFox1", 5, type = "discrete")

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
df1$`Gaussian target:` <- factor(df1$TargetDistribution,labels = c("Ill-conditioned","Spherical"))

gg <- ggplot(data = df1, aes(x=Dimension,y=Acceptance,color=`Criterion:`,
                             shape=`Gaussian target:`)) +
  annotate(label="0.675",x=20,y=0.69,color="black",geom="text") +
  geom_jitter() +
  stat_smooth(se=FALSE) +
  geom_hline(yintercept = 0.675) +
  theme_bw() +
  ylab("Target acceptance") +
  ggtitle("Emprically optimal acceptance rates")+
  scale_color_manual(labels=c("Min ESS","Mean ESS"),
                     values=c(pal[3],pal[5])) 

gg

gg2 <- ggplot(data = df1, aes(x=Dimension,y=lambda,color=`Criterion:`,
              shape=`Gaussian target:`)) +
  geom_jitter() +
  stat_smooth(se=FALSE) +
  theme_bw() +
  ylab("Lambda") +
  ggtitle("Empirically optimal edge lengths")+
  scale_color_manual(labels=c("Min ESS","Mean ESS"),
                     values=c(pal[3],pal[5])) 
gg2


df <- read_table2("inst/output/nPropsResults.txt", 
                  col_names = FALSE)
df <- df[,1:3]
colnames(df) <- c("nProps","Mean ESS","Min ESS")
df <- reshape2::melt(df,measure.vars=2:3)
df$nProps <- factor(df$nProps)
df$variable <- factor(df$variable)

gg3 <- ggplot(df,aes(x=nProps,y=value,fill=variable,color=variable)) +
  geom_boxplot()+
  ylab("Effective sample size (ESS)") +
  ggtitle("Varying P for 100-D spherical target") +
  xlab("Number of proposals (P)") +
  scale_x_discrete(breaks = c(10,20,30,40,50,60,70,80,90,100)) +
  scale_fill_manual(values=c(pal[5],pal[3])) +
  scale_color_manual(values=c(pal[5],pal[3])) +
  theme_bw()
gg3

library(grid)
library(gridExtra)

source("inst/grid_arrange.R")
ggsave(grid_arrange_shared_legend(gg,gg2,gg3),file="inst/figures/acceptFigOrig.pdf",
       width = 12,height = 4)

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


