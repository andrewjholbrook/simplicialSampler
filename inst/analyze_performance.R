setwd("~/simplicialSampler/")

library(readr)
library(ggplot2)
library(reshape2)
library(wesanderson)
pal <- wes_palette("FantasticFox1", 5, type = "discrete")

df <- read_table2("inst/output/rwmComparison.txt",
                            col_names = FALSE)
df <- df[,-ncol(df)]
df$X2 <- df$X2 / df$X5
df$X3 <- df$X3 / df$X6
df$X4 <- df$X4 / df$X7
df <- df[,1:4]
colnames(df) <- c("Dimension","Mean","Min","Slowdown")
df <- melt(df, measure.vars=2:3,value.name = "ESS",variable.name = "Statistic:")
df$ESSs <- df$ESS / df$Slowdown
df <- df[,-2]
df <- melt(df, measure.vars=3:4,value.name = "Relative improvement",variable.name = "Criterion:")
df$Statistic <- factor(df$Statistic)
df$Criterion <- factor(df$Criterion)

gg <- ggplot(df, aes(x=Dimension,y=`Relative improvement`,color=`Criterion:`, shape=`Statistic:`)) +
  geom_point() +
  geom_smooth() +
  xlab(NULL) +
  scale_y_continuous(trans = "log2",breaks = c(1,2,4,8,16,32),limits=c(1,32)) +
  ylab("Relative improvement (Simpl/RWM)") +
  ggtitle("Vanilla algorithms / spherical target") +
  scale_color_manual(values=c(pal[3],pal[5])) +
  theme_bw()
gg

df <- read_table2("inst/output/rwmComparison.txt",
                  col_names = FALSE)
df <- df[,-ncol(df)]
df$X2 <- df$X2 / df$X8
df$X3 <- df$X3 / df$X9
df$X4 <- df$X4 / df$X10
df <- df[,1:4]
colnames(df) <- c("Dimension","Mean","Min","Slowdown")
df <- melt(df, measure.vars=2:3,value.name = "ESS",variable.name = "Statistic:")
df$ESSs <- df$ESS / df$Slowdown
df <- df[,-2]
df <- melt(df, measure.vars=3:4,value.name = "Relative improvement",variable.name = "Criterion:")
df$Statistic <- factor(df$Statistic)
df$Criterion <- factor(df$Criterion)

gg4 <- ggplot(df, aes(x=Dimension,y=`Relative improvement`,color=`Criterion:`, shape=`Statistic:`)) +
  geom_point() +
  geom_smooth() +
  scale_y_continuous(trans = "log2",breaks = c(1,2,4,8,16,32,64,128,256,512),limits=c(0.75,512)) +
  ylab("Relative improvement (Simpl/MTM)") +
  scale_color_manual(values=c(pal[3],pal[5])) +
  theme_bw()
gg4


# source("inst/grid_arrange.R")
# ggsave(gg,file="inst/figures/sphereNormFig.pdf",
#        width = 4.5,height = 3)
# 
# system2(command = "pdfcrop",
#         args    = c("~/simplicialSampler/inst/figures/sphereNormFig.pdf",
#                     "~/simplicialSampler/inst/figures/sphereNormFig.pdf")
# )


#
#######
################ ill gaussian
#######
#

library(readr)
library(ggplot2)
library(reshape2)

df <- read_table2("inst/output/scaledRwmComparison.txt",
                  col_names = FALSE)
df <- df[,-ncol(df)]
df$X2 <- df$X2 / df$X5
df$X3 <- df$X3 / df$X6
df$X4 <- df$X4 / df$X7
df <- df[,1:4]
colnames(df) <- c("Dimension","Mean","Min","Slowdown")
df <- melt(df, measure.vars=2:3,value.name = "ESS",variable.name = "Statistic:")
df$ESSs <- df$ESS / df$Slowdown
df <- df[,-2]
df <- melt(df, measure.vars=3:4,value.name = "Relative improvement",variable.name = "Criterion:")
df$Statistic <- factor(df$Statistic)
df$Criterion <- factor(df$Criterion)


gg2 <- ggplot(df, aes(x=Dimension,y=`Relative improvement`,color=`Criterion:`, shape=`Statistic:`)) +
  geom_point() +
  geom_smooth() +
  scale_y_continuous(trans = "log2",breaks = c(1,2,4,8,16,32),limits=c(1,32)) +
  ylab("") +
  xlab(NULL) +
  ggtitle("PC algorithms / diagonal \"ill\" target") +
  scale_color_manual(values=c(pal[3],pal[5])) +
  theme_bw()
gg2

df <- read_table2("inst/output/scaledRwmComparison.txt",
                  col_names = FALSE)
df <- df[,-ncol(df)]
df$X2 <- df$X2 / df$X8
df$X3 <- df$X3 / df$X9
df$X4 <- df$X4 / df$X10
df <- df[,1:4]
colnames(df) <- c("Dimension","Mean","Min","Slowdown")
df <- melt(df, measure.vars=2:3,value.name = "ESS",variable.name = "Statistic:")
df$ESSs <- df$ESS / df$Slowdown
df <- df[,-2]
df <- melt(df, measure.vars=3:4,value.name = "Relative improvement",variable.name = "Criterion:")
df$Statistic <- factor(df$Statistic)
df$Criterion <- factor(df$Criterion)


gg5 <- ggplot(df, aes(x=Dimension,y=`Relative improvement`,color=`Criterion:`, shape=`Statistic:`)) +
  geom_point() +
  geom_smooth() +
  scale_y_continuous(trans = "log2",breaks = c(1,2,4,8,16,32,64,128,256,512),limits=c(0.75,512)) +
  ylab("") +
  scale_color_manual(values=c(pal[3],pal[5])) +
  theme_bw()
gg5


#
#######
################ fully ill gaussian
#######
#

library(readr)
library(ggplot2)
library(reshape2)

df <- read_table2("inst/output/fullyIllGaussianSequential.txt",
                  col_names = FALSE)
df <- df[,-ncol(df)]
df$X2 <- df$X2 / df$X5
df$X3 <- df$X3 / df$X6
df$X4 <- df$X4 / df$X7
df <- df[,1:4]
colnames(df) <- c("Dimension","Mean","Min","Slowdown")
df <- melt(df, measure.vars=2:3,value.name = "ESS",variable.name = "Statistic:")
df$ESSs <- df$ESS / df$Slowdown
df <- df[,-2]
df <- melt(df, measure.vars=3:4,value.name = "Relative improvement",variable.name = "Criterion:")
df$Statistic <- factor(df$Statistic)
df$Criterion <- factor(df$Criterion)


gg3 <- ggplot(df, aes(x=Dimension,y=`Relative improvement`,color=`Criterion:`, shape=`Statistic:`)) +
  geom_point() +
  geom_smooth() +
  scale_y_continuous(trans = "log2",breaks = c(1,2,4,8,16,32),limits=c(1,32)) +
  ylab("") +
  xlab(NULL) +
  ggtitle("PC algorithms / full \"ill\" target") +
  scale_color_manual(values=c(pal[3],pal[5])) +
  theme_bw()
gg3

df <- read_table2("inst/output/fullyIllGaussianSequential.txt",
                  col_names = FALSE)
df <- df[,-ncol(df)]
df$X2 <- df$X2 / df$X8
df$X3 <- df$X3 / df$X9
df$X4 <- df$X4 / df$X10
df <- df[,1:4]
colnames(df) <- c("Dimension","Mean","Min","Slowdown")
df <- melt(df, measure.vars=2:3,value.name = "ESS",variable.name = "Statistic:")
df$ESSs <- df$ESS / df$Slowdown
df <- df[,-2]
df <- melt(df, measure.vars=3:4,value.name = "Relative improvement",variable.name = "Criterion:")
df$Statistic <- factor(df$Statistic)
df$Criterion <- factor(df$Criterion)


gg6 <- ggplot(df, aes(x=Dimension,y=`Relative improvement`,color=`Criterion:`, shape=`Statistic:`)) +
  geom_point() +
  geom_smooth() +
  scale_y_continuous(trans = "log2",breaks = c(1,2,4,8,16,32,64,128,256,512),limits=c(0.75,512)) +
  ylab("") +
  scale_color_manual(values=c(pal[3],pal[5])) +
  theme_bw()
gg6


library(grid)
library(gridExtra)

source("inst/grid_arrange.R")
ggsave(grid_arrange_shared_legend(gg,gg2,gg3,gg4,gg5,gg6,ncol = 3,nrow=2),
       file="inst/figures/sphereNormFigOrig.pdf",
       width = 12,height = 7)

system2(command = "pdftk",
        args    = c("~/simplicialSampler/inst/figures/sphereNormFigOrig.pdf",
                    "cat 2-end",
                    "output ~/simplicialSampler/inst/figures/rwmComparisonFig.pdf")
)
system2(command = "pdfcrop",
        args    = c("~/simplicialSampler/inst/figures/rwmComparisonFig.pdf",
                    "~/simplicialSampler/inst/figures/rwmComparisonFig.pdf")
)
system2(command = "rm",
        args    = c("~/simplicialSampler/inst/figures/sphereNormFigOrig.pdf")
)

