setwd("~/simplicialSampler/")

library(readr)
library(ggplot2)
library(reshape2)

# df <- read_table2("inst/output/rwmComparison.txt", 
#                             col_names = FALSE)
# df <- df[,-ncol(df)]
# df <- rbind(as.matrix(df[,1:4]),as.matrix(df[,c(1,5:7)]))
# df <- data.frame(df)
# colnames(df) <- c("Dimension", "Mean","Min","Seconds")
# df$Algorithm <- c(rep("uSS",nrow(df)/2),rep("uRWM",nrow(df)/2))
# df <- melt(df,measure.vars=2:3,value.name = "ESS",variable.name = "Criterion")
# df$ESSs <- df$ESS / df$Seconds
# df$Criterion <- factor(df$Criterion)
# df$Algorithm <- factor(df$Algorithm)
# 
# gg <- ggplot(df, aes(x=Dimension,y=ESS,color=Algorithm, shape=Criterion)) +
#   geom_jitter() +
#   geom_line() +
#   ylab("Effective sample size") +
#   ggtitle("Unscaled methods for spherical Gaussian target") +
#   theme_bw()
# gg
# 
# gg2 <- ggplot(df, aes(x=Dimension,y=ESSs,color=Algorithm, shape=Criterion)) +
#   geom_jitter() +
#   geom_line() +
#   ylab("Effective sample size per second") +
#   theme_bw()
# gg2
# 
# library(grid)
# library(gridExtra)
# 
# source("inst/grid_arrange.R")
# ggsave(grid_arrange_shared_legend(gg,gg2),file="inst/figures/sphereNormFigOrig.pdf",
#        width = 9,height = 4)
# 
# system2(command = "pdftk", 
#         args    = c("~/simplicialSampler/inst/figures/sphereNormFigOrig.pdf",
#                     "cat 2-end",
#                     "output ~/simplicialSampler/inst/figures/sphereNormFig.pdf") 
# )
# system2(command = "pdfcrop", 
#         args    = c("~/simplicialSampler/inst/figures/sphereNormFig.pdf", 
#                     "~/simplicialSampler/inst/figures/sphereNormFig.pdf") 
# )
# system2(command = "rm", 
#         args    = c("~/simplicialSampler/inst/figures/sphereNormFigOrig.pdf") 
# )


setwd("~/simplicialSampler/")

library(readr)
library(ggplot2)
library(reshape2)

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
  geom_jitter() +
  geom_line() +
  scale_y_continuous(trans = "log2",breaks = c(1,2,4,8,16,32),limits=c(1,32)) +
  ylab("Relative improvement (SS/RWM)") +
  ggtitle("Unscaled algorithms target spherical Gaussian") +
  theme_bw()
gg


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
  geom_jitter() +
  geom_line() +
  scale_y_continuous(trans = "log2",breaks = c(1,2,4,8,16,32),limits=c(1,32)) +
  ylab("") +
  ggtitle("Scaled algorithms target \"ill\" Gaussian") +
  theme_bw()
gg2


library(grid)
library(gridExtra)

source("inst/grid_arrange.R")
ggsave(grid_arrange_shared_legend(gg,gg2),file="inst/figures/sphereNormFigOrig.pdf",
       width = 9,height = 4)

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

