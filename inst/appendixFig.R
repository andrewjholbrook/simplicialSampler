library(mvtnorm)
library(coda)
library(ggplot2)
library(wesanderson)

N <- 3
set.seed(1)
pal <- wes_palette("Zissou1", 5, type = "continuous")


setwd("~/simplicialSampler/")

source("R/simplicialSampler.R")

N     <- 3

output <- simplicialSampler(N=N, x0=rep(0,N), maxIt = 110000, adaptStepSize=TRUE,
                            target = "sphericalGaussian",targetAccept = 0.5,burnin = 10000)


df1 <- data.frame(Samples=c(output[[1]][,1],output[[1]][,2],output[[1]][,3]),
                  Dimension=c(rep(1,100000),rep(2,100000),rep(3,100000)))
df1$Dimension <- factor(df1$Dimension)
gg <- ggplot(df1,aes(sample=Samples,color=Dimension)) +
  geom_abline(slope=1) + stat_qq() + 
  scale_color_manual(values=c(pal[5],pal[3],pal[1])) +
  xlim(c(-5,5)) + ylim(c(-5,5)) +
  ylab("Simplicial sampler quantiles") + xlab("Theoretical quantiles") + 
  ggtitle("3D spherical Gaussian target") +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  theme_bw()
gg

N <- 3
# get samples directly from mixture of gaussians
DirectSamples <- matrix(0,100000,N)
for(i in 1:100000) {
  DirectSamples[i,] <- diag(sqrt(1:3))%*%rnorm(N) 
}

scaleFUN <- function(x) sprintf("%.1f", x)

output <- simplicialSampler(N=N, x0=rep(0,N),
                            maxIt = 110000,
                            adaptStepSize=TRUE,
                            adaptScales = TRUE,
                            target = "diagGaussian",
                            targetAccept = 0.5,
                            burnin=10000)

output[[1]][,1] <- output[[1]][order(output[[1]][,1]),1]
output[[1]][,2] <- output[[1]][order(output[[1]][,2]),2]
output[[1]][,3] <- output[[1]][order(output[[1]][,3]),3]
df1 <- data.frame(Samples=c(output[[1]][,3],output[[1]][,2],output[[1]][,1]),
                  Dimension=c(rep(3,100000),rep(2,100000),rep(1,100000)))
df1$Dimension <- factor(df1$Dimension)
#df1$Dimension <- forcats::fct_rev(df1$Dimension)
DirectSamples[,1] <- DirectSamples[order(DirectSamples[,1]),1]
DirectSamples[,2] <- DirectSamples[order(DirectSamples[,2]),2]
DirectSamples[,3] <- DirectSamples[order(DirectSamples[,3]),3]
df1$Direct <- c(DirectSamples[,3],DirectSamples[,2],DirectSamples[,1])
gg2 <- ggplot(df1,aes(x=Direct,y=Samples,color=Dimension)) +
  geom_abline(slope=1) +
  geom_point() +
 # xlim(c(-5,10)) + ylim(c(-5,10)) +
  ylab("Simplicial sampler quantiles") + xlab("Direct sample quantiles") + 
  ggtitle("3D non-spherical Gaussian target") +
  scale_y_continuous(labels = scaleFUN) +
  scale_color_manual(values=c(pal[5],pal[3],pal[1])) +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  theme_bw()
gg2



library(ggpubr)

ggsave(file="~/simplicialSampler/inst/figures/accuracyFig.pdf",
       ggarrange(gg,gg2, ncol = 2, nrow = 1, common.legend = TRUE, 
                 legend = "bottom"),device = "pdf",width = 12,height=4)

system2(command = "pdfcrop",
        args    = c("~/simplicialSampler/inst/figures/accuracyFig.pdf",
                    "~/simplicialSampler/inst/figures/accuracyFig.pdf")
)














