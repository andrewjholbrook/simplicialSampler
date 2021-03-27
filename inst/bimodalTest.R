setwd("~/simplicialSampler/")

source("R/simplicialSampler.R")

library(coda)
library(ggplot2)
library(wesanderson)
pal <- wes_palette("FantasticFox1", 5, type = "discrete")


# set seed

n <- 3
N <- 10000

sss <- vector()
srw <- vector()
uss <- vector()
urw <- vector()

set.seed(666)

for(i in 1:3){

results2 <- randomWalk(N=n,x0=rep(0,n), maxIt = N,
                       target = "bimodalGaussian", adaptCov = TRUE)
srw <- c(srw,results2[[1]][,1])


results <- simplicialSampler(N=n,x0=rep(0,n), maxIt = N,lambda=2.4/sqrt(n),
                             adaptStepSize = TRUE, targetAccept = results2[[2]],
                             target = "bimodalGaussian", adaptScales = TRUE,
                             Gaussians = TRUE)
sss <- c(sss,results[[1]][,1])

results3 <- randomWalk(N=n,x0=rep(0,n), maxIt = N,
                       target = "bimodalGaussian", adaptCov = FALSE)
urw <- c(urw,results3[[1]][,1])


results4 <- simplicialSampler(N=n,x0=rep(0,n), maxIt = N,lambda=2.4/sqrt(n),
                             adaptStepSize = TRUE, targetAccept = results3[[2]],
                             target = "bimodalGaussian", adaptScales = FALSE,
                             Gaussians = TRUE)
uss <- c(uss,results4[[1]][,1])
}

df1 <- data.frame(States=srw,Algorithm="sRWM",Chain=rep(1:3,each=N))
df1$States <- df1$States + 10
df2 <- data.frame(States=urw,Algorithm="uRWM",Chain=rep(1:3,each=N))
df3 <- data.frame(States=uss,Algorithm="uSS",Chain=rep(1:3,each=N))
df3$States <- df3$States + 20
df4 <- data.frame(States=sss,Algorithm="sSS",Chain=rep(1:3,each=N))
df4$States <- df4$States + 30


df <- rbind(df4,df3,df2,df1)
df$Algorithm <- factor(df$Algorithm,levels = c("sSS","uSS","sRWM","uRWM"))
df$Iteration <- rep(1:N,6)

gg <- ggplot(df[df$Chain==1 ,],aes(y=States,x=Iteration,color=Algorithm)) +
  geom_line() +
  ylab("Values") +
  ggtitle("Chain 1") +
  scale_color_manual(values=c(pal[2],pal[3],pal[4],pal[5])) +
  scale_y_continuous(breaks=c(0,5,10,15,20,25,30,35), labels=c(0,5,0,5,0,5,0,5)) +
  theme_bw()
gg

gg2 <- ggplot(df[df$Chain==2 ,],aes(y=States,x=Iteration,color=Algorithm)) +
  geom_line() +
  ggtitle("Chain 2") +
  ylab(NULL) +
  scale_color_manual(values=c(pal[2],pal[3],pal[4],pal[5])) +
  scale_y_continuous(breaks=c(0,5,10,15,20,25,30,35), labels=c(0,5,0,5,0,5,0,5)) +
  theme_bw()
gg2

gg3<- ggplot(df[df$Chain==3 ,],aes(y=States,x=Iteration,color=Algorithm)) +
  geom_line() +
  ggtitle("Chain 3") +
  ylab(NULL) +
  scale_color_manual(values=c(pal[2],pal[3],pal[4],pal[5])) +
  scale_y_continuous(breaks=c(0,5,10,15,20,25,30,35), labels=c(0,5,0,5,0,5,0,5)) +
  theme_bw()
gg3

library(grid)
library(gridExtra)

source("inst/grid_arrange.R")
ggsave(grid_arrange_shared_legend(gg,gg2,gg3,ncol = 3),
       file="inst/figures/sphereNormFigOrig.pdf",
       width = 12,height = 4)

system2(command = "pdftk",
        args    = c("~/simplicialSampler/inst/figures/sphereNormFigOrig.pdf",
                    "cat 2-end",
                    "output ~/simplicialSampler/inst/figures/multiModFig.pdf")
)
system2(command = "pdfcrop",
        args    = c("~/simplicialSampler/inst/figures/multiModFig.pdf",
                    "~/simplicialSampler/inst/figures/multiModFig.pdf")
)
system2(command = "rm",
        args    = c("~/simplicialSampler/inst/figures/sphereNormFigOrig.pdf")
)




