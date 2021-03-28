setwd("~/simplicialSampler/")

library(ggplot2)
library(readr)
library(wesanderson)
pal <- wes_palette("FantasticFox1", 5, type = "discrete")


df <- read_table2("inst/output/bimodalComparison.txt", col_names = FALSE)
df <- df[,-10]
df$X3 <- df$X2 / df$X3
df$X5 <- df$X4 / df$X5
df$X7 <- df$X6 / df$X7
df$X9 <- df$X8 / df$X9
df <- df[df$X1<5,]

df1 <- data.frame(df$X1,df$X2,df$X4,df$X6,df$X8)
df2 <- data.frame(df$X1,df$X3,df$X5,df$X7,df$X9)
colnames(df1) <- c("Dimension","sSS","sRWM","uRWM","uSS")
colnames(df2) <- c("Dimension","sSS","sRWM","uRWM","uSS")
df1 <- reshape2::melt(df1,measure.vars=2:5,
                      value.name="Jumps",variable.name="Algorithm")
df1$Algorithm <- factor(df1$Algorithm,levels = c("sSS","uSS","sRWM","uRWM"))
df1$Dimension <- factor(df1$Dimension)

df2 <- reshape2::melt(df2,measure.vars=2:5,
                      value.name="Jumps",variable.name="Algorithm")
df2$Algorithm <- factor(df2$Algorithm,levels = c("sSS","uSS","sRWM","uRWM"))
df2$Dimension <- factor(df2$Dimension)

gg <- ggplot(df1,aes(x=Dimension,y=Jumps,fill=Algorithm)) +
  geom_boxplot() +
  scale_fill_manual(values=c(pal[2],pal[3],pal[4],pal[5]))+
  ylab("Intermodal jumps") +
  theme_bw()
gg


gg2 <- ggplot(df2,aes(x=Dimension,y=Jumps,fill=Algorithm)) +
  geom_boxplot() +
  ylab("Intermodal jumps per second") +
  scale_fill_manual(values=c(pal[2],pal[3],pal[4],pal[5]))+
  theme_bw()
gg2

library(grid)
library(gridExtra)

source("inst/grid_arrange.R")
ggsave(grid_arrange_shared_legend(gg,gg2,ncol = 2,position = "right"),
       file="inst/figures/sphereNormFigOrig.pdf",
       width = 10,height = 4)

system2(command = "pdftk",
        args    = c("~/simplicialSampler/inst/figures/sphereNormFigOrig.pdf",
                    "cat 2-end",
                    "output ~/simplicialSampler/inst/figures/multiModFig2.pdf")
)
system2(command = "pdfcrop",
        args    = c("~/simplicialSampler/inst/figures/multiModFig2.pdf",
                    "~/simplicialSampler/inst/figures/multiModFig2.pdf")
)
system2(command = "rm",
        args    = c("~/simplicialSampler/inst/figures/sphereNormFigOrig.pdf")
)









