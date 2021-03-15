setwd("~/simplicialSampler/")

library(readr)
library(ggplot2)
library(reshape2)

df <- read_table2("inst/output/bananaComparisons.txt",
                            col_names = FALSE)
df <- df[,-ncol(df)]
df <- df[,-1]
df <- df[,1:6]/df[,7:12]
df <- data.frame( rbind(as.matrix(df[,1:3]),as.matrix(df[,4:6])) )
df$Scaled <- c( rep("Unscaled",100),rep("Scaled",100) )

colnames(df) <- c("Mean","Min","Slowdown","Scaled")
df <- melt(df, measure.vars=1:2,value.name = "ESS",variable.name = "Statistic")
df$ESSs <- df$ESS / df$Slowdown
df <- df[,-1]
df <- melt(df, measure.vars=3:4,value.name = "Relative improvement",variable.name = "Criterion:")
# df$Statistic <- factor(df$`Statistic:`)
# df$Criterion <- factor(df$Criterion)
df$Category <- with(df, interaction(Statistic, Scaled), drop = TRUE )

gg <- ggplot(df, aes(y=`Relative improvement`,x=Category,fill=`Criterion:`)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log2",
                     labels=c(0.125,0.5,2,8,32,128,512),
                     breaks=c(0.125,0.5,2,8,32,128,512),
                     limits=c(0.1,512)) +
  ylab("Relative improvement (uSS/uRWM)") +
  xlab("Criterion") +
  ggtitle("Unscaled algorithms target banana distribution") +
  scale_fill_manual(values=c("maroon1","green3")) +
  theme_bw()
gg


ggsave(gg,file="inst/figures/bananaCompareFig.pdf",
       width = 5,height = 3)


system2(command = "pdfcrop",
        args    = c("~/simplicialSampler/inst/figures/bananaCompareFig.pdf",
                    "~/simplicialSampler/inst/figures/bananaCompareFig.pdf")
)


