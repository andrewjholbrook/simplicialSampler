setwd("~/simplicialSampler/")

library(readr)
library(ggplot2)
library(wesanderson)
pal <- wes_palette("FantasticFox1", 5, type = "discrete")


df <- read_table2("inst/output/p1P2Comparison.txt", 
                                    col_names = FALSE)

df <- df[,1:4]
df$X4 <- 10/11 * df$X4
df$X5 <- df$X3 / df$X4
colnames(df) <- c("Sampler","Dimension", "ESS", "Seconds", "ESSs")
df$Seconds <- df$Seconds /10000

df$Sampler <- factor(df$Sampler, labels = c("Slow P1, 2D props",
                                            "Slow P1, D props",
                                            "Fast P1, 2D props",
                                            "Fast P1, D props",
                                            "Simpl"))

gg <- ggplot(df,aes(x=Dimension,y=ESS,color=Sampler)) +
  geom_smooth(se=FALSE) +
  ylab("Mean ESS") +
  scale_color_manual(values = pal) +
  theme_bw()
gg

gg2 <- ggplot(df,aes(x=Dimension,y=Seconds,color=Sampler)) +
  geom_smooth(se=FALSE) +
  scale_color_manual(values = pal) +
  ylab("Seconds per iteration") +
  theme_bw()
gg2

gg3 <- ggplot(df,aes(x=Dimension,y=ESSs,color=Sampler)) +
  geom_smooth(se=FALSE) +
  scale_color_manual(values = pal) +
  ylab("Mean ESS per second") +
  theme_bw()
gg3

# library(grid)
# library(gridExtra)
# 
# source("inst/grid_arrange.R")

library(ggpubr)

ggarrange(gg, gg2, gg3, ncol = 3, nrow = 1, common.legend = TRUE, 
          legend = "bottom")

ggsave(ggarrange(gg, gg2, gg3, ncol = 3, nrow = 1, common.legend = TRUE, 
                 legend = "bottom"),file="inst/figures/p1P2Fig.pdf",
       device = "pdf",
       width = 12,height = 4)

