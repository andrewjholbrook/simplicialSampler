setwd("~/simplicialSampler/")

library(readr)
library(ggplot2)
library(wesanderson)
pal <- wes_palette("FantasticFox1", 5, type = "discrete")

lambdaSearch <- read_table2("inst/output/mtmSigmaSearch.txt", 
                              col_names = FALSE)
lambdaSearch <- lambdaSearch[,c(1,2,3,4,5)]
colnames(lambdaSearch) <- c("Dimension","sigma","meanESS",
                            "Acceptance","TargetAcceptance")
#
###
####### spherical gaussian unscaled
###
#
df1 <- lambdaSearch

keep <- vector() # get max mean ESSs
for (i in unique(df1$Dimension)) {
  keep <- c(keep,max(df1$meanESS[df1$Dimension==i]))
}
df1$isMaxMeanESS <- df1$meanESS %in% keep

df1 <- df1[df1$isMaxMeanESS,]

# for (i in 1:dim(df1)[1]) {
#   if(df1$isMaxMinESS[i] == TRUE &
#      df1$isMaxMeanESS[i] == TRUE) {
#     dummyRow <- df1[i,]
#     dummyRow$isMaxMinESS <- FALSE
#     df1 <- rbind(df1, dummyRow)
#   }
# }
 df1
#saveRDS(df1,file="inst/output/optimalAcceptanceMTM.rds")

gg <- ggplot(data = df1, aes(x=Dimension,y=Acceptance)) +
  geom_jitter() +
  stat_smooth(se=TRUE) +
  theme_bw() +
  ylab("Target acceptance") +
  ggtitle("Emprically optimal acceptance rates, multiple-try Metropolis")

gg

ggsave("mtmScaling.pdf",device="pdf",path="inst/figures/",width = 6,height=3)

system2(command = "pdfcrop",
        args    = c("~/simplicialSampler/inst/figures/mtmScaling.pdf",
                    "~/simplicialSampler/inst/figures/mtmScaling.pdf")
)

