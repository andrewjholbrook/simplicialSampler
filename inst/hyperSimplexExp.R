library(mvtnorm)
library(coda)
library(ggplot2)
library(wesanderson)

N <- 3
set.seed(1)
pal <- wes_palette("Zissou1", 5, type = "continuous")


setwd("~/simplicialSampler/")

source("R/simplicialSampler.R")

maxIt <- 10000
N     <- 3

proposal <- function(N,x,lambda,distrib,gaussians=FALSE) {
    
    Dim <- length(x)
    x <- c(x,rep(0,N-Dim))
    
    v_Nplus1 <- rep((1+sqrt(N+1))/N,N)
    M <- rbind(diag(N), v_Nplus1)
    M <- M - matrix(v_Nplus1,N+1,N)
    M4 <- M * lambda /sqrt(2) 
    if(gaussians) M4 <- M4 * sqrt(rchisq(n=1,df=N))
    U <- pracma::randortho(n=N)
    M4 <- M4 %*% U
    M4 <- M4 + matrix(x,N+1,N,byrow = TRUE)
    M4 <- M4[,1:Dim]
    output <- M4[sample(1:(N+1),1,prob = target(M4,distrib)),]
    
    return(output)
}

simplicialSampler <- function(nProps,N,x0,lambda=1,maxIt=10000,adaptStepSize=TRUE,
                              targetAccept=0.5,target=NULL,
                              Gaussians=FALSE) {
    #if(N!=length(x0)) stop("Dimension mismatch.")
    
    chain <- matrix(0,maxIt,N)
    
    Acceptances = 0 # total acceptances within adaptation run (<= SampBound)
    SampBound = 5   # current total samples before adapting radius
    SampCount = 0   # number of samples collected (adapt when = SampBound)
    Proposed = 0
    
    accept <- rep(0,maxIt)
    chain[1,] <- x0
    for (i in 2:maxIt){
        Proposed = Proposed + 1
        
        
        prpsl <- proposal(nProps,chain[i-1,],lambda,target,
                          gaussians = Gaussians)
        
        chain[i,] <- prpsl
        SampCount <- SampCount + 1
        if(any(chain[i,] != chain[i-1,])){
            accept[i] <- 1
            Acceptances = Acceptances + 1
        }  
        
        if (SampCount == SampBound) { # adjust lambda at increasing intervals
            AcceptRatio <- Acceptances / SampBound
            AdaptRatio  <- AcceptRatio / targetAccept
            if (AdaptRatio>2) AdaptRatio <- 2
            if (AdaptRatio<0.5) AdaptRatio <- 0.5
            lambda <- lambda * AdaptRatio
            
            SampCount <- 0
            SampBound <- ceiling(SampBound^1.01)
            Acceptances <- 0
        }
        
        
        if(i %% 1000 == 0) cat(i,"\n")
    }
    ratio <- sum(accept)/(maxIt-1)
    cat("Acceptance ratio: ", ratio,"\n")
    return(list(chain,ratio,lambda))
}

# spherical target / no adapt scales
output <- simplicialSampler(nProps=100,N=N, x0=rep(0,N), maxIt = 100000, adaptStepSize=TRUE,
                            target = "sphericalGaussian",targetAccept = 0.5)
effectiveSize(as.mcmc(output[[1]]))
for(i in 1:N){
    qqplot(output[[1]][,i],rnorm(maxIt,sd=1))
    qqline(output[[1]][,i], distribution = function(p) qnorm(p,sd=1 ))
}

df1 <- data.frame(Samples=c(output[[1]][,1],output[[1]][,2],output[[1]][,3]),
                  Dimension=c(rep(1,100000),rep(2,100000),rep(3,100000)))
df1$Dimension <- factor(df1$Dimension)
gg <- ggplot(df1,aes(sample=Samples,color=Dimension)) +
    geom_abline(slope=1) + stat_qq() + 
    scale_color_manual(values=c(pal[5],pal[3],pal[1])) +
    xlim(c(-5,5)) + ylim(c(-5,5)) +
    ylab("Sample quantiles") + xlab("Theoretical quantiles") + 
    ggtitle("3D Gaussian target: 100 proposals / 0.5 target rate") +
    guides(color = guide_legend(override.aes = list(size = 10))) +
    theme_bw()
gg

###################################################
N <- 3
# get samples directly from mixture of gaussians
DirectSamples <- matrix(0,100000,N)
for(i in 1:100000) {
    DirectSamples[i,] <- rnorm(N) + sample(c(0,5),size=1)
}

scaleFUN <- function(x) sprintf("%.1f", x)

# spherical target / no adapt scales
output <- simplicialSampler(nProps=100,N=N, x0=rep(0,N), maxIt = 100000, adaptStepSize=TRUE,
                            target = "bimodalGaussian",targetAccept = 0.15)

output[[1]][,1] <- output[[1]][order(output[[1]][,1]),1]
output[[1]][,2] <- output[[1]][order(output[[1]][,2]),2]
output[[1]][,3] <- output[[1]][order(output[[1]][,3]),3]
df1 <- data.frame(Samples=c(output[[1]][,1],output[[1]][,2],output[[1]][,3]),
                  Dimension=c(rep(1,100000),rep(2,100000),rep(3,100000)))
df1$Dimension <- factor(df1$Dimension)
DirectSamples[,1] <- DirectSamples[order(DirectSamples[,1]),1]
DirectSamples[,2] <- DirectSamples[order(DirectSamples[,2]),2]
DirectSamples[,3] <- DirectSamples[order(DirectSamples[,3]),3]
df1$Direct <- c(DirectSamples[,1],DirectSamples[,2],DirectSamples[,3])
gg2 <- ggplot(df1,aes(x=Direct,y=Samples,color=Dimension)) +
    geom_abline(slope=1) +
    geom_point() +
    xlim(c(-5,10)) + ylim(c(-5,10)) +
    ylab("Sample quantiles") + xlab("Theoretical quantiles") + 
    ggtitle("3D mixture of Gaussians target: 100 proposals / 0.15 target rate") +
    scale_y_continuous(labels = scaleFUN) +
    scale_color_manual(values=c(pal[5],pal[3],pal[1])) +
    guides(color = guide_legend(override.aes = list(size = 10))) +
    theme_bw()
gg2

library(grid)
library(gridExtra)
source("inst/grid_arrange.R")
ggsave(grid_arrange_shared_legend(gg,gg2,ncol = 2,nrow=1,position = "bottom"),
       file="inst/figures/sphereNormFigOrig.pdf",
       width = 12,height = 4)

system2(command = "pdftk",
        args    = c("~/simplicialSampler/inst/figures/sphereNormFigOrig.pdf",
                    "cat 2-end",
                    "output ~/simplicialSampler/inst/figures/accuracyFig.pdf")
)
system2(command = "pdfcrop",
        args    = c("~/simplicialSampler/inst/figures/accuracyFig.pdf",
                    "~/simplicialSampler/inst/figures/accuracyFig.pdf")
)
system2(command = "rm",
        args    = c("~/simplicialSampler/inst/figures/sphereNormFigOrig.pdf")
)












