library(mvtnorm)
library(coda)
library(ggplot2)
library(cowplot)

setwd("~/simplicialSampler/")

source("R/simplicialSampler.R")

set.seed(1)

N <- 1000
Dim <- 2
distrib <- "sphericalGaussian"
x <- rep(-10,Dim)
x <- c(x,rep(0,N-Dim))
  
v_Nplus1 <- rep((1+sqrt(N+1))/N,N)
M <- rbind(diag(N), v_Nplus1)
M <- M - matrix(v_Nplus1,N+1,N)
M <- M * 60
U <- pracma::randortho(n=N)
M2 <- M %*% U
M <- M + matrix(x,N+1,N,byrow = TRUE)
M2 <- M2 + matrix(x,N+1,N,byrow = TRUE)
M2 <- M2[,1:Dim]

M3 <- M2[sample(1:(N+1),1,prob = target(M2,distrib)),]
M <- data.frame(x=M[,1],y=M[,2])
M2 <- data.frame(x=M2[,1],y=M2[,2])
M3 <- data.frame(x=M3[1],y=M3[2])


x <- matrix(rnorm(100000*2),100000,2)
M4 <- data.frame(x=x[,1],y=x[,2])


gg <- ggplot(M,aes(x=x,y=y)) +
  geom_point(data=M2,aes(x=x,y=y),color="lightblue") +
  geom_point() +
  geom_point(data=M3,color="red") +
  geom_density2d(data=M4,color="plum2",adjust=2) +
  geom_segment(aes(x=-10,y=-10,xend=M[1,2],yend=M[1,1]),color="black",linetype="dashed")+
  geom_segment(aes(x=-10,y=-10,xend=M[2,2],yend=M[2,1]),color="black",linetype="dashed")+
  geom_segment(aes(x=-10,y=-10,xend=M[3,2],yend=M[3,1]),color="black",linetype="dashed")+
  geom_segment(aes(x=M[1,2],y=M[1,1],xend=M[3,2],yend=M[3,1]),color="black",linetype="dashed")+
  geom_segment(aes(x=M[2,2],y=M[2,1],xend=M[3,2],yend=M[3,1]),color="black",linetype="dashed")+
  geom_segment(aes(x=M[2,2],y=M[2,1],xend=M[1,2],yend=M[1,1]),color="black",linetype="dashed")+
  geom_point(data = M[N+1,],color="tan1") +
  xlim(c(-20,60)) +
  ylim(c(-20,60)) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm")) 
gg <- as_grob(gg)
  
x <- c(M3[,1],M3[,2])
x <- c(x,rep(0,N-Dim))

v_Nplus1 <- rep((1+sqrt(N+1))/N,N)
M <- rbind(diag(N), v_Nplus1)
M <- M - matrix(v_Nplus1,N+1,N)
M <- M *60
U <- pracma::randortho(n=N)
M2 <- M %*% U
M <- M + matrix(x,N+1,N,byrow = TRUE)
M2 <- M2 + matrix(x,N+1,N,byrow = TRUE)
M2 <- M2[,1:Dim]

M3 <- M2[sample(1:(N+1),1,prob = target(M2,distrib)),]
M <- data.frame(x=M[,1],y=M[,2])
M2 <- data.frame(x=M2[,1],y=M2[,2])
M3 <- data.frame(x=M3[1],y=M3[2])


M4 <- matrix(rnorm(100000*2),100000,2)
M4 <- data.frame(x=M4[,1],y=M4[,2])


gg2 <- ggplot(M,aes(x=x,y=y)) +
  geom_point(data=M2,aes(x=x,y=y),color="lightblue") +
  geom_point() +
  geom_point(data=M3,color="red") +
  geom_density2d(data=M4,color="plum2",adjust=2) +
  geom_segment(aes(x=M[N+1,1],y=M[N+1,2],xend=M[1,1],yend=M[1,2]),color="black",linetype="dashed")+
  geom_segment(aes(x=M[N+1,1],y=M[N+1,2],xend=M[2,1],yend=M[2,2]),color="black",linetype="dashed")+
  geom_segment(aes(x=M[N+1,1],y=M[N+1,2],xend=M[3,1],yend=M[3,2]),color="black",linetype="dashed")+
  geom_segment(aes(x=M[1,1],y=M[1,2],xend=M[3,1],yend=M[3,2]),color="black",linetype="dashed")+
  geom_segment(aes(x=M[2,1],y=M[2,2],xend=M[3,1],yend=M[3,2]),color="black",linetype="dashed")+
  geom_segment(aes(x=M[2,1],y=M[2,2],xend=M[1,1],yend=M[1,2]),color="black",linetype="dashed")+
  geom_point(data = M[N+1,],color="tan1") +
  xlim(c(-20,60)) +
  ylim(c(-20,60)) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm")) 
gg2 <- as_grob(gg2)


x <- c(M3[,1],M3[,2])
x <- c(x,rep(0,N-Dim))

v_Nplus1 <- rep((1+sqrt(N+1))/N,N)
M <- rbind(diag(N), v_Nplus1)
M <- M - matrix(v_Nplus1,N+1,N)
M <- M *60
U <- pracma::randortho(n=N)
M2 <- M %*% U
M <- M + matrix(x,N+1,N,byrow = TRUE)
M2 <- M2 + matrix(x,N+1,N,byrow = TRUE)
M2 <- M2[,1:Dim]

M3 <- M2[sample(1:(N+1),1,prob = target(M2,distrib)),]
M <- data.frame(x=M[,1],y=M[,2])
M2 <- data.frame(x=M2[,1],y=M2[,2])
M3 <- data.frame(x=M3[1],y=M3[2])


M4 <- matrix(rnorm(100000*2),100000,2)
M4 <- data.frame(x=M4[,1],y=M4[,2])


gg3 <- ggplot(M,aes(x=x,y=y)) +
  geom_point(data=M2,aes(x=x,y=y),color="lightblue") +
  geom_point() +
  geom_density2d(data=M4,color="plum2",adjust=2) +
  geom_point(data=M3,color="red") +
  geom_segment(aes(x=M[N+1,1],y=M[N+1,2],xend=M[1,1],yend=M[1,2]),color="black",linetype="dashed")+
  geom_segment(aes(x=M[N+1,1],y=M[N+1,2],xend=M[2,1],yend=M[2,2]),color="black",linetype="dashed")+
  geom_segment(aes(x=M[N+1,1],y=M[N+1,2],xend=M[3,1],yend=M[3,2]),color="black",linetype="dashed")+
  geom_segment(aes(x=M[1,1],y=M[1,2],xend=M[3,1],yend=M[3,2]),color="black",linetype="dashed")+
  geom_segment(aes(x=M[2,1],y=M[2,2],xend=M[3,1],yend=M[3,2]),color="black",linetype="dashed")+
  geom_segment(aes(x=M[2,1],y=M[2,2],xend=M[1,1],yend=M[1,2]),color="black",linetype="dashed")+
  geom_point(data = M[N+1,],color="tan1") +
  xlim(c(-20,60)) +
  ylim(c(-20,60)) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm")) 
gg3 <- as_grob(gg3)


library(grid)
library(gridExtra)

ggsave(file="inst/figures/hyperSimplex.pdf",plot=arrangeGrob(gg,gg2,gg3,ncol=3,widths = ),
       width = 12,height = 4,device = pdf)
system2(command = "pdfcrop",
        args    = c("~/simplicialSampler/inst/figures/hyperSimplex.pdf",
                    "~/simplicialSampler/inst/figures/hyperSimplex.pdf")
)
