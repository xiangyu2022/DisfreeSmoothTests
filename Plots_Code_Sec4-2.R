#########################  Histogram #####################  
source("C:\\Users\\zhan6\\Downloads\\hjs.R")
library(nloptr)
library(LaplacesDemon)
library(expm)
library(truncnorm)
library(orthopolynom)
library(pracma)
library(tolerance)
load("Sec4-2.RData")

# Orginal data generation
set.seed(1)
n=100; p1=0.3; p2=0.5
a=-10;b=10 #truncated
indx = rmultinom(n=1,size=n,prob=c(1-p1-p2, p1,p2))
xx_q= c(runif(indx[1],a,b), rtruncnorm(indx[2],a,b,mean=-5,sd=3),
          rtrunc(indx[3],spec="laplace",location=5,scale=3,a=a,b=b))
hist(xx_q,breaks=40,prob=T,ylim=c(0,0.2), main="",
     cex=1.5,cex.axis=1.5,lwd=4.5,cex.lab=1.5,xlim=c(-10,10),xlab=c("x"))
M1 = 6
densq = function(x) (1-p1-p2)*dunif(x,a,b)+p1*dtruncnorm(x,a,b,mean=-5,sd=3)+
  p2*dtrunc(x,spec="laplace",location=5,scale=3,a=a,b=b)
curve(densq(x),xlim=c(a,b),ylim=c(0,1),add=T,col="red",lty=5,lwd=4.5)

# Reference distribution is a truncated normal distribution with unknown mean and sd
densf <- function(x,pars) dtruncnorm(x,a=a,b=b,mean=pars[1],sd=pars[2])
likef <- function(pars) -sum(log(densf(xx_q,pars)))
MLEf <- constrOptim(theta=c(mean(xx_q),sd(xx_q)),f=likef,ui=c(0,1),ci=c(0),method = "Nelder-Mead")$par
curve(densf(x,MLEf),xlim=c(a,b),ylim=c(0,1),add=T,col="blue3",lty=1,lwd=4.5)

# Postulated model g1 is a mixture of half unif[a,b], half truncated normal with unknown mean and sd
densg1 <- function(x,pars) 0.5*dunif(x,a,b) + 0.5*dtruncnorm(x,a=a,b=b,mean=pars[1],sd=pars[2])
likeg1 <- function(pars) -sum(log(densg1(xx_q,pars)))
MLEg1 <- constrOptim(theta=c(mean(xx_q),sd(xx_q)),f=likeg1,ui=c(0,1),ci=c(0),method = "Nelder-Mead")$par
curve(densg1(x,MLEg1),xlim=c(a,b),ylim=c(0,1),col="orange",add=T,lty=2,lwd=4.5)
G_g1 <- function(x) integrate(function(x)densg1(x,MLEg1),a,x)$value
G_g1 <- Vectorize(G_g1)

# Postulated model g2 is a Laplace distribution with unknown location and scale
densg2 <- function(x,pars) 0.3*dtruncnorm(x,a=a,b=b,mean=-5,sd=1) + 0.7*dtruncnorm(x,a=a,b=b,mean=pars[1],sd=pars[2]) 
likeg2 <- function(pars) -sum(log(densg2(xx_q,pars)))
MLEg2 <- constrOptim(theta=c(mean(xx_q),sd(xx_q)),f=likeg2,ui=c(0,1),ci=c(0),method = "Nelder-Mead")$par
curve(densg2(x,MLEg2),xlim=c(a,b),ylim=c(0,1),col="grey10",add=T,lty=3,lwd=4.5)
G_g2 <- function(x) integrate(function(x)densg2(x,MLEg2),a,x)$value
G_g2 <- Vectorize(G_g2)

# g3 is a mixture of unif[a,b], truncated normal and truncated laplace, with the wrong parameters and unknown proportions 
densg3 <- function(x,pars) (1-pars[1]-pars[2])*dunif(x,a,b) + pars[1]*dtruncnorm(x,a,b,mean=-4,sd=1) + pars[2]*dtrunc(x,spec="laplace",location=4,scale=1,a=a,b=b)
likeg3 <- function(pars) -sum(log(densg3(xx_q,pars)))
MLEg3 <- constrOptim(theta=c(0.3,0.5), f=likeg3, ui= rbind(diag(2),-diag(2),c(1,1),c(-1,-1)),
                     ci=c(rep(0,2),rep(-1,2),0,-1), method = "Nelder-Mead", control = list(abstol=1e-15,reltol=1e-15))$par
curve(densg3(x,MLEg3),xlim=c(a,b),ylim=c(0,1),col="orchid3",add=T,cex=1.5,lty=4,lwd=4.5)
G_g3 <- function(x) integrate(function(x)densg3(x,MLEg3),a,x)$value
G_g3 <- Vectorize(G_g3)


#### This generates the Figure 2 in Section 4.2
legend("topleft", legend=c("Q", expression(F[gamma]), expression(G[beta * ',' * 1]), expression(G[beta * ',' * 2]), expression(G[beta * ',' * 3])), lty=c(5,1:4), col= c("red","blue3","orange","grey","orchid3"), lwd=4.5, cex=1.3, bty = "n")

setEPS()
postscript("hist.eps",width = 9, height = 6)
# par(mar = c(5, 6.5, 2, 2))
hist(xx_q,breaks=40,prob=T,ylim=c(0,0.16), main="",cex=1.5,cex.axis=1.5,lwd=4.5,cex.lab=1.5,xlim=c(-10,10),xlab=c("x"))
curve(densq(x),xlim=c(a,b),ylim=c(0,1),add=T,col="red",lty=5,lwd=4.5)
curve(densf(x,MLEf),xlim=c(a,b),ylim=c(0,1),add=T,col="blue3",lty=1,lwd=4.5)
curve(densg1(x,MLEg1),xlim=c(a,b),ylim=c(0,1),col="orange",add=T,lty=2,lwd=4.5)
curve(densg2(x,MLEg2),xlim=c(a,b),ylim=c(0,1),col="grey10",add=T,lty=3,lwd=4.5)
curve(densg3(x,MLEg3),xlim=c(a,b),ylim=c(0,1),col="orchid3",add=T,cex=1.5,lty=4,lwd=4.5)
legend("topleft", legend=c("Q", expression(F[gamma]), expression(G[beta * ',' * 1]), expression(G[beta * ',' * 2]), expression(G[beta * ',' * 3])), lty=c(5,1:4), col= c("red","blue3","orange","grey10","orchid3"), lwd=4.5, cex=1.3, bty = "n")
dev.off()


#### This generates the upper left panel of Figure 3 in Section 4.2
cc=seq(0,8,length=1000)
prob_tn_boot_f <- prob_tn_boot_g1 <- prob_tn_boot_g2 <- prob_tn_boot_g3 <-numeric(1000)
prob_tn_order_boot_f <- prob_tn_order_boot_g1 <- prob_tn_order_boot_g2 <- prob_tn_order_boot_g3 <-numeric(1000)
for (k in 1:1000){
  # f
  prob_tn_order_boot_f[k] <- sum(res_tn_order_boot_f<cc[k])/length(res_tn_order_boot_f)
  # g1
  prob_tn_boot_g1[k] <- sum(res_tn_boot_g1<cc[k])/length(res_tn_boot_g1)
  prob_tn_order_boot_g1[k] <- sum(res_tn_order_boot_g1<cc[k])/length(res_tn_order_boot_f)
  # g2
  prob_tn_boot_g2[k] <- sum(res_tn_boot_g2<cc[k])/length(res_tn_boot_g2)
  prob_tn_order_boot_g2[k] <- sum(res_tn_order_boot_g2<cc[k])/length(res_tn_order_boot_f)
  # g3
  prob_tn_boot_g3[k] <- sum(res_tn_boot_g3<cc[k])/length(res_tn_boot_g3)
  prob_tn_order_boot_g3[k] <- sum(res_tn_order_boot_g3<cc[k])/length(res_tn_order_boot_f)
}
postscript("tnbefore.eps",width = 9, height = 6)
par(mar = c(5, 6.5, 2, 2))
plot(cc,prob_tn_boot_f,type = "l", lty=1,col="blue3",xlab="c",cex=2.5,cex.axis=2.5,lwd=4,cex.lab=2.3,
     ylab=expression('P('*hat(T)["n"]<='c)'))
lines(cc,prob_tn_boot_g1,type = "l",col="orange",lty=3,lwd=4)
lines(cc,prob_tn_boot_g2,type = "l",col="grey10",lty=2,lwd=4)
lines(cc,prob_tn_boot_g3,type = "l",col="orchid3",lty=4,lwd=4)
legend("bottomright", legend=c(expression(F[gamma]), expression(G[beta * ',' * 1]), expression(G[beta * ',' * 2]), expression(G[beta * ',' * 3])),lty=1:4,col  = c("blue3","orange","grey10","orchid3"), lwd=4,cex=2,bty = "n")
dev.off()

#### This generates the upper right panel of Figure 3 in Section 4.2
cc=seq(0,20,length=1000)
for (k in 1:1000){
  # q
  prob_tn_order_boot_f[k] <- sum(res_tn_order_boot_f<cc[k])/length(res_tn_order_boot_f)
  # g1
  prob_tn_order_boot_g1[k] <- sum(res_tn_order_boot_g1<cc[k])/length(res_tn_order_boot_f)
  # g2
  prob_tn_order_boot_g2[k] <- sum(res_tn_order_boot_g2<cc[k])/length(res_tn_order_boot_f)
  # g3
  prob_tn_order_boot_g3[k] <- sum(res_tn_order_boot_g3<cc[k])/length(res_tn_order_boot_f)
}
postscript("subtnbefore.eps",width = 9, height = 6)
par(mar = c(5, 6.5, 2, 2))
plot(cc,prob_tn_order_boot_f,type = "l", lty=1,col="blue3",xlab="c",cex=2.5,cex.axis=2.5,lwd=4,cex.lab=2.3,
     ylab=expression('P('*tilde(T)["n"]<='c)'))
lines(cc,prob_tn_order_boot_g1,type = "l",col="orange",lty=3,lwd=4)
lines(cc,prob_tn_order_boot_g2,type = "l",col="grey10",lty=2,lwd=4)
lines(cc,prob_tn_order_boot_g3,type = "l",col="orchid3",lty=4,lwd=4)
legend("bottomright", legend=c(expression(F[gamma]), expression(G[beta * ',' * 1]), expression(G[beta * ',' * 2]), expression(G[beta * ',' * 3])),lty=1:4,col  = c("blue3","orange","grey10","orchid3"), lwd=4,cex=2,bty = "n")
dev.off()


#### This generates the lower left panel of Figure 3 in Section 4.2
cc=seq(0,4,length=1000)
prob_tn_k2_f <- prob_tn_k2_g1 <- prob_tn_k2_g2 <- prob_tn_k2_g3<- numeric(1000)
for (k in 1:1000){
  prob_tn_k2_f[k] <- sum(res_tn_boot_f<cc[k])/length(res_tn_boot_f)
  prob_tn_k2_g1[k] <- sum(res_tn_boot_K2_g1<cc[k])/length(res_tn_boot_K2_g1)
  prob_tn_k2_g2[k] <- sum(res_tn_boot_K2_g2<cc[k])/length(res_tn_boot_K2_g2)
  prob_tn_k2_g3[k] <- sum(res_tn_boot_K2_g3<cc[k])/length(res_tn_boot_K2_g3)
}
postscript("tnafter.eps",width = 9, height = 6)
par(mar = c(5, 6.5, 2, 2))
plot(cc,prob_tn_k2_f,type = "l", lty=1,col="blue3",xlab="c",cex=2.5,cex.axis=2.5,lwd=4,cex.lab=2.3,
     ylab=expression('P('*hat(T)["n"]^{scriptscriptstyle(K)}<='c)'))
lines(cc,prob_tn_k2_g1,type = "l",col="orange",lty=3,lwd=4)
lines(cc,prob_tn_k2_g2,type = "l",col="grey10",lty=2,lwd=4)
lines(cc,prob_tn_k2_g3,type = "l",col="orchid3",lty=4,lwd=4)
legend("bottomright", legend=c(expression(F[gamma]), expression(G[beta * ',' * 1]), expression(G[beta * ',' * 2]), expression(G[beta * ',' * 3])),lty=1:4,col  = c("blue3","orange","grey10","orchid3"), lwd=4,cex=2,bty = "n")
dev.off()

#### This generates the lower right panel of Figure 3 in Section 4.2
cc=seq(0,12,length=1000)
prob_tn_order_k2_f <- prob_tn_order_k2_g1 <- prob_tn_order_k2_g2 <- prob_tn_order_k2_g3<- numeric(1000)
for (k in 1:1000){
  prob_tn_order_k2_f[k] <- sum(res_tn_order_boot_f<cc[k])/length(res_tn_order_boot_f)
  prob_tn_order_k2_g1[k] <- sum(res_tn_order_boot_K2_g1<cc[k])/length(res_tn_order_boot_K2_g1)
  prob_tn_order_k2_g2[k] <- sum(res_tn_order_boot_K2_g2<cc[k])/length(res_tn_order_boot_K2_g2)
  prob_tn_order_k2_g3[k] <- sum(res_tn_order_boot_K2_g3<cc[k])/length(res_tn_order_boot_K2_g3)
}
postscript("subtnafter.eps",width = 9, height = 6)
par(mar = c(5, 6.5, 2, 2))
plot(cc,prob_tn_order_k2_f,type = "l", lty=1,col="blue3",xlab="c",cex=2.5,cex.axis=2.5,lwd=4,
     cex.lab=2.3,ylab=expression('P('*tilde(T)["n"]^{scriptscriptstyle(K)}<='c)'))
lines(cc,prob_tn_order_k2_g1,type = "l",col="orange",lty=3,lwd=4)
lines(cc,prob_tn_order_k2_g2,type = "l",col="grey10",lty=2,lwd=4)
lines(cc,prob_tn_order_k2_g3,type = "l",col="orchid3",lty=4,lwd=4)
legend("bottomright", legend=c(expression(F[gamma]), expression(G[beta * ',' * 1]), expression(G[beta * ',' * 2]), expression(G[beta * ',' * 3])),lty=1:4,col  = c("blue3","orange","grey10","orchid3"), lwd=4,cex=2,bty = "n")
dev.off()



