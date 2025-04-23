#########################  Histogram #####################  
source("hjs.R")
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
xx_k2 = c(runif(indx[1],a,b), rtruncnorm(indx[2],a,b,mean=-5,sd=3),
          rtrunc(indx[3],spec="laplace",location=5,scale=3,a=a,b=b))
hist(xx_k2,breaks=40,prob=T,ylim=c(0,0.2), main="",
     cex=1.5,cex.axis=1.5,lwd=4.5,cex.lab=1.5,xlim=c(-10,10),xlab=c("x"))
M1 = 6

densp = function(x) (1-p1-p2)*dunif(x,a,b)+p1*dtruncnorm(x,a,b,mean=-5,sd=3)+
  p2*dtrunc(x,spec="laplace",location=5,scale=3,a=a,b=b)

curve(densp(x),xlim=c(a,b),ylim=c(0,1),add=T,col="red",lty=5,lwd=4.5)

# Reference distribution is a truncated normal distribution with unknown mean and sd
densq <- function(x,pars) dtruncnorm(x,a=a,b=b,mean=pars[1],sd=pars[2])
likeq <- function(pars) -sum(log(densq(xx_k2,pars)))
MLEq <- constrOptim(theta=c(mean(xx_k2),sd(xx_k2)),f=likeq,ui=c(0,1),ci=c(0),method = "Nelder-Mead")$par
curve(densq(x,MLEq),xlim=c(a,b),ylim=c(0,1),add=T,col="blue3",lty=1,lwd=4.5)
G_obsq = function(x) ptruncnorm(x,a=a,b=b,mean=MLEq[1],sd=MLEq[2])

# Postulated model F1 is a mixture of half unif[a,b], half truncated normal with unknown mean and sd
densf1 <- function(x,pars) 0.5*dunif(x,a,b) + 0.5*dtruncnorm(x,a=a,b=b,mean=pars[1],sd=pars[2])
likef1 <- function(pars) -sum(log(densf1(xx_k2,pars)))
MLEf1 <- constrOptim(theta=c(mean(xx_k2),sd(xx_k2)),f=likef1,ui=c(0,1),ci=c(0),method = "Nelder-Mead")$par
curve(densf1(x,MLEf1),xlim=c(a,b),ylim=c(0,1),col="orange",add=T,lty=2,lwd=4.5)
G_f1 <- function(x) integrate(function(x)densf1(x,MLEf1),a,x)$value
G_f1 <- Vectorize(G_f1)

# Postulated model F2 is a Laplace distribution with unknown location and scale
densf2 <- function(x,pars) 0.3*dtruncnorm(x,a=a,b=b,mean=-5,sd=1) + 0.7*dtruncnorm(x,a=a,b=b,mean=pars[1],sd=pars[2]) 
likef2 <- function(pars) -sum(log(densf2(xx_k2,pars)))
MLEf2 <- constrOptim(theta=c(mean(xx_k2),sd(xx_k2)),f=likef2,ui=c(0,1),ci=c(0),method = "Nelder-Mead")$par
curve(densf2(x,MLEf2),xlim=c(a,b),ylim=c(0,1),col="grey10",add=T,lty=3,lwd=4.5)
G_f2 <- function(x) integrate(function(x)densf2(x,MLEf2),a,x)$value
G_f2 <- Vectorize(G_f2)

# F3 is a mixture of unif[a,b], truncated normal and truncated laplace, with the wrong parameters and unknown proportions 
densf3 <- function(x,pars) (1-pars[1]-pars[2])*dunif(x,a,b) + pars[1]*dtruncnorm(x,a,b,mean=-4,sd=1) + pars[2]*dtrunc(x,spec="laplace",location=4,scale=1,a=a,b=b)
likef3 <- function(pars) -sum(log(densf3(xx_k2,pars)))
MLEf3 <- constrOptim(theta=c(0.3,0.5), f=likef3, ui= rbind(diag(2),-diag(2),c(1,1),c(-1,-1)),
                     ci=c(rep(0,2),rep(-1,2),0,-1), method = "Nelder-Mead", control = list(abstol=1e-15,reltol=1e-15))$par
curve(densf3(x,MLEf3),xlim=c(a,b),ylim=c(0,1),col="orchid3",add=T,cex=1.5,lty=4,lwd=4.5)
G_f3 <- function(x) integrate(function(x)densf3(x,MLEf3),a,x)$value
G_f3 <- Vectorize(G_f3)

legend("topleft", legend=c("Q", expression(F[gamma]), expression(G[beta * ',' * 1]), expression(G[beta * ',' * 2]), expression(G[beta * ',' * 3])), lty=c(5,1:4), col= c("red","blue3","orange","grey","orchid3"), lwd=4.5, cex=1.3, bty = "n")

setEPS()
postscript("hist.eps",width = 9, height = 6)
# par(mar = c(5, 6.5, 2, 2))
hist(xx_k2,breaks=40,prob=T,ylim=c(0,0.16), main="",cex=1.5,cex.axis=1.5,lwd=4.5,cex.lab=1.5,xlim=c(-10,10),xlab=c("x"))
curve(densp(x),xlim=c(a,b),ylim=c(0,1),add=T,col="red",lty=5,lwd=4.5)
curve(densq(x,MLEq),xlim=c(a,b),ylim=c(0,1),add=T,col="blue3",lty=1,lwd=4.5)
curve(densf1(x,MLEf1),xlim=c(a,b),ylim=c(0,1),col="orange",add=T,lty=2,lwd=4.5)
curve(densf2(x,MLEf2),xlim=c(a,b),ylim=c(0,1),col="grey10",add=T,lty=3,lwd=4.5)
curve(densf3(x,MLEf3),xlim=c(a,b),ylim=c(0,1),col="orchid3",add=T,cex=1.5,lty=4,lwd=4.5)
legend("topleft", legend=c("Q", expression(F[gamma]), expression(G[beta * ',' * 1]), expression(G[beta * ',' * 2]), expression(G[beta * ',' * 3])), lty=c(5,1:4), col= c("red","blue3","orange","grey10","orchid3"), lwd=4.5, cex=1.3, bty = "n")
dev.off()

################################## Before K2 ################################## 

cc=seq(0,8,length=1000)
prob_tn_boot_q <- prob_tn_boot_f1 <- prob_tn_boot_f2 <- prob_tn_boot_f3 <-numeric(1000)
prob_tn_order_boot_q <- prob_tn_order_boot_f1 <- prob_tn_order_boot_f2 <- prob_tn_order_boot_f3 <-numeric(1000)
for (k in 1:1000){
  # q
  prob_tn_boot_q[k] <- sum(res_tn_boot_q<cc[k])/length(res_tn_boot_q)
  prob_tn_order_boot_q[k] <- sum(res_tn_order_boot_q<cc[k])/length(res_tn_order_boot_q)
  # f1
  # prob_ks_boot_f1[k] <- sum(res_ks_boot_f1<cc[k])/length(res_ks_boot_f1)
  prob_tn_boot_f1[k] <- sum(res_tn_boot_f1<cc[k])/length(res_tn_boot_f1)
  prob_tn_order_boot_f1[k] <- sum(res_tn_order_boot_f1<cc[k])/length(res_tn_order_boot_q)
  # f2
  # prob_ks_boot_f2[k] <- sum(res_ks_boot_f2<cc[k])/length(res_ks_boot_f2)
  prob_tn_boot_f2[k] <- sum(res_tn_boot_f2<cc[k])/length(res_tn_boot_f2)
  prob_tn_order_boot_f2[k] <- sum(res_tn_order_boot_f2<cc[k])/length(res_tn_order_boot_q)
  # f3
  # prob_ks_boot_f3[k] <- sum(res_ks_boot_f3<cc[k])/length(res_ks_boot_f3)
  prob_tn_boot_f3[k] <- sum(res_tn_boot_f3<cc[k])/length(res_tn_boot_f3)
  prob_tn_order_boot_f3[k] <- sum(res_tn_order_boot_f3<cc[k])/length(res_tn_order_boot_q)
}
postscript("tnbefore.eps",width = 9, height = 6)
par(mar = c(5, 6.5, 2, 2))
plot(cc,prob_tn_boot_q,type = "l", lty=1,col="blue3",xlab="c",cex=2.5,cex.axis=2.5,lwd=4,cex.lab=2.5,
     ylab=expression('P('*hat(T)["n"]<='c)'))
lines(cc,prob_tn_boot_f1,type = "l",col="orange",lty=3,lwd=4)
lines(cc,prob_tn_boot_f2,type = "l",col="grey10",lty=2,lwd=4)
lines(cc,prob_tn_boot_f3,type = "l",col="orchid3",lty=4,lwd=4)
legend("bottomright", legend=c(expression(F[gamma]), expression(G[beta * ',' * 1]), expression(G[beta * ',' * 2]), expression(G[beta * ',' * 3])),lty=1:4,col  = c("blue3","orange","grey10","orchid3"), lwd=4,cex=2,bty = "n")
dev.off()

cc=seq(0,20,length=1000)
for (k in 1:1000){
  # q
  prob_tn_order_boot_q[k] <- sum(res_tn_order_boot_q<cc[k])/length(res_tn_order_boot_q)
  # f1
  prob_tn_order_boot_f1[k] <- sum(res_tn_order_boot_f1<cc[k])/length(res_tn_order_boot_q)
  # f2
  prob_tn_order_boot_f2[k] <- sum(res_tn_order_boot_f2<cc[k])/length(res_tn_order_boot_q)
  # f3
  prob_tn_order_boot_f3[k] <- sum(res_tn_order_boot_f3<cc[k])/length(res_tn_order_boot_q)
}
postscript("subtnbefore.eps",width = 9, height = 6)
par(mar = c(5, 6.5, 2, 2))
plot(cc,prob_tn_order_boot_q,type = "l", lty=1,col="blue3",xlab="c",cex=2.5,cex.axis=2.5,lwd=4,cex.lab=2.5,
     ylab=expression('P('*tilde(T)["n"]<='c)'))
lines(cc,prob_tn_order_boot_f1,type = "l",col="orange",lty=3,lwd=4)
lines(cc,prob_tn_order_boot_f2,type = "l",col="grey10",lty=2,lwd=4)
lines(cc,prob_tn_order_boot_f3,type = "l",col="orchid3",lty=4,lwd=4)
legend("bottomright", legend=c(expression(F[gamma]), expression(G[beta * ',' * 1]), expression(G[beta * ',' * 2]), expression(G[beta * ',' * 3])),lty=1:4,col  = c("blue3","orange","grey10","orchid3"), lwd=4,cex=2,bty = "n")
dev.off()


################################## After K2 ################################## 
cc=seq(0,4,length=1000)
prob_tn_k2_q <- prob_tn_k2_f1 <- prob_tn_k2_f2 <- prob_tn_k2_f3<- numeric(1000)
for (k in 1:1000){
  prob_tn_k2_q[k] <- sum(res_tn_boot_q<cc[k])/length(res_tn_boot_q)
  prob_tn_k2_f1[k] <- sum(res_tn_boot_K2_f1<cc[k])/length(res_tn_boot_K2_f1)
  prob_tn_k2_f2[k] <- sum(res_tn_boot_K2_f2<cc[k])/length(res_tn_boot_K2_f2)
  prob_tn_k2_f3[k] <- sum(res_tn_boot_K2_f3<cc[k])/length(res_tn_boot_K2_f3)
}
postscript("tnafter.eps",width = 9, height = 6)
par(mar = c(5, 6.5, 2, 2))
plot(cc,prob_tn_k2_q,type = "l", lty=1,col="blue3",xlab="c",cex=2.5,cex.axis=2.5,lwd=4,cex.lab=2.5,
     ylab=expression('P('*hat(T)["n"]<='c)'))
lines(cc,prob_tn_k2_f1,type = "l",col="orange",lty=3,lwd=4)
lines(cc,prob_tn_k2_f2,type = "l",col="grey10",lty=2,lwd=4)
lines(cc,prob_tn_k2_f3,type = "l",col="orchid3",lty=4,lwd=4)
legend("bottomright", legend=c(expression(F[gamma]), expression(G[beta * ',' * 1]), expression(G[beta * ',' * 2]), expression(G[beta * ',' * 3])),lty=1:4,col  = c("blue3","orange","grey10","orchid3"), lwd=4,cex=2,bty = "n")
dev.off()

cc=seq(0,12,length=1000)
prob_tn_order_k2_q <- prob_tn_order_k2_f1 <- prob_tn_order_k2_f2 <- prob_tn_order_k2_f3<- numeric(1000)
for (k in 1:1000){
  prob_tn_order_k2_q[k] <- sum(res_tn_order_boot_q<cc[k])/length(res_tn_order_boot_q)
  prob_tn_order_k2_f1[k] <- sum(res_tn_order_boot_K2_f1<cc[k])/length(res_tn_order_boot_K2_f1)
  prob_tn_order_k2_f2[k] <- sum(res_tn_order_boot_K2_f2<cc[k])/length(res_tn_order_boot_K2_f2)
  prob_tn_order_k2_f3[k] <- sum(res_tn_order_boot_K2_f3<cc[k])/length(res_tn_order_boot_K2_f3)
}
postscript("subtnafter.eps",width = 9, height = 6)
par(mar = c(5, 6.5, 2, 2))
plot(cc,prob_tn_order_k2_q,type = "l", lty=1,col="blue3",xlab="c",cex=2.5,cex.axis=2.5,lwd=4,
     cex.lab=2.5,ylab=expression('P('*tilde(T)["n"]<='c)'))
lines(cc,prob_tn_order_k2_f1,type = "l",col="orange",lty=3,lwd=4)
lines(cc,prob_tn_order_k2_f2,type = "l",col="grey10",lty=2,lwd=4)
lines(cc,prob_tn_order_k2_f3,type = "l",col="orchid3",lty=4,lwd=4)
legend("bottomright", legend=c(expression(F[gamma]), expression(G[beta * ',' * 1]), expression(G[beta * ',' * 2]), expression(G[beta * ',' * 3])),lty=1:4,col  = c("blue3","orange","grey10","orchid3"), lwd=4,cex=2,bty = "n")
dev.off()

