source("hjs.R")
library(nloptr)
library(LaplacesDemon)
library(expm)
library(truncnorm)
library(orthopolynom)
library(pracma)
library(tolerance)

# Orginal data generation
set.seed(1)
n=100; p1=0.3; p2=0.5
a=-10;b=10 #truncated
indx = rmultinom(n=1,size=n,prob=c(1-p1-p2, p1,p2))
xx_k2 = c(runif(indx[1],a,b), rtruncnorm(indx[2],a,b,mean=-5,sd=3),
       rtrunc(indx[3],spec="laplace",location=5,scale=3,a=a,b=b))
hist(xx_k2,breaks=40,prob=T,ylim=c(0,0.3))
M1 = 6

# Reference distribution is a truncated normal distribution with unknown mean and sd
densq <- function(x,pars) dtruncnorm(x,a=a,b=b,mean=pars[1],sd=pars[2])
likeq <- function(pars) -sum(log(densq(xx_k2,pars)))
MLEq <- constrOptim(theta=c(mean(xx_k2),sd(xx_k2)),f=likeq,ui=c(0,1),ci=c(0),method = "Nelder-Mead")$par
curve(densq(x,MLEq),xlim=c(a,b),ylim=c(0,1),add=T,col="blue3")
G_obsq = function(x) ptruncnorm(x,a=a,b=b,mean=MLEq[1],sd=MLEq[2])

# Postulated model F1 is a mixture of half unif[a,b], half truncated normal with unknown mean and sd
densf1 <- function(x,pars) 0.5*dunif(x,a,b) + 0.5*dtruncnorm(x,a=a,b=b,mean=pars[1],sd=pars[2])
likef1 <- function(pars) -sum(log(densf1(xx_k2,pars)))
MLEf1 <- constrOptim(theta=c(mean(xx_k2),sd(xx_k2)),f=likef1,ui=c(0,1),ci=c(0),method = "Nelder-Mead")$par
curve(densf1(x,MLEf1),xlim=c(a,b),ylim=c(0,1),col="orange",add=T)
G_f1 <- function(x) integrate(function(x)densf1(x,MLEf1),a,x)$value
G_f1 <- Vectorize(G_f1)

# Postulated model F2 is a Laplace distribution with unknown location and scale
densf2 <- function(x,pars) 0.3*dtruncnorm(x,a=a,b=b,mean=-5,sd=1) + 0.7*dtruncnorm(x,a=a,b=b,mean=pars[1],sd=pars[2]) 
likef2 <- function(pars) -sum(log(densf2(xx_k2,pars)))
MLEf2 <- constrOptim(theta=c(mean(xx_k2),sd(xx_k2)),f=likef2,ui=c(0,1),ci=c(0),method = "Nelder-Mead")$par
curve(densf2(x,MLEf2),xlim=c(a,b),ylim=c(0,1),col="grey10",add=T)
G_f2 <- function(x) integrate(function(x)densf2(x,MLEf2),a,x)$value
G_f2 <- Vectorize(G_f2)

# F3 is a mixture of unif[a,b], truncated normal and truncated laplace, with the wrong parameters and unknown proportions 
densf3 <- function(x,pars) (1-pars[1]-pars[2])*dunif(x,a,b) + pars[1]*dtruncnorm(x,a,b,mean=-4,sd=1) + pars[2]*dtrunc(x,spec="laplace",location=4,scale=1,a=a,b=b)
likef3 <- function(pars) -sum(log(densf3(xx_k2,pars)))
MLEf3 <- constrOptim(theta=c(0.3,0.5), f=likef3, ui= rbind(diag(2),-diag(2),c(1,1),c(-1,-1)),
                     ci=c(rep(0,2),rep(-1,2),0,-1), method = "Nelder-Mead", control = list(abstol=1e-15,reltol=1e-15))$par
curve(densf3(x,MLEf3),xlim=c(a,b),ylim=c(0,1),col="darkgreen",add=T)
G_f3 <- function(x) integrate(function(x)densf3(x,MLEf3),a,x)$value
G_f3 <- Vectorize(G_f3)


################## Orthornormalzied Score functions for Q, F1, F2 and F3 ##################
norm_score_q <- function(pars){
  mu<-pars[1]; sigma <- pars[2]
  grad1 <- function(x) (x-mu)/sigma^2 + (dnorm((b-mu)/sigma)/sigma-dnorm((a-mu)/sigma)/sigma)/(pnorm(b,mu,sigma)-pnorm(a,mu,sigma))
  grad2 <- function(x) (x-mu)^2/sigma^3 -1/sigma + ((b-mu)/sigma^2 *dnorm((b-mu)/sigma) - (a-mu)/sigma^2 *dnorm((a-mu)/sigma))/(pnorm(b,mu,sigma)-pnorm(a,mu,sigma))
  fish.exact <- matrix(NA, 2, 2)
  fish.exact[1, 1] <- integrate(function(x) {grad1(x)^2 *densq(x,pars)}, a, b)$value
  fish.exact[2, 2] <- integrate(function(x) {grad2(x)^2 *densq(x,pars)}, a, b)$value
  fish.exact[1, 2] <- integrate(function(x) {grad1(x)*grad2(x) *densq(x,pars)}, a, b)$value
  fish.exact[2, 1] <- fish.exact[1, 2]
  
  ## Sqaure root Inverse Fisher times score functions 
  fisher_sqroot_inverse <- pinv(expm::sqrtm(fish.exact))
  norm_score_func = function(x) fisher_sqroot_inverse %*% c(grad1(x),grad2(x)) 
  norm_score_func1 = function(x) as.numeric(fisher_sqroot_inverse[1,] %*% c(grad1(x),grad2(x)))
  norm_score_func2 = function(x) as.numeric(fisher_sqroot_inverse[2,] %*% c(grad1(x),grad2(x)))
  norm_score_func = Vectorize(norm_score_func)
  norm_score_func1 = Vectorize(norm_score_func1)
  norm_score_func2 = Vectorize(norm_score_func2)
  return (list(norm_score_func=norm_score_func,
               norm_score_func1=norm_score_func1,
               norm_score_func2=norm_score_func2))
}
norm_score_q1 = norm_score_q(MLEq)$norm_score_func1
norm_score_q2 = norm_score_q(MLEq)$norm_score_func2

norm_score_f1 <- function(pars){
  mu<-pars[1]; sigma <- pars[2]
  grad1 <- function(x) ((x-mu)/sigma^2 + (dnorm((b-mu)/sigma)/sigma-dnorm((a-mu)/sigma)/sigma)/(pnorm(b,mu,sigma)-pnorm(a,mu,sigma)))*0.5*dtruncnorm(x,a=a,b=b,mean=mu,sd=sigma)/densf1(x,pars)
  grad2 <- function(x) ((x-mu)^2/sigma^3 -1/sigma + ((b-mu)/sigma^2 *dnorm((b-mu)/sigma) - (a-mu)/sigma^2 *dnorm((a-mu)/sigma))/(pnorm(b,mu,sigma)-pnorm(a,mu,sigma)))*0.5*dtruncnorm(x,a=a,b=b,mean=mu,sd=sigma)/densf1(x,pars)
  
  fish.exact <- matrix(NA, 2, 2)
  fish.exact[1, 1] <- integrate(function(x) {grad1(x)^2 *densf1(x,pars)}, a, b)$value
  fish.exact[2, 2] <- integrate(function(x) {grad2(x)^2 *densf1(x,pars)}, a, b)$value
  fish.exact[1, 2] <- integrate(function(x) {grad1(x)*grad2(x) *densf1(x,pars)}, a, b)$value
  fish.exact[2, 1] <- fish.exact[1, 2]
  
  ## Sqaure root Inverse Fisher times score functions 
  fisher_sqroot_inverse <- pinv(expm::sqrtm(fish.exact))
  norm_score_func = function(x) fisher_sqroot_inverse %*% c(grad1(x),grad2(x)) 
  norm_score_func1 = function(x) as.numeric(fisher_sqroot_inverse[1,] %*% c(grad1(x),grad2(x)))
  norm_score_func2 = function(x) as.numeric(fisher_sqroot_inverse[2,] %*% c(grad1(x),grad2(x)))
  norm_score_func = Vectorize(norm_score_func)
  norm_score_func1 = Vectorize(norm_score_func1)
  norm_score_func2 = Vectorize(norm_score_func2)
  return (list(norm_score_func=norm_score_func,
               norm_score_func1=norm_score_func1,
               norm_score_func2=norm_score_func2))
}
norm_score_f11 = norm_score_f1(MLEf1)$norm_score_func1
norm_score_f12 = norm_score_f1(MLEf1)$norm_score_func2

norm_score_f2 <- function(pars){
  mu<-pars[1]; sigma <- pars[2]
  grad1 <- function(x) ((x-mu)/sigma^2 + (dnorm((b-mu)/sigma)/sigma-dnorm((a-mu)/sigma)/sigma)/(pnorm(b,mu,sigma)-pnorm(a,mu,sigma)))*0.7*dtruncnorm(x,a=a,b=b,mean=mu,sd=sigma)/densf2(x,pars)
  grad2 <- function(x) ((x-mu)^2/sigma^3 -1/sigma + ((b-mu)/sigma^2 *dnorm((b-mu)/sigma) - (a-mu)/sigma^2 *dnorm((a-mu)/sigma))/(pnorm(b,mu,sigma)-pnorm(a,mu,sigma)))*0.7*dtruncnorm(x,a=a,b=b,mean=mu,sd=sigma)/densf2(x,pars)
  
  fish.exact <- matrix(NA, 2, 2)
  fish.exact[1, 1] <- integrate(function(x) {grad1(x)^2 *densf2(x,pars)}, a, b, rel.tol = 1e-12)$value
  fish.exact[2, 2] <- integrate(function(x) {grad2(x)^2 *densf2(x,pars)}, a, b, rel.tol = 1e-12)$value
  fish.exact[1, 2] <- integrate(function(x) {grad1(x)*grad2(x) *densf2(x,pars)}, a, b, rel.tol = 1e-12)$value
  fish.exact[2, 1] <- fish.exact[1, 2]
  
  ## Sqaure root Inverse Fisher times score functions 
  fisher_sqroot_inverse <- pinv(expm::sqrtm(fish.exact))
  norm_score_func = function(x) fisher_sqroot_inverse %*% c(grad1(x),grad2(x)) 
  norm_score_func1 = function(x) as.numeric(fisher_sqroot_inverse[1,] %*% c(grad1(x),grad2(x)))
  norm_score_func2 = function(x) as.numeric(fisher_sqroot_inverse[2,] %*% c(grad1(x),grad2(x)))
  norm_score_func = Vectorize(norm_score_func)
  norm_score_func1 = Vectorize(norm_score_func1)
  norm_score_func2 = Vectorize(norm_score_func2)
  return (list(norm_score_func=norm_score_func,
               norm_score_func1=norm_score_func1,
               norm_score_func2=norm_score_func2))
}
norm_score_f21 = norm_score_f2(MLEf2)$norm_score_func1
norm_score_f22 = norm_score_f2(MLEf2)$norm_score_func2

norm_score_f3 <- function(pars){
  grad1 <- function(x) (-dunif(x,a,b) + dtruncnorm(x,a,b,mean=-4,sd=1))/densf3(x,pars)
  grad2 <- function(x) (-dunif(x,a,b) + dtrunc(x,spec="laplace",location=4,scale=1,a=a,b=b))/densf3(x,pars)
  fish.exact <- matrix(NA, 2, 2)
  fish.exact[1, 1] <- integrate(function(x) {grad1(x)^2 *densf3(x,pars)}, a, b)$value
  fish.exact[2, 2] <- integrate(function(x) {grad2(x)^2 *densf3(x,pars)}, a, b)$value
  fish.exact[1, 2] <- integrate(function(x) {grad1(x)*grad2(x) *densf3(x,pars)}, a, b)$value
  fish.exact[2, 1] <- fish.exact[1, 2]
  
  ## Sqaure root Inverse Fisher times score functions 
  fisher_sqroot_inverse <- pinv(expm::sqrtm(fish.exact))
  norm_score_func = function(x) fisher_sqroot_inverse %*% c(grad1(x),grad2(x)) 
  norm_score_func1 = function(x) as.numeric(fisher_sqroot_inverse[1,] %*% c(grad1(x),grad2(x)))
  norm_score_func2 = function(x) as.numeric(fisher_sqroot_inverse[2,] %*% c(grad1(x),grad2(x)))
  norm_score_func = Vectorize(norm_score_func)
  norm_score_func1 = Vectorize(norm_score_func1)
  norm_score_func2 = Vectorize(norm_score_func2)
  return (list(norm_score_func=norm_score_func,
               norm_score_func1=norm_score_func1,
               norm_score_func2=norm_score_func2))
}
norm_score_f31 = norm_score_f3(MLEf3)$norm_score_func1
norm_score_f32 = norm_score_f3(MLEf3)$norm_score_func2

vj=list()
vj[[1]] = function(x) return (1)
for (i in 1:(M1)){
  call1 = paste("inner",i,"_1=integrate(function(x) {hj[[",i+1
                ,"]](G_obsq(x))*norm_score_q1(x)*densq(x,MLEq) }, a, b)$value",sep="")
  call2 = paste("inner",i,"_2=integrate(function(x) {hj[[",i+1
                ,"]](G_obsq(x))*norm_score_q2(x)*densq(x,MLEq) }, a, b)$value",sep="")
  call3 = paste("vj[[",i+1,"]]= function(x) {hj[[",i+1,"]](G_obsq(x)) - norm_score_q1(x) * inner",i,
                "_1 - norm_score_q2(x) * inner",i,"_2}",sep="")
  eval(parse(text=call1))
  eval(parse(text=call2))
  eval(parse(text=call3))
  print(i)
}

vj_f1=list(); vj_f1[[1]] = function(x) return (1)
for (i in 1:(M1)){
  call1 = paste("inner_f1",i,"_1=integrate(function(x) {hj[[",i+1
                ,"]](G_f1(x))*norm_score_f11(x)*densf1(x,MLEf1) }, a, b)$value",sep="")
  call2 = paste("inner_f1",i,"_2=integrate(function(x) {hj[[",i+1
                ,"]](G_f1(x))*norm_score_f12(x)*densf1(x,MLEf1) }, a, b)$value",sep="")
  call3 = paste("vj_f1[[",i+1,"]]= function(x) {hj[[",i+1,"]](G_f1(x)) - norm_score_f11(x) * inner_f1",i,
                "_1 - norm_score_f12(x) * inner_f1",i,"_2}",sep="")
  eval(parse(text=call1))
  eval(parse(text=call2))
  eval(parse(text=call3))
  print(i)
}
vj_f2=list(); vj_f2[[1]] = function(x) return (1)
for (i in 1:(M1)){
  call1 = paste("inner_f2",i,"_1=integrate(function(x) {hj[[",i+1
                ,"]](G_f2(x))*norm_score_f21(x)*densf2(x,MLEf2) }, a, b)$value",sep="")
  call2 = paste("inner_f2",i,"_2=integrate(function(x) {hj[[",i+1
                ,"]](G_f2(x))*norm_score_f22(x)*densf2(x,MLEf2) }, a, b)$value",sep="")
  call3 = paste("vj_f2[[",i+1,"]]= function(x) {hj[[",i+1,"]](G_f2(x)) - norm_score_f21(x) * inner_f2",i,
                "_1 - norm_score_f22(x) * inner_f2",i,"_2}",sep="")
  eval(parse(text=call1))
  eval(parse(text=call2))
  eval(parse(text=call3))
  print(i)
}
vj_f3=list(); vj_f3[[1]] = function(x) return (1)
for (i in 1:(M1)){
  call1 = paste("inner_f3",i,"_1=integrate(function(x) {hj[[",i+1
                ,"]](G_f3(x))*norm_score_f31(x)*densf3(x,MLEf3) }, a, b)$value",sep="")
  call2 = paste("inner_f3",i,"_2=integrate(function(x) {hj[[",i+1
                ,"]](G_f3(x))*norm_score_f32(x)*densf3(x,MLEf3) }, a, b)$value",sep="")
  call3 = paste("vj_f3[[",i+1,"]]= function(x) {hj[[",i+1,"]](G_f3(x)) - norm_score_f31(x) * inner_f3",i,
                "_1 - norm_score_f32(x) * inner_f3",i,"_2}",sep="")
  eval(parse(text=call1))
  eval(parse(text=call2))
  eval(parse(text=call3))
  print(i)
}


################### Steps before K2 ################### 
B=100000
# Starting from here, we do the bootstrap 
res_tn_boot_q <- res_tn_order_boot_q <- numeric(B); 
vj_boot_q <- numeric(M1)
for (i in 1:B){
  xx_q= rtruncnorm(n,a=a,b=b,mean=MLEq[1],sd=MLEq[2])
  # likeq <- function(pars) -sum(log(densq(xx_q,pars)))
  # MLEq_boot <- constrOptim(theta=c(mean(xx_q),sd(xx_q)),f=likeq,ui=c(0,1),ci=c(0),method = "Nelder-Mead")$par
  # G_boot_q <- function(x) ptruncnorm(x,mean=MLEq_boot[1],sd=MLEq_boot[2],a=a,b=b)
  for (j in 1:M1){
    vj_boot_q[j] <- sqrt(n)*mean(vj[[j+1]](xx_q))
  }
  res_tn_boot_q[i] <- max(cumsum((vj_boot_q)^2)/1:M1)
  res_tn_order_boot_q[i] <- max(cumsum(sort((vj_boot_q)^2,decreasing=TRUE)/1:M1))
  print(i)
}
res_tn_boot_f1 <-  res_tn_order_boot_f1 <- numeric(B)
vj_boot_f1 <- numeric(M1)
for (i in 1:B){
  indx_f1 = rmultinom(n=1,size=n,prob=c(1-0.5, 0.5))
  xx_f1 = c(runif(indx_f1[1],a,b), rtruncnorm(indx_f1[2],a=a, b=b,mean=MLEf1[1],sd=MLEf1[2])
            #rtrunc(indx_f1[2],spec="laplace",location=MLEf1[1],scale=MLEf1[2],a=a,b=b)
            )
  # likef1_boot <- function(pars) -sum(log(densf1(xx_f1,pars)))
  # MLEf1_boot <- constrOptim(theta=c(median(xx_f1),sd(xx_f1)),f=likef1_boot,ui=c(0,1),ci=c(0),method = "Nelder-Mead")$par
  # G_boot_f1 <- function(x) integrate(function(x)densf1(x,MLEf1_boot),a,x)$value
  # G_boot_f1 <- Vectorize(G_boot_f1)
  for (j in 1:M1){
    vj_boot_f1[j] <- sqrt(n)*mean(vj_f1[[j+1]](xx_f1))
  }
  res_tn_boot_f1[i] <- max(cumsum((vj_boot_f1)^2)/1:M1)
  res_tn_order_boot_f1[i] <- max(cumsum(sort((vj_boot_f1)^2, decreasing=TRUE)/1:M1))
  print(i)
}
res_tn_boot_f2 <- res_tn_order_boot_f2 <- numeric(B)
vj_boot_f2 <- numeric(M1)
for (i in 1:B){
  indx_f2 = rmultinom(n=1,size=n,prob=c(0.3, 0.7))
  xx_f2 = c(rtruncnorm(indx_f2[1],a=a,b=b,mean=-5,sd=1),
            rtruncnorm(indx_f2[1],a=a,b=b,mean=MLEf2[1],sd=MLEf2[2]))
  # likef2_boot <- function(pars) -sum(log(densf2(xx_f2,pars)))
  # MLEf2_boot <- constrOptim(theta=c(mean(xx_k2),sd(xx_k2)),f=likef2_boot,ui=c(0,1),ci=c(0),method = "Nelder-Mead")$par
  # G_boot_f2 <- function(x) integrate(function(x)densf2(x,MLEf2_boot),a,x)$value
  # G_boot_f2 <- Vectorize(G_boot_f2)
  for (j in 1:M1){
    vj_boot_f2[j] <- sqrt(n)*mean(vj_f2[[j+1]](xx_f2))
  }
  res_tn_boot_f2[i] <- max(cumsum((vj_boot_f2)^2)/1:M1)
  res_tn_order_boot_f2[i] <- max(cumsum(sort((vj_boot_f2)^2, decreasing=TRUE)/1:M1))
  print(i)
}
res_tn_boot_f3 <- res_tn_order_boot_f3 <- numeric(B)
vj_boot_f3 <- numeric(M1)
for (i in 1:B){
  indx_f3 = rmultinom(n=1,size=n,prob=c(1-MLEf3[1]-MLEf3[2], MLEf3[1], MLEf3[2]))
  xx_f3 = c(runif(indx_f3[1],a,b), rtruncnorm(indx_f3[2],a,b,mean=-4,sd=1),
            rtrunc(indx_f3[3],spec="laplace",location=4,scale=1,a=a,b=b))
  # likef3_boot <- function(pars) -sum(log(densf3(xx_f3,pars)))
  # MLEf3_boot <- constrOptim(theta=MLEf3, f=likef3_boot, ui= rbind(diag(2),-diag(2),c(1,1),c(-1,-1)),
  #                           ci=c(rep(0,2),rep(-1,2),0,-1), method = "Nelder-Mead", control = list(abstol=1e-12,reltol=1e-12))$par
  # G_boot_f3 <- function(x) integrate(function(x) densf3(x,MLEf3_boot),a,x)$value
  # G_boot_f3 <- Vectorize(G_boot_f3)
  for (j in 1:M1){
    vj_boot_f3[j] <- sqrt(n)*mean(vj_f3[[j+1]](xx_f3))
  }
  res_tn_boot_f3[i] <- max(cumsum((vj_boot_f3)^2)/1:M1)
  res_tn_order_boot_f3[i] <- max(cumsum(sort((vj_boot_f3)^2, decreasing=TRUE)/1:M1))
  
  print(i)
}


cc=seq(0,5,length=1000)
prob_tn_boot_q <- prob_tn_boot_f1 <- prob_tn_boot_f2 <- prob_tn_boot_f3 <-numeric(1000)
prob_tn_order_boot_q <- prob_tn_order_boot_f1 <- prob_tn_order_boot_f2 <- prob_tn_order_boot_f3 <-numeric(1000)
for (k in 1:1000){
  # q
  prob_tn_boot_q[k] <- sum(res_tn_boot_q<cc[k])/length(res_tn_boot_q)
  prob_tn_order_boot_q[k] <- sum(res_tn_order_boot_q<cc[k])/length(res_tn_order_boot_q)
  # f1
  prob_tn_boot_f1[k] <- sum(res_tn_boot_f1<cc[k])/length(res_tn_boot_f1)
  prob_tn_order_boot_f1[k] <- sum(res_tn_order_boot_f1<cc[k])/length(res_tn_order_boot_q)
  # f2
  prob_tn_boot_f2[k] <- sum(res_tn_boot_f2<cc[k])/length(res_tn_boot_f2)
  prob_tn_order_boot_f2[k] <- sum(res_tn_order_boot_f2<cc[k])/length(res_tn_order_boot_q)
  # f3
  prob_tn_boot_f3[k] <- sum(res_tn_boot_f3<cc[k])/length(res_tn_boot_f3)
  prob_tn_order_boot_f3[k] <- sum(res_tn_order_boot_f3<cc[k])/length(res_tn_order_boot_q)
}

plot(cc,prob_tn_boot_q,type = "l", lty=1,xlab="c",cex=2,ylim=c(0,1),cex.axis=1.3,lwd=1.5,cex.lab=1.5,ylab="P(K<=c)")
lines(cc,prob_tn_boot_f1,type = "l",col="red",lty=2,lwd=1.5)
lines(cc,prob_tn_boot_f2,type = "l",col="blue",lty=3,lwd=1.5)
lines(cc,prob_tn_boot_f3,type = "l",col="darkgreen",lty=4,lwd=1.5)
legend("bottomright", legend=c("Normal",
                               "0.5Uniform+0.5Normal",
                               "Laplace",
                               "True"),lty=1:4,col = c("black","red","blue","darkgreen"),lwd=1.3,cex=1.3,bty = "n")

cc=seq(0,10,length=1000)
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
plot(cc,prob_tn_order_boot_q,type = "l", lty=1,xlab="c",cex=2,ylim=c(0,1),cex.axis=1.3,lwd=1.5,cex.lab=1.5,ylab="P(K<=c)")
lines(cc,prob_tn_order_boot_f1,type = "l",col="red",lty=2,lwd=1.5)
lines(cc,prob_tn_order_boot_f2,type = "l",col="blue",lty=3,lwd=1.5)
lines(cc,prob_tn_order_boot_f3,type = "l",col="darkgreen",lty=4,lwd=1.5)
legend("bottomright", legend=c("Normal",
                               "0.5Uniform+0.5Normal",
                               "Laplace",
                               "True"),lty=1:4,col = c("black","red","blue","darkgreen"),lwd=1.5,cex=1.3,bty = "n")

##################################K2-STEP1##################################
isometry1=function(x) sqrt(densq(x,MLEq)/densf1(x,MLEf1))
isometry2=function(x) sqrt(densq(x,MLEq)/densf2(x,MLEf2))
isometry3=function(x) sqrt(densq(x,MLEq)/densf3(x,MLEf3))


l1h=list();l2h=list();l3h=list()
for (i in 1:(M1)){
  for (j in 1:3){
    call = paste("l",j,"h[[",i,"]]<-function(x) isometry",j,"(x) * hj[[",i+1,"]](G_obsq(x))",sep="")
    eval(parse(text=call))
  }
}

l1b=list();l2b=list();l3b=list()
for (i in 1:2){
  for (j in 1:3){
    call = paste("l",j,"b[[",i,"]]<-function(x) isometry",j,"(x) *norm_score_q",i,"(x)",sep="")
    eval(parse(text=call))
  }
}

innerl1 = integrate(function(x) {sqrt(densq(x,MLEq)*densf1(x,MLEf1))}, lower=a, upper=b)$value
innerl2 = integrate(function(x) {sqrt(densq(x,MLEq)*densf2(x,MLEf2))}, lower=a, upper=b)$value
innerl3 = integrate(function(x) {sqrt(densq(x,MLEq)*densf3(x,MLEf3))}, lower=a, upper=b)$value

inner1lh = inner2lh = inner3lh = c()
for (i in 1:(M1)){
  inner1lh[i]=integrate(function(x) {l1h[[i]](x)*densf1(x,MLEf1)},lower=a,upper=b)$value
  inner2lh[i]=integrate(function(x) {l2h[[i]](x)*densf2(x,MLEf2)},lower=a,upper=b)$value
  inner3lh[i]=integrate(function(x) {l3h[[i]](x)*densf3(x,MLEf3)},lower=a,upper=b)$value
  }

inner1lb = inner2lb = inner3lb = c()
for (i in 1:2){
  inner1lb[i]=integrate(function(x) {l1b[[i]](x)*densf1(x,MLEf1)},lower=a,upper=b)$value
  inner2lb[i]=integrate(function(x) {l2b[[i]](x)*densf2(x,MLEf2)},lower=a,upper=b)$value
  inner3lb[i]=integrate(function(x) {l3b[[i]](x)*densf3(x,MLEf3)},lower=a,upper=b)$value
  }

Kl1h = Kl2h = Kl3h = list()
for (i in 1:(M1)){
  for (j in 1:3){
    call = paste("Kl",j,"h[[",i,"]]<-function(x) l",j,"h[[",i,"]](x)- (1-isometry",j,"(x))/(1-innerl",j,")*inner",j,"lh[",i,"]",sep="")
    eval(parse(text=call))
  }
}
Kl1b = Kl2b = Kl3b = list()
for (i in 1:2){
  for (j in 1:3){
    call = paste("Kl",j,"b[[",i,"]]<-function(x) l",j,"b[[",i,"]](x)- (1-isometry",j,"(x))/(1-innerl",j,")*inner",j,"lb[",i,"]",sep="")
    eval(parse(text=call))
  }
}

inner_step3_1_1=integrate(function(x) {(norm_score_f11(x)-Kl1b[[1]](x))*Kl1b[[2]](x)*densf1(x,MLEf1)},lower=a,upper=b)$value
inner_step3_2_1=integrate(function(x) {norm_score_f11(x)*Kl1b[[1]](x)*densf1(x,MLEf1)},lower=a,upper=b)$value

inner_step3_1_2=integrate(function(x) {(norm_score_f21(x)-Kl2b[[1]](x))*Kl2b[[2]](x)*densf2(x,MLEf2)},lower=a,upper=b)$value
inner_step3_2_2=integrate(function(x) {norm_score_f21(x)*Kl2b[[1]](x)*densf2(x,MLEf2)},lower=a,upper=b)$value

inner_step3_1_3=integrate(function(x) {(norm_score_f31(x)-Kl3b[[1]](x))*Kl3b[[2]](x)*densf3(x,MLEf3)},lower=a,upper=b)$value
inner_step3_2_3=integrate(function(x) {norm_score_f31(x)*Kl3b[[1]](x)*densf3(x,MLEf3)},lower=a,upper=b)$value


tildec2_f1 = function(x){return (Kl1b[[2]](x)-inner_step3_1_1/(1-inner_step3_2_1)*(norm_score_f11(x)-Kl1b[[1]](x)))} # Same with below
tildec2_f2 = function(x){return (Kl2b[[2]](x)-inner_step3_1_2/(1-inner_step3_2_2)*(norm_score_f21(x)-Kl2b[[1]](x)))} # Same with below
tildec2_f3 = function(x){return (Kl3b[[2]](x)-inner_step3_1_3/(1-inner_step3_2_3)*(norm_score_f31(x)-Kl3b[[1]](x)))} # Same with below


inner1_1h = inner1_2h = inner1_3h = inner_U1_h = inner_U2_h = inner_U3_h = c()
U1_h = U2_h = U3_h = list()
for (i in 1:(M1)){
  inner1_1h[i]=integrate(function(x) {(norm_score_f11(x)-Kl1b[[1]](x))*Kl1h[[i]](x)*densf1(x,MLEf1)},lower=a,upper=b)$value
  inner1_2h[i]=integrate(function(x) {(norm_score_f21(x)-Kl2b[[1]](x))*Kl2h[[i]](x)*densf2(x,MLEf2)},lower=a,upper=b)$value
  inner1_3h[i]=integrate(function(x) {(norm_score_f31(x)-Kl3b[[1]](x))*Kl3h[[i]](x)*densf3(x,MLEf3)},lower=a,upper=b)$value
}
for (i in 1:(M1)){
  call003 = paste("U1_h[[",i,"]]<-function(x) {Kl1h[[",i,"]](x)-inner1_1h[",i,"]/(1-inner_step3_2_1)*(norm_score_f11(x)-Kl1b[[1]](x))}",sep="")
  call004 = paste("U2_h[[",i,"]]<-function(x) {Kl2h[[",i,"]](x)-inner1_2h[",i,"]/(1-inner_step3_2_2)*(norm_score_f21(x)-Kl2b[[1]](x))}",sep="")
  call005 = paste("U3_h[[",i,"]]<-function(x) {Kl3h[[",i,"]](x)-inner1_3h[",i,"]/(1-inner_step3_2_3)*(norm_score_f31(x)-Kl3b[[1]](x))}",sep="")
  eval(parse(text=call003))
  eval(parse(text=call004))
  eval(parse(text=call005))
  inner_U1_h[i] = integrate(function(x) {(norm_score_f12(x)-tildec2_f1(x))*U1_h[[i]](x)*densf1(x,MLEf1)},lower=a,upper=b)$value
  inner_U2_h[i] = integrate(function(x) {(norm_score_f22(x)-tildec2_f2(x))*U2_h[[i]](x)*densf2(x,MLEf2)},lower=a,upper=b)$value
  inner_U3_h[i] = integrate(function(x) {(norm_score_f32(x)-tildec2_f3(x))*U3_h[[i]](x)*densf3(x,MLEf3)},lower=a,upper=b)$value
}
inner_U1=integrate(function(x) {norm_score_f12(x)*tildec2_f1(x)*densf1(x,MLEf1)},lower=a,upper=b)$value
inner_U2=integrate(function(x) {norm_score_f22(x)*tildec2_f2(x)*densf2(x,MLEf2)},lower=a,upper=b)$value
inner_U3=integrate(function(x) {norm_score_f32(x)*tildec2_f3(x)*densf3(x,MLEf3)},lower=a,upper=b)$value

UKl1h=list();UKl2h=list();UKl3h=list()
for (i in 1:(M1)){
  call01 = paste("UKl1h[[",i,"]]<-function(x) U1_h[[",i,"]](x)- inner_U1_h[",i,"]/(1-inner_U1)*(norm_score_f12(x) - tildec2_f1(x))",sep="")
  eval(parse(text=call01))
  call02 = paste("UKl2h[[",i,"]]<-function(x) U2_h[[",i,"]](x)- inner_U2_h[",i,"]/(1-inner_U2)*(norm_score_f22(x) - tildec2_f2(x))",sep="")
  eval(parse(text=call02))
  call03 = paste("UKl3h[[",i,"]]<-function(x) U3_h[[",i,"]](x)- inner_U3_h[",i,"]/(1-inner_U3)*(norm_score_f32(x) - tildec2_f3(x))",sep="")
  eval(parse(text=call03))
}


##################################Calculating tilde h or tilde psi##################################


Kvj1=list()
Kvj1[[1]] = function(x) return (1)
for (i in 1:(M1)){
  call1 = paste("inner",i,"_11=integrate(function(x) {hj[[",i+1
                ,"]](G_obsq(x))*norm_score_q1(x)*densq(x,MLEq) }, a, b)$value",sep="")
  call2 = paste("inner",i,"_12=integrate(function(x) {hj[[",i+1
               ,"]](G_obsq(x))*norm_score_q2(x)*densq(x,MLEq) }, a, b)$value",sep="")
  call3 = paste("Kvj1[[",i+1,"]]= function(x) {UKl1h[[",i,"]](x) - norm_score_f11(x) * inner",i,
                "_11 - norm_score_f12(x) * inner",i,"_12}",sep="")
  eval(parse(text=call1))
  eval(parse(text=call2))
  eval(parse(text=call3))
  print(i)
}

Kvj2=list()
Kvj2[[1]] = function(x) return (1)
for (i in 1:M1){
  call1 = paste("inner",i,"_21=integrate(function(x) {hj[[",i+1
                ,"]](G_obsq(x))*norm_score_q1(x)*densq(x,MLEq) }, a, b)$value",sep="")
  call2 =paste("inner",i,"_22=integrate(function(x) {hj[[",i+1
               ,"]](G_obsq(x))*norm_score_q2(x)*densq(x,MLEq) }, a, b)$value",sep="")
  call3 = paste("Kvj2[[",i+1,"]]= function(x) {UKl2h[[",i,"]](x) - norm_score_f21(x) * inner",i,
                "_21 - norm_score_f22(x) * inner",i,"_22}",sep="")
  eval(parse(text=call1))
  eval(parse(text=call2))
  eval(parse(text=call3))
  print(i)
}

Kvj3=list()
Kvj3[[1]] = function(x) return (1)
for (i in 1:M1){
  call1 = paste("inner",i,"_31=integrate(function(x) {hj[[",i+1
                ,"]](G_obsq(x))*norm_score_q1(x)*densq(x,MLEq) }, a, b)$value",sep="")
  call2 =paste("inner",i,"_32=integrate(function(x) {hj[[",i+1
               ,"]](G_obsq(x))*norm_score_q2(x)*densq(x,MLEq) }, a, b)$value",sep="")
  call3 = paste("Kvj3[[",i+1,"]]= function(x) {UKl3h[[",i,"]](x) - norm_score_f31(x) * inner",i,
                "_31 - norm_score_f32(x) * inner",i,"_32}",sep="")
  eval(parse(text=call1))
  eval(parse(text=call2))
  eval(parse(text=call3))
  print(i)
}

varcov_q=varcov_f1=varcov_f2=varcov_f3=matrix(nrow=M1+1,ncol=M1+1)
for (m in 1:(M1+1)){
  for (s in 1:(M1+1))
    varcov_q[m,s]= integrate(function(x) {vj[[m]](x)*vj[[s]](x)*densq(x,MLEq)}, a, b, rel.tol=1e-8)$value
}
for (m in 1:(M1+1)){
  for (s in 1:(M1+1))
    varcov_f1[m,s]= integrate(function(x) {Kvj1[[m]](x)*Kvj1[[s]](x)*densf1(x,MLEf1)}, a, b, rel.tol=1e-8)$value
}
for (m in 1:(M1+1)){
  for (s in 1:(M1+1))
    varcov_f2[m,s]= integrate(function(x) {Kvj2[[m]](x)*Kvj2[[s]](x)*densf2(x,MLEf2)}, a, b, rel.tol=1e-8)$value
}
for (m in 1:(M1+1)){
  for (s in 1:(M1+1))
    varcov_f3[m,s]= integrate(function(x) {Kvj3[[m]](x)*Kvj3[[s]](x)*densf3(x,MLEf3)}, a, b, rel.tol=1e-8)$value
}


##################################Starting Simulation##################################

start_time21 <- Sys.time()
res_tn_boot_K2_f1 <- res_tn_order_boot_K2_f1<- numeric(B)
vn_boot_K2_f1 <- numeric(M1)
set.seed(1)
for (k in 1:B){
  indx_f1 = rmultinom(n=1,size=n,prob=c(1-0.5, 0.5))
  xx_f1k2 = c(runif(indx_f1[1],a,b), rtruncnorm(indx_f1[2], a=a, b=b, mean=MLEf1[1], sd=MLEf1[2]))
  for (j in 1:M1){
    vn_boot_K2_f1[j] <- sqrt(n)*mean(Kvj1[[j+1]](xx_f1k2))
  }
  res_tn_boot_K2_f1[k] <- max(cumsum((vn_boot_K2_f1)^2)/1:M1)
  res_tn_order_boot_K2_f1[k] <- max(cumsum(sort((vn_boot_K2_f1)^2, decreasing=TRUE)/1:M1))
  print(k)
}
end_time21 <- Sys.time()
end_time21-start_time21

start_time22 <- Sys.time()
res_tn_boot_K2_f2 <- res_tn_order_boot_K2_f2<- numeric(B)
vn_boot_K2_f2 <- numeric(M1)
set.seed(1)
for (k in 1:B){
  indx_f2 = rmultinom(n=1,size=n,prob=c(0.3, 0.7))
  xx_f2k2 = c(rtruncnorm(indx_f2[1],a=a,b=b,mean=-5,sd=1), 
              rtruncnorm(indx_f2[2],a=a,b=b,mean=MLEf2[1],sd=MLEf2[2]))
  for (j in 1:M1){
    vn_boot_K2_f2[j] <- sqrt(n)*mean(Kvj2[[j+1]](xx_f2k2))
  }
  res_tn_boot_K2_f2[k] <- max(cumsum((vn_boot_K2_f2)^2)/1:M1)
  res_tn_order_boot_K2_f2[k] <- max(cumsum(sort((vn_boot_K2_f2)^2, decreasing=TRUE)/1:M1))
  print(k)
}
end_time22 <- Sys.time()
end_time22-start_time22

start_time23 <- Sys.time()
res_tn_boot_K2_f3 <- res_tn_order_boot_K2_f3<- numeric(B)
vn_boot_K2_f3 <- numeric(M1)
set.seed(1)
for (k in 1:B){
  indx_f3 = rmultinom(n=1,size=n,prob=c(1-MLEf3[1]-MLEf3[2], MLEf3[1], MLEf3[2]))
  xx_f3k2 = c(runif(indx_f3[1],a,b), rtruncnorm(indx_f3[2],a,b,mean=-4,sd=1),
            rtrunc(indx_f3[3],spec="laplace",location=4,scale=1,a=a,b=b))
  for (j in 1:M1){
    vn_boot_K2_f3[j] <- sqrt(n)*mean(Kvj3[[j+1]](xx_f3k2))
  }
  res_tn_boot_K2_f3[k] <- max(cumsum((vn_boot_K2_f3)^2)/1:M1)
  res_tn_order_boot_K2_f3[k] <- max(cumsum(sort((vn_boot_K2_f3)^2, decreasing=TRUE)/1:M1))
  print(k)
}
end_time23 <- Sys.time()
end_time23-start_time23

save.image("Sec4-2.RData")


