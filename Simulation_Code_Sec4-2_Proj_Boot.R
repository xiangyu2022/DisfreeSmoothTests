source("C:\\Users\\zhan6\\Downloads\\hjs.R")
library(nloptr)
library(LaplacesDemon)
library(expm)
library(truncnorm)
library(orthopolynom)
library(pracma)
library(tolerance)

# Data generation
set.seed(1)
n=100; p1=0.3; p2=0.5
a=-10;b=10
indx = rmultinom(n=1,size=n,prob=c(1-p1-p2, p1,p2))
xx_q = c(runif(indx[1],a,b), rtruncnorm(indx[2],a,b,mean=-5,sd=3),
         rtrunc(indx[3],spec="laplace",location=5,scale=3,a=a,b=b))
M1 = 6

# Reference distribution is a truncated normal distribution with unknown mean and sd
densf <- function(x,pars) dtruncnorm(x,a=a,b=b,mean=pars[1],sd=pars[2])
likef <- function(pars) -sum(log(densf(xx_q,pars)))
MLEf <- constrOptim(theta=c(mean(xx_q),sd(xx_q)),f=likef,ui=c(0,1),ci=c(0),method = "Nelder-Mead")$par
G_obsf <- function(x) ptruncnorm(x,a=a,b=b,mean=MLEf[1],sd=MLEf[2])
G_obsf <- Vectorize(G_obsf)

# Postulated model G1 is a mixture of half unif[a,b], half truncated normal with unknown mean and sd
densg1 <- function(x,pars) 0.5*dunif(x,a,b) + 0.5*dtruncnorm(x,a=a,b=b,mean=pars[1],sd=pars[2])
likeg1 <- function(pars) -sum(log(densg1(xx_q,pars)))
MLEg1 <- constrOptim(theta=c(mean(xx_q),sd(xx_q)),f=likeg1,ui=c(0,1),ci=c(0),method = "Nelder-Mead")$par
G_g1 <- function(x) integrate(function(x)densg1(x,MLEg1),a,x)$value
G_g1 <- Vectorize(G_g1)

# Postulated model g2 is a Laplace distribution with unknown location and scale
densg2 <- function(x,pars) 0.3*dtruncnorm(x,a=a,b=b,mean=-5,sd=1) + 0.7*dtruncnorm(x,a=a,b=b,mean=pars[1],sd=pars[2]) 
likeg2 <- function(pars) -sum(log(densg2(xx_q,pars)))
MLEg2 <- constrOptim(theta=c(mean(xx_q),sd(xx_q)),f=likeg2,ui=c(0,1),ci=c(0),method = "Nelder-Mead")$par
G_g2 <- function(x) integrate(function(x)densg2(x,MLEg2),a,x)$value
G_g2 <- Vectorize(G_g2)

# g3 is a mixture of unif[a,b], truncated normal and truncated laplace, with the wrong parameters and unknown proportions 
densg3 <- function(x,pars) (1-pars[1]-pars[2])*dunif(x,a,b) + pars[1]*dtruncnorm(x,a,b,mean=-4,sd=1) + pars[2]*dtrunc(x,spec="laplace",location=4,scale=1,a=a,b=b)
likeg3 <- function(pars) -sum(log(densg3(xx_q,pars)))
MLEg3 <- constrOptim(theta=c(0.3,0.5), f=likeg3, ui= rbind(diag(2),-diag(2),c(1,1),c(-1,-1)),
                     ci=c(rep(0,2),rep(-1,2),0,-1), method = "Nelder-Mead", control = list(abstol=1e-15,reltol=1e-15))$par
G_g3 <- function(x) integrate(function(x)densg3(x,MLEg3),a,x)$value
G_g3 <- Vectorize(G_g3)

################## Orthornormalzied Score functions for Q, G1, G2 and G3 ##################
norm_score_f <- function(pars){
  mu<-pars[1]; sigma <- pars[2]
  grad1 <- function(x) (x-mu)/sigma^2 + (dnorm((b-mu)/sigma)/sigma-dnorm((a-mu)/sigma)/sigma)/(pnorm(b,mu,sigma)-pnorm(a,mu,sigma))
  grad2 <- function(x) (x-mu)^2/sigma^3 -1/sigma + ((b-mu)/sigma^2 *dnorm((b-mu)/sigma) - (a-mu)/sigma^2 *dnorm((a-mu)/sigma))/(pnorm(b,mu,sigma)-pnorm(a,mu,sigma))
  fish.exact <- matrix(NA, 2, 2)
  fish.exact[1, 1] <- integrate(function(x) {grad1(x)^2 *densf(x,pars)}, a, b)$value
  fish.exact[2, 2] <- integrate(function(x) {grad2(x)^2 *densf(x,pars)}, a, b)$value
  fish.exact[1, 2] <- integrate(function(x) {grad1(x)*grad2(x) *densf(x,pars)}, a, b)$value
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
norm_score_f1 = norm_score_f(MLEf)$norm_score_func1
norm_score_f2 = norm_score_f(MLEf)$norm_score_func2

norm_score_gg1 <- function(pars){
  mu<-pars[1]; sigma <- pars[2]
  grad1 <- function(x) ((x-mu)/sigma^2 + (dnorm((b-mu)/sigma)/sigma-dnorm((a-mu)/sigma)/sigma)/(pnorm(b,mu,sigma)-pnorm(a,mu,sigma)))*0.5*dtruncnorm(x,a=a,b=b,mean=mu,sd=sigma)/densg1(x,pars)
  grad2 <- function(x) ((x-mu)^2/sigma^3 -1/sigma + ((b-mu)/sigma^2 *dnorm((b-mu)/sigma) - (a-mu)/sigma^2 *dnorm((a-mu)/sigma))/(pnorm(b,mu,sigma)-pnorm(a,mu,sigma)))*0.5*dtruncnorm(x,a=a,b=b,mean=mu,sd=sigma)/densg1(x,pars)
  
  fish.exact <- matrix(NA, 2, 2)
  fish.exact[1, 1] <- integrate(function(x) {grad1(x)^2 *densg1(x,pars)}, a, b)$value
  fish.exact[2, 2] <- integrate(function(x) {grad2(x)^2 *densg1(x,pars)}, a, b)$value
  fish.exact[1, 2] <- integrate(function(x) {grad1(x)*grad2(x) *densg1(x,pars)}, a, b)$value
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
norm_score_g11 = norm_score_gg1(MLEg1)$norm_score_func1
norm_score_g12 = norm_score_gg1(MLEg1)$norm_score_func2

norm_score_gg2 <- function(pars){
  mu<-pars[1]; sigma <- pars[2]
  grad1 <- function(x) ((x-mu)/sigma^2 + (dnorm((b-mu)/sigma)/sigma-dnorm((a-mu)/sigma)/sigma)/(pnorm(b,mu,sigma)-pnorm(a,mu,sigma)))*0.7*dtruncnorm(x,a=a,b=b,mean=mu,sd=sigma)/densg2(x,pars)
  grad2 <- function(x) ((x-mu)^2/sigma^3 -1/sigma + ((b-mu)/sigma^2 *dnorm((b-mu)/sigma) - (a-mu)/sigma^2 *dnorm((a-mu)/sigma))/(pnorm(b,mu,sigma)-pnorm(a,mu,sigma)))*0.7*dtruncnorm(x,a=a,b=b,mean=mu,sd=sigma)/densg2(x,pars)
  
  fish.exact <- matrix(NA, 2, 2)
  fish.exact[1, 1] <- integrate(function(x) {grad1(x)^2 *densg2(x,pars)}, a, b, rel.tol = 1e-12)$value
  fish.exact[2, 2] <- integrate(function(x) {grad2(x)^2 *densg2(x,pars)}, a, b, rel.tol = 1e-12)$value
  fish.exact[1, 2] <- integrate(function(x) {grad1(x)*grad2(x) *densg2(x,pars)}, a, b, rel.tol = 1e-12)$value
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
norm_score_g21 = norm_score_gg2(MLEg2)$norm_score_func1
norm_score_g22 = norm_score_gg2(MLEg2)$norm_score_func2

norm_score_gg3 <- function(pars){
  grad1 <- function(x) (-dunif(x,a,b) + dtruncnorm(x,a,b,mean=-4,sd=1))/densg3(x,pars)
  grad2 <- function(x) (-dunif(x,a,b) + dtrunc(x,spec="laplace",location=4,scale=1,a=a,b=b))/densg3(x,pars)
  fish.exact <- matrix(NA, 2, 2)
  fish.exact[1, 1] <- integrate(function(x) {grad1(x)^2 *densg3(x,pars)}, a, b)$value
  fish.exact[2, 2] <- integrate(function(x) {grad2(x)^2 *densg3(x,pars)}, a, b)$value
  fish.exact[1, 2] <- integrate(function(x) {grad1(x)*grad2(x) *densg3(x,pars)}, a, b)$value
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
norm_score_g31 = norm_score_gg3(MLEg3)$norm_score_func1
norm_score_g32 = norm_score_gg3(MLEg3)$norm_score_func2


## Tilde hs for F, G1, G2, and G3 using the composite of shifted Legendre polynomials and CDFs. 
vj=list()
vj[[1]] = function(x) return (1)
for (i in 1:(M1)){
  call1 = paste("inner",i,"_1=integrate(function(x) {hj[[",i+1
                ,"]](G_obsf(x))*norm_score_f1(x)*densf(x,MLEf) }, a, b)$value",sep="")
  call2 = paste("inner",i,"_2=integrate(function(x) {hj[[",i+1
                ,"]](G_obsf(x))*norm_score_f2(x)*densf(x,MLEf) }, a, b)$value",sep="")
  call3 = paste("vj[[",i+1,"]]= function(x) {hj[[",i+1,"]](G_obsf(x)) - norm_score_f1(x) * inner",i,
                "_1 - norm_score_f2(x) * inner",i,"_2}",sep="")
  eval(parse(text=call1))
  eval(parse(text=call2))
  eval(parse(text=call3))
  print(i)
}

vj_g1=list(); vj_g1[[1]] = function(x) return (1)
for (i in 1:(M1)){
  call1 = paste("inner_g1",i,"_1=integrate(function(x) {hj[[",i+1
                ,"]](G_g1(x))*norm_score_g11(x)*densg1(x,MLEg1) }, a, b)$value",sep="")
  call2 = paste("inner_g1",i,"_2=integrate(function(x) {hj[[",i+1
                ,"]](G_g1(x))*norm_score_g12(x)*densg1(x,MLEg1) }, a, b)$value",sep="")
  call3 = paste("vj_g1[[",i+1,"]]= function(x) {hj[[",i+1,"]](G_g1(x)) - norm_score_g11(x) * inner_g1",i,
                "_1 - norm_score_g12(x) * inner_g1",i,"_2}",sep="")
  eval(parse(text=call1))
  eval(parse(text=call2))
  eval(parse(text=call3))
  print(i)
}
vj_g2=list(); vj_g2[[1]] = function(x) return (1)
for (i in 1:(M1)){
  call1 = paste("inner_g2",i,"_1=integrate(function(x) {hj[[",i+1
                ,"]](G_g2(x))*norm_score_g21(x)*densg2(x,MLEg2) }, a, b)$value",sep="")
  call2 = paste("inner_g2",i,"_2=integrate(function(x) {hj[[",i+1
                ,"]](G_g2(x))*norm_score_g22(x)*densg2(x,MLEg2) }, a, b)$value",sep="")
  call3 = paste("vj_g2[[",i+1,"]]= function(x) {hj[[",i+1,"]](G_g2(x)) - norm_score_g21(x) * inner_g2",i,
                "_1 - norm_score_g22(x) * inner_g2",i,"_2}",sep="")
  eval(parse(text=call1))
  eval(parse(text=call2))
  eval(parse(text=call3))
  print(i)
}
vj_g3=list(); vj_g3[[1]] = function(x) return (1)
for (i in 1:(M1)){
  call1 = paste("inner_g3",i,"_1=integrate(function(x) {hj[[",i+1
                ,"]](G_g3(x))*norm_score_g31(x)*densg3(x,MLEg3) }, a, b)$value",sep="")
  call2 = paste("inner_g3",i,"_2=integrate(function(x) {hj[[",i+1
                ,"]](G_g3(x))*norm_score_g32(x)*densg3(x,MLEg3) }, a, b)$value",sep="")
  call3 = paste("vj_g3[[",i+1,"]]= function(x) {hj[[",i+1,"]](G_g3(x)) - norm_score_g31(x) * inner_g3",i,
                "_1 - norm_score_g32(x) * inner_g3",i,"_2}",sep="")
  eval(parse(text=call1))
  eval(parse(text=call2))
  eval(parse(text=call3))
  print(i)
}


################### Steps before K2 ################### 
B=100000
# Starting from here, we do the bootstrap obtaining the original distributions 
res_tn_boot_f <- res_tn_order_boot_f <- numeric(B); vj_boot_f <- numeric(M1)
for (i in 1:B){
  xx_f = rtruncnorm(n,a=a,b=b,mean=MLEf[1],sd=MLEf[2])
  for (j in 1:M1){
    vj_boot_f[j] <- sqrt(n)*mean(vj[[j+1]](xx_f))
  }
  res_tn_boot_f[i] <- max(cumsum((vj_boot_f)^2)/1:M1)
  res_tn_order_boot_f[i] <- max(cumsum(sort((vj_boot_f)^2,decreasing=TRUE)/1:M1))
  print(i)
}
res_tn_boot_g1 <-  res_tn_order_boot_g1 <- numeric(B);vj_boot_g1 <- numeric(M1)
for (i in 1:B){
  indx_g1 = rmultinom(n=1,size=n,prob=c(1-0.5, 0.5))
  xx_g1 = c(runif(indx_g1[1],a,b), rtruncnorm(indx_g1[2],a=a, b=b,mean=MLEg1[1],sd=MLEg1[2]))
  for (j in 1:M1){
    vj_boot_g1[j] <- sqrt(n)*mean(vj_g1[[j+1]](xx_g1))
  }
  res_tn_boot_g1[i] <- max(cumsum((vj_boot_g1)^2)/1:M1)
  res_tn_order_boot_g1[i] <- max(cumsum(sort((vj_boot_g1)^2, decreasing=TRUE)/1:M1))
  print(i)
}
res_tn_boot_g2 <- res_tn_order_boot_g2 <- numeric(B); vj_boot_g2 <- numeric(M1)
for (i in 1:B){
  indx_g2 = rmultinom(n=1,size=n,prob=c(0.3, 0.7))
  xx_g2 = c(rtruncnorm(indx_g2[1],a=a,b=b,mean=-5,sd=1),
            rtruncnorm(indx_g2[1],a=a,b=b,mean=MLEg2[1],sd=MLEg2[2]))
  for (j in 1:M1){
    vj_boot_g2[j] <- sqrt(n)*mean(vj_g2[[j+1]](xx_g2))
  }
  res_tn_boot_g2[i] <- max(cumsum((vj_boot_g2)^2)/1:M1)
  res_tn_order_boot_g2[i] <- max(cumsum(sort((vj_boot_g2)^2, decreasing=TRUE)/1:M1))
  print(i)
}
res_tn_boot_g3 <- res_tn_order_boot_g3 <- numeric(B); vj_boot_g3 <- numeric(M1)
for (i in 1:B){
  indx_g3 = rmultinom(n=1,size=n,prob=c(1-MLEg3[1]-MLEg3[2], MLEg3[1], MLEg3[2]))
  xx_g3 = c(runif(indx_g3[1],a,b), rtruncnorm(indx_g3[2],a,b,mean=-4,sd=1),
            rtrunc(indx_g3[3],spec="laplace",location=4,scale=1,a=a,b=b))
  for (j in 1:M1){
    vj_boot_g3[j] <- sqrt(n)*mean(vj_g3[[j+1]](xx_g3))
  }
  res_tn_boot_g3[i] <- max(cumsum((vj_boot_g3)^2)/1:M1)
  res_tn_order_boot_g3[i] <- max(cumsum(sort((vj_boot_g3)^2, decreasing=TRUE)/1:M1))
  
  print(i)
}

##################################K2-STEP1##################################

### Calculating the isometry \ell which maps functions from  $L^2(F)$ to $L^2(G)$ 
isometry1=function(x) sqrt(densf(x,MLEf)/densg1(x,MLEg1))
isometry2=function(x) sqrt(densf(x,MLEf)/densg2(x,MLEg2))
isometry3=function(x) sqrt(densf(x,MLEf)/densg3(x,MLEg3))

### The isometry l applied to the original orthonormal basis h and orthonormalized score functions b
l1h=list();l2h=list();l3h=list()
for (i in 1:(M1)){
  for (j in 1:3){
    call = paste("l",j,"h[[",i,"]]<-function(x) isometry",j,"(x) * hj[[",i+1,"]](G_obsf(x))",sep="")
    eval(parse(text=call))
  }
}

l1b=list();l2b=list();l3b=list()
for (i in 1:2){
  for (j in 1:3){
    call = paste("l",j,"b[[",i,"]]<-function(x) isometry",j,"(x) *norm_score_f",i,"(x)",sep="")
    eval(parse(text=call))
  }
}

### Calculating necessary inner products for defining K
innerl1 = integrate(function(x) {sqrt(densf(x,MLEf)*densg1(x,MLEg1))}, lower=a, upper=b)$value
innerl2 = integrate(function(x) {sqrt(densf(x,MLEf)*densg2(x,MLEg2))}, lower=a, upper=b)$value
innerl3 = integrate(function(x) {sqrt(densf(x,MLEf)*densg3(x,MLEg3))}, lower=a, upper=b)$value

inner1lh = inner2lh = inner3lh = c()
for (i in 1:(M1)){
  inner1lh[i]=integrate(function(x) {l1h[[i]](x)*densg1(x,MLEg1)},lower=a,upper=b)$value
  inner2lh[i]=integrate(function(x) {l2h[[i]](x)*densg2(x,MLEg2)},lower=a,upper=b)$value
  inner3lh[i]=integrate(function(x) {l3h[[i]](x)*densg3(x,MLEg3)},lower=a,upper=b)$value
}

inner1lb = inner2lb = inner3lb = c()
for (i in 1:2){
  inner1lb[i]=integrate(function(x) {l1b[[i]](x)*densg1(x,MLEg1)},lower=a,upper=b)$value
  inner2lb[i]=integrate(function(x) {l2b[[i]](x)*densg2(x,MLEg2)},lower=a,upper=b)$value
  inner3lb[i]=integrate(function(x) {l3b[[i]](x)*densg3(x,MLEg3)},lower=a,upper=b)$value
}

### Appliying the unitary operator K to lh and lb. 
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

### Calculating necessary inner products for defining tilde c
inner_step3_1_1=integrate(function(x) {(norm_score_g11(x)-Kl1b[[1]](x))*Kl1b[[2]](x)*densg1(x,MLEg1)},lower=a,upper=b)$value
inner_step3_2_1=integrate(function(x) {norm_score_g11(x)*Kl1b[[1]](x)*densg1(x,MLEg1)},lower=a,upper=b)$value

inner_step3_1_2=integrate(function(x) {(norm_score_g21(x)-Kl2b[[1]](x))*Kl2b[[2]](x)*densg2(x,MLEg2)},lower=a,upper=b)$value
inner_step3_2_2=integrate(function(x) {norm_score_g21(x)*Kl2b[[1]](x)*densg2(x,MLEg2)},lower=a,upper=b)$value

inner_step3_1_3=integrate(function(x) {(norm_score_g31(x)-Kl3b[[1]](x))*Kl3b[[2]](x)*densg3(x,MLEg3)},lower=a,upper=b)$value
inner_step3_2_3=integrate(function(x) {norm_score_g31(x)*Kl3b[[1]](x)*densg3(x,MLEg3)},lower=a,upper=b)$value

### Calculating tilde c
tildec2_g1 = function(x){return (Kl1b[[2]](x)-inner_step3_1_1/(1-inner_step3_2_1)*(norm_score_g11(x)-Kl1b[[1]](x)))} 
tildec2_g2 = function(x){return (Kl2b[[2]](x)-inner_step3_1_2/(1-inner_step3_2_2)*(norm_score_g21(x)-Kl2b[[1]](x)))} 
tildec2_g3 = function(x){return (Kl3b[[2]](x)-inner_step3_1_3/(1-inner_step3_2_3)*(norm_score_g31(x)-Kl3b[[1]](x)))} 

### Calculating necessary inner products for defining U
inner1_1h = inner1_2h = inner1_3h = inner_U1_h = inner_U2_h = inner_U3_h = c()
U1_h = U2_h = U3_h = list()
for (i in 1:(M1)){
  inner1_1h[i]=integrate(function(x) {(norm_score_g11(x)-Kl1b[[1]](x))*Kl1h[[i]](x)*densg1(x,MLEg1)},lower=a,upper=b)$value
  inner1_2h[i]=integrate(function(x) {(norm_score_g21(x)-Kl2b[[1]](x))*Kl2h[[i]](x)*densg2(x,MLEg2)},lower=a,upper=b)$value
  inner1_3h[i]=integrate(function(x) {(norm_score_g31(x)-Kl3b[[1]](x))*Kl3h[[i]](x)*densg3(x,MLEg3)},lower=a,upper=b)$value
}
for (i in 1:(M1)){
  call003 = paste("U1_h[[",i,"]]<-function(x) {Kl1h[[",i,"]](x)-inner1_1h[",i,"]/(1-inner_step3_2_1)*(norm_score_g11(x)-Kl1b[[1]](x))}",sep="")
  call004 = paste("U2_h[[",i,"]]<-function(x) {Kl2h[[",i,"]](x)-inner1_2h[",i,"]/(1-inner_step3_2_2)*(norm_score_g21(x)-Kl2b[[1]](x))}",sep="")
  call005 = paste("U3_h[[",i,"]]<-function(x) {Kl3h[[",i,"]](x)-inner1_3h[",i,"]/(1-inner_step3_2_3)*(norm_score_g31(x)-Kl3b[[1]](x))}",sep="")
  eval(parse(text=call003))
  eval(parse(text=call004))
  eval(parse(text=call005))
  inner_U1_h[i] = integrate(function(x) {(norm_score_g12(x)-tildec2_g1(x))*U1_h[[i]](x)*densg1(x,MLEg1)},lower=a,upper=b)$value
  inner_U2_h[i] = integrate(function(x) {(norm_score_g22(x)-tildec2_g2(x))*U2_h[[i]](x)*densg2(x,MLEg2)},lower=a,upper=b)$value
  inner_U3_h[i] = integrate(function(x) {(norm_score_g32(x)-tildec2_g3(x))*U3_h[[i]](x)*densg3(x,MLEg3)},lower=a,upper=b)$value
}
inner_U1=integrate(function(x) {norm_score_g12(x)*tildec2_g1(x)*densg1(x,MLEg1)},lower=a,upper=b)$value
inner_U2=integrate(function(x) {norm_score_g22(x)*tildec2_g2(x)*densg2(x,MLEg2)},lower=a,upper=b)$value
inner_U3=integrate(function(x) {norm_score_g32(x)*tildec2_g3(x)*densg3(x,MLEg3)},lower=a,upper=b)$value

### Appliying the unitary operator U to Klh and Klb. 
UKl1h=list();UKl2h=list();UKl3h=list()
for (i in 1:(M1)){
  call01 = paste("UKl1h[[",i,"]]<-function(x) U1_h[[",i,"]](x)- inner_U1_h[",i,"]/(1-inner_U1)*(norm_score_g12(x) - tildec2_g1(x))",sep="")
  eval(parse(text=call01))
  call02 = paste("UKl2h[[",i,"]]<-function(x) U2_h[[",i,"]](x)- inner_U2_h[",i,"]/(1-inner_U2)*(norm_score_g22(x) - tildec2_g2(x))",sep="")
  eval(parse(text=call02))
  call03 = paste("UKl3h[[",i,"]]<-function(x) U3_h[[",i,"]](x)- inner_U3_h[",i,"]/(1-inner_U3)*(norm_score_g32(x) - tildec2_g3(x))",sep="")
  eval(parse(text=call03))
}


##################################Calculating tilde h or tilde psi##################################
Kvj1=list()
Kvj1[[1]] = function(x) return (1)
for (i in 1:(M1)){
  call1 = paste("inner",i,"_11=integrate(function(x) {hj[[",i+1
                ,"]](G_obsf(x))*norm_score_f1(x)*densf(x,MLEf) }, a, b)$value",sep="")
  call2 = paste("inner",i,"_12=integrate(function(x) {hj[[",i+1
                ,"]](G_obsf(x))*norm_score_f2(x)*densf(x,MLEf) }, a, b)$value",sep="")
  call3 = paste("Kvj1[[",i+1,"]]= function(x) {UKl1h[[",i,"]](x) - norm_score_g11(x) * inner",i,
                "_11 - norm_score_g12(x) * inner",i,"_12}",sep="")
  eval(parse(text=call1))
  eval(parse(text=call2))
  eval(parse(text=call3))
  print(i)
}

Kvj2=list()
Kvj2[[1]] = function(x) return (1)
for (i in 1:M1){
  call1 = paste("inner",i,"_21=integrate(function(x) {hj[[",i+1
                ,"]](G_obsf(x))*norm_score_f1(x)*densf(x,MLEf) }, a, b)$value",sep="")
  call2 =paste("inner",i,"_22=integrate(function(x) {hj[[",i+1
               ,"]](G_obsf(x))*norm_score_f2(x)*densf(x,MLEf) }, a, b)$value",sep="")
  call3 = paste("Kvj2[[",i+1,"]]= function(x) {UKl2h[[",i,"]](x) - norm_score_g21(x) * inner",i,
                "_21 - norm_score_g22(x) * inner",i,"_22}",sep="")
  eval(parse(text=call1))
  eval(parse(text=call2))
  eval(parse(text=call3))
  print(i)
}

Kvj3=list()
Kvj3[[1]] = function(x) return (1)
for (i in 1:M1){
  call1 = paste("inner",i,"_31=integrate(function(x) {hj[[",i+1
                ,"]](G_obsf(x))*norm_score_f1(x)*densf(x,MLEf) }, a, b)$value",sep="")
  call2 =paste("inner",i,"_32=integrate(function(x) {hj[[",i+1
               ,"]](G_obsf(x))*norm_score_f2(x)*densf(x,MLEf) }, a, b)$value",sep="")
  call3 = paste("Kvj3[[",i+1,"]]= function(x) {UKl3h[[",i,"]](x) - norm_score_g31(x) * inner",i,
                "_31 - norm_score_g32(x) * inner",i,"_32}",sep="")
  eval(parse(text=call1))
  eval(parse(text=call2))
  eval(parse(text=call3))
  print(i)
}

##################################Starting simulation for K2-based test statistics ##################################
start_time21 <- Sys.time()
res_tn_boot_K2_g1 <- res_tn_order_boot_K2_g1<- numeric(B)
vn_boot_K2_g1 <- numeric(M1)
set.seed(1)
for (k in 1:B){
  indx_g1 = rmultinom(n=1,size=n,prob=c(1-0.5, 0.5))
  xx_g1k2 = c(runif(indx_g1[1],a,b), rtruncnorm(indx_g1[2], a=a, b=b, mean=MLEg1[1], sd=MLEg1[2]))
  for (j in 1:M1){
    vn_boot_K2_g1[j] <- sqrt(n)*mean(Kvj1[[j+1]](xx_g1k2))
  }
  res_tn_boot_K2_g1[k] <- max(cumsum((vn_boot_K2_g1)^2)/1:M1)
  res_tn_order_boot_K2_g1[k] <- max(cumsum(sort((vn_boot_K2_g1)^2, decreasing=TRUE)/1:M1))
  print(k)
}
end_time21 <- Sys.time()
end_time21-start_time21

start_time22 <- Sys.time()
res_tn_boot_K2_g2 <- res_tn_order_boot_K2_g2<- numeric(B)
vn_boot_K2_g2 <- numeric(M1)
set.seed(1)
for (k in 1:B){
  indx_g2 = rmultinom(n=1,size=n,prob=c(0.3, 0.7))
  xx_g2k2 = c(rtruncnorm(indx_g2[1],a=a,b=b,mean=-5,sd=1), 
              rtruncnorm(indx_g2[2],a=a,b=b,mean=MLEg2[1],sd=MLEg2[2]))
  for (j in 1:M1){
    vn_boot_K2_g2[j] <- sqrt(n)*mean(Kvj2[[j+1]](xx_g2k2))
  }
  res_tn_boot_K2_g2[k] <- max(cumsum((vn_boot_K2_g2)^2)/1:M1)
  res_tn_order_boot_K2_g2[k] <- max(cumsum(sort((vn_boot_K2_g2)^2, decreasing=TRUE)/1:M1))
  print(k)
}
end_time22 <- Sys.time()
end_time22-start_time22

start_time23 <- Sys.time()
res_tn_boot_K2_g3 <- res_tn_order_boot_K2_g3<- numeric(B)
vn_boot_K2_g3 <- numeric(M1)
set.seed(1)
for (k in 1:B){
  indx_g3 = rmultinom(n=1,size=n,prob=c(1-MLEg3[1]-MLEg3[2], MLEg3[1], MLEg3[2]))
  xx_g3k2 = c(runif(indx_g3[1],a,b), rtruncnorm(indx_g3[2],a,b,mean=-4,sd=1),
              rtrunc(indx_g3[3],spec="laplace",location=4,scale=1,a=a,b=b))
  for (j in 1:M1){
    vn_boot_K2_g3[j] <- sqrt(n)*mean(Kvj3[[j+1]](xx_g3k2))
  }
  res_tn_boot_K2_g3[k] <- max(cumsum((vn_boot_K2_g3)^2)/1:M1)
  res_tn_order_boot_K2_g3[k] <- max(cumsum(sort((vn_boot_K2_g3)^2, decreasing=TRUE)/1:M1))
  print(k)
}
end_time23 <- Sys.time()
end_time23-start_time23

save.image("Sec4-2.RData")

#################################### Power #################################### 
set.seed(1)
power_f = power_order_f = power_g1 = power_order_g1 = power_g2 = power_order_g2 = power_g3 = power_order_g3 = c()
power_K2_g1 = power_order_K2_g1 = power_K2_g2 = power_order_K2_g2 = power_K2_g3 = power_order_K2_g3 = c()

for (o in 1:10000){
  indx = rmultinom(n=1,size=n,prob=c(1-p1-p2, p1,p2))
  xx_q = c(runif(indx[1],a,b), rtruncnorm(indx[2],a,b,mean=-5,sd=3),
           rtrunc(indx[3],spec="laplace",location=5,scale=3,a=a,b=b))
  
  # Reference distribution is a truncated normal distribution with unknown mean and sd
  densf <- function(x,pars) dtruncnorm(x,a=a,b=b,mean=pars[1],sd=pars[2])
  likef <- function(pars) -sum(log(densf(xx_q,pars)))
  MLEf <- constrOptim(theta=c(mean(xx_q),sd(xx_q)),f=likef,ui=c(0,1),ci=c(0),method = "Nelder-Mead")$par
  G_obsf = function(x) ptruncnorm(x,a=a,b=b,mean=MLEf[1],sd=MLEf[2])
  G_obsf <- Vectorize(G_obsf)
  
  # Postulated model g1 is a mixture of half unif[a,b], half truncated normal with unknown mean and sd
  densg1 <- function(x,pars) 0.5*dunif(x,a,b) + 0.5*dtruncnorm(x,a=a,b=b,mean=pars[1],sd=pars[2])
  likeg1 <- function(pars) -sum(log(densg1(xx_q,pars)))
  MLEg1 <- constrOptim(theta=c(mean(xx_q),sd(xx_q)),f=likeg1,ui=c(0,1),ci=c(0),method = "Nelder-Mead")$par
  G_g1 <- function(x) integrate(function(x)densg1(x,MLEg1),a,x)$value
  G_g1 <- Vectorize(G_g1)
  
  # Postulated model g2 is a Laplace distribution with unknown location and scale
  densg2 <- function(x,pars) 0.3*dtruncnorm(x,a=a,b=b,mean=-5,sd=1) + 0.7*dtruncnorm(x,a=a,b=b,mean=pars[1],sd=pars[2]) 
  likeg2 <- function(pars) -sum(log(densg2(xx_q,pars)))
  MLEg2 <- constrOptim(theta=c(mean(xx_q),sd(xx_q)),f=likeg2,ui=c(0,1),ci=c(0),method = "Nelder-Mead")$par
  G_g2 <- function(x) integrate(function(x)densg2(x,MLEg2),a,x)$value
  G_g2 <- Vectorize(G_g2)
  
  # g3 is a mixture of unif[a,b], truncated normal and truncated laplace, with the wrong parameters and unknown proportions 
  densg3 <- function(x,pars) (1-pars[1]-pars[2])*dunif(x,a,b) + pars[1]*dtruncnorm(x,a,b,mean=-4,sd=1) + pars[2]*dtrunc(x,spec="laplace",location=4,scale=1,a=a,b=b)
  likeg3 <- function(pars) -sum(log(densg3(xx_q,pars)))
  MLEg3 <- constrOptim(theta=c(0.3,0.5), f=likeg3, ui= rbind(diag(2),-diag(2),c(1,1),c(-1,-1)),
                       ci=c(rep(0,2),rep(-1,2),0,-1), method = "Nelder-Mead", control = list(abstol=1e-15,reltol=1e-15))$par
  G_g3 <- function(x) integrate(function(x)densg3(x,MLEg3),a,x)$value
  G_g3 <- Vectorize(G_g3)
  
  norm_score_f <- function(pars){
    mu<-pars[1]; sigma <- pars[2]
    grad1 <- function(x) (x-mu)/sigma^2 + (dnorm((b-mu)/sigma)/sigma-dnorm((a-mu)/sigma)/sigma)/(pnorm(b,mu,sigma)-pnorm(a,mu,sigma))
    grad2 <- function(x) (x-mu)^2/sigma^3 -1/sigma + ((b-mu)/sigma^2 *dnorm((b-mu)/sigma) - (a-mu)/sigma^2 *dnorm((a-mu)/sigma))/(pnorm(b,mu,sigma)-pnorm(a,mu,sigma))
    fish.exact <- matrix(NA, 2, 2)
    fish.exact[1, 1] <- integrate(function(x) {grad1(x)^2 *densf(x,pars)}, a, b)$value
    fish.exact[2, 2] <- integrate(function(x) {grad2(x)^2 *densf(x,pars)}, a, b)$value
    fish.exact[1, 2] <- integrate(function(x) {grad1(x)*grad2(x) *densf(x,pars)}, a, b)$value
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
  norm_score_f1 = norm_score_f(MLEf)$norm_score_func1
  norm_score_f2 = norm_score_f(MLEf)$norm_score_func2
  
  norm_score_gg1 <- function(pars){
    mu<-pars[1]; sigma <- pars[2]
    grad1 <- function(x) ((x-mu)/sigma^2 + (dnorm((b-mu)/sigma)/sigma-dnorm((a-mu)/sigma)/sigma)/(pnorm(b,mu,sigma)-pnorm(a,mu,sigma)))*0.5*dtruncnorm(x,a=a,b=b,mean=mu,sd=sigma)/densg1(x,pars)
    grad2 <- function(x) ((x-mu)^2/sigma^3 -1/sigma + ((b-mu)/sigma^2 *dnorm((b-mu)/sigma) - (a-mu)/sigma^2 *dnorm((a-mu)/sigma))/(pnorm(b,mu,sigma)-pnorm(a,mu,sigma)))*0.5*dtruncnorm(x,a=a,b=b,mean=mu,sd=sigma)/densg1(x,pars)
    
    fish.exact <- matrix(NA, 2, 2)
    fish.exact[1, 1] <- integrate(function(x) {grad1(x)^2 *densg1(x,pars)}, a, b)$value
    fish.exact[2, 2] <- integrate(function(x) {grad2(x)^2 *densg1(x,pars)}, a, b)$value
    fish.exact[1, 2] <- integrate(function(x) {grad1(x)*grad2(x) *densg1(x,pars)}, a, b)$value
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
  norm_score_g11 = norm_score_gg1(MLEg1)$norm_score_func1
  norm_score_g12 = norm_score_gg1(MLEg1)$norm_score_func2
  
  norm_score_gg2 <- function(pars){
    mu<-pars[1]; sigma <- pars[2]
    grad1 <- function(x) ((x-mu)/sigma^2 + (dnorm((b-mu)/sigma)/sigma-dnorm((a-mu)/sigma)/sigma)/(pnorm(b,mu,sigma)-pnorm(a,mu,sigma)))*0.7*dtruncnorm(x,a=a,b=b,mean=mu,sd=sigma)/densg2(x,pars)
    grad2 <- function(x) ((x-mu)^2/sigma^3 -1/sigma + ((b-mu)/sigma^2 *dnorm((b-mu)/sigma) - (a-mu)/sigma^2 *dnorm((a-mu)/sigma))/(pnorm(b,mu,sigma)-pnorm(a,mu,sigma)))*0.7*dtruncnorm(x,a=a,b=b,mean=mu,sd=sigma)/densg2(x,pars)
    
    fish.exact <- matrix(NA, 2, 2)
    fish.exact[1, 1] <- integrate(function(x) {grad1(x)^2 *densg2(x,pars)}, a, b, rel.tol = 1e-12)$value
    fish.exact[2, 2] <- integrate(function(x) {grad2(x)^2 *densg2(x,pars)}, a, b, rel.tol = 1e-12)$value
    fish.exact[1, 2] <- integrate(function(x) {grad1(x)*grad2(x) *densg2(x,pars)}, a, b, rel.tol = 1e-12)$value
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
  norm_score_g21 = norm_score_gg2(MLEg2)$norm_score_func1
  norm_score_g22 = norm_score_gg2(MLEg2)$norm_score_func2
  
  norm_score_gg3 <- function(pars){
    grad1 <- function(x) (-dunif(x,a,b) + dtruncnorm(x,a,b,mean=-4,sd=1))/densg3(x,pars)
    grad2 <- function(x) (-dunif(x,a,b) + dtrunc(x,spec="laplace",location=4,scale=1,a=a,b=b))/densg3(x,pars)
    fish.exact <- matrix(NA, 2, 2)
    fish.exact[1, 1] <- integrate(function(x) {grad1(x)^2 *densg3(x,pars)}, a, b)$value
    fish.exact[2, 2] <- integrate(function(x) {grad2(x)^2 *densg3(x,pars)}, a, b)$value
    fish.exact[1, 2] <- integrate(function(x) {grad1(x)*grad2(x) *densg3(x,pars)}, a, b)$value
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
  norm_score_g31 = norm_score_gg3(MLEg3)$norm_score_func1
  norm_score_g32 = norm_score_gg3(MLEg3)$norm_score_func2
  
  vj = list(); vj1=list(); vj2=list(); vj3=list()
  vj[[1]] = function(x) return (1); vj1[[1]] = function(x) return (1); vj2[[1]] = function(x) return (1); vj3[[1]] = function(x) return (1)
  
  for (i in 1:(M1)){
    call1 = paste("inner",i,"_1=integrate(function(x) {hj[[",i+1
                  ,"]](G_obsf(x))*norm_score_f1(x)*densf(x,MLEf) }, a, b)$value",sep="")
    call2 = paste("inner",i,"_2=integrate(function(x) {hj[[",i+1
                  ,"]](G_obsf(x))*norm_score_f2(x)*densf(x,MLEf) }, a, b)$value",sep="")
    call3 = paste("vj[[",i+1,"]]= function(x) {hj[[",i+1,"]](G_obsf(x)) - norm_score_f1(x) * inner",i,
                  "_1 - norm_score_f2(x) * inner",i,"_2}",sep="")
    eval(parse(text=call1))
    eval(parse(text=call2))
    eval(parse(text=call3))
    call1 = paste("inner_g1",i,"_1=integrate(function(x) {hj[[",i+1
                  ,"]](G_g1(x))*norm_score_g11(x)*densg1(x,MLEg1) }, a, b)$value",sep="")
    call2 = paste("inner_g1",i,"_2=integrate(function(x) {hj[[",i+1
                  ,"]](G_g1(x))*norm_score_g12(x)*densg1(x,MLEg1) }, a, b)$value",sep="")
    call3 = paste("vj_g1[[",i+1,"]]= function(x) {hj[[",i+1,"]](G_g1(x)) - norm_score_g11(x) * inner_g1",i,
                  "_1 - norm_score_g12(x) * inner_g1",i,"_2}",sep="")
    eval(parse(text=call1))
    eval(parse(text=call2))
    eval(parse(text=call3))
    call1 = paste("inner_g2",i,"_1=integrate(function(x) {hj[[",i+1
                  ,"]](G_g2(x))*norm_score_g21(x)*densg2(x,MLEg2) }, a, b)$value",sep="")
    call2 = paste("inner_g2",i,"_2=integrate(function(x) {hj[[",i+1
                  ,"]](G_g2(x))*norm_score_g22(x)*densg2(x,MLEg2) }, a, b)$value",sep="")
    call3 = paste("vj_g2[[",i+1,"]]= function(x) {hj[[",i+1,"]](G_g2(x)) - norm_score_g21(x) * inner_g2",i,
                  "_1 - norm_score_g22(x) * inner_g2",i,"_2}",sep="")
    eval(parse(text=call1))
    eval(parse(text=call2))
    eval(parse(text=call3))
    call1 = paste("inner_g3",i,"_1=integrate(function(x) {hj[[",i+1
                  ,"]](G_g3(x))*norm_score_g31(x)*densg3(x,MLEg3) }, a, b)$value",sep="")
    call2 = paste("inner_g3",i,"_2=integrate(function(x) {hj[[",i+1
                  ,"]](G_g3(x))*norm_score_g32(x)*densg3(x,MLEg3) }, a, b)$value",sep="")
    call3 = paste("vj_g3[[",i+1,"]]= function(x) {hj[[",i+1,"]](G_g3(x)) - norm_score_g31(x) * inner_g3",i,
                  "_1 - norm_score_g32(x) * inner_g3",i,"_2}",sep="")
    eval(parse(text=call1))
    eval(parse(text=call2))
    eval(parse(text=call3))
  }
  
  ##################################K2-STEP1##################################
  isometry1=function(x) sqrt(densf(x,MLEf)/densg1(x,MLEg1))
  isometry2=function(x) sqrt(densf(x,MLEf)/densg2(x,MLEg2))
  isometry3=function(x) sqrt(densf(x,MLEf)/densg3(x,MLEg3))
  
  
  l1h=list();l2h=list();l3h=list()
  for (i in 1:(M1)){
    for (j in 1:3){
      call = paste("l",j,"h[[",i,"]]<-function(x) isometry",j,"(x) * hj[[",i+1,"]](G_obsf(x))",sep="")
      eval(parse(text=call))
    }
  }
  
  l1b=list();l2b=list();l3b=list()
  for (i in 1:2){
    for (j in 1:3){
      call = paste("l",j,"b[[",i,"]]<-function(x) isometry",j,"(x) *norm_score_f",i,"(x)",sep="")
      eval(parse(text=call))
    }
  }
  
  innerl1 = integrate(function(x) {sqrt(densf(x,MLEf)*densg1(x,MLEg1))}, lower=a, upper=b)$value
  innerl2 = integrate(function(x) {sqrt(densf(x,MLEf)*densg2(x,MLEg2))}, lower=a, upper=b)$value
  innerl3 = integrate(function(x) {sqrt(densf(x,MLEf)*densg3(x,MLEg3))}, lower=a, upper=b)$value
  
  inner1lh = inner2lh = inner3lh = c()
  for (i in 1:(M1)){
    inner1lh[i]=integrate(function(x) {l1h[[i]](x)*densg1(x,MLEg1)},lower=a,upper=b)$value
    inner2lh[i]=integrate(function(x) {l2h[[i]](x)*densg2(x,MLEg2)},lower=a,upper=b)$value
    inner3lh[i]=integrate(function(x) {l3h[[i]](x)*densg3(x,MLEg3)},lower=a,upper=b)$value
  }
  
  inner1lb = inner2lb = inner3lb = c()
  for (i in 1:2){
    inner1lb[i]=integrate(function(x) {l1b[[i]](x)*densg1(x,MLEg1)},lower=a,upper=b)$value
    inner2lb[i]=integrate(function(x) {l2b[[i]](x)*densg2(x,MLEg2)},lower=a,upper=b)$value
    inner3lb[i]=integrate(function(x) {l3b[[i]](x)*densg3(x,MLEg3)},lower=a,upper=b)$value
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
  
  inner_step3_1_1=integrate(function(x) {(norm_score_g11(x)-Kl1b[[1]](x))*Kl1b[[2]](x)*densg1(x,MLEg1)},lower=a,upper=b)$value
  inner_step3_2_1=integrate(function(x) {norm_score_g11(x)*Kl1b[[1]](x)*densg1(x,MLEg1)},lower=a,upper=b)$value
  
  inner_step3_1_2=integrate(function(x) {(norm_score_g21(x)-Kl2b[[1]](x))*Kl2b[[2]](x)*densg2(x,MLEg2)},lower=a,upper=b)$value
  inner_step3_2_2=integrate(function(x) {norm_score_g21(x)*Kl2b[[1]](x)*densg2(x,MLEg2)},lower=a,upper=b)$value
  
  inner_step3_1_3=integrate(function(x) {(norm_score_g31(x)-Kl3b[[1]](x))*Kl3b[[2]](x)*densg3(x,MLEg3)},lower=a,upper=b)$value
  inner_step3_2_3=integrate(function(x) {norm_score_g31(x)*Kl3b[[1]](x)*densg3(x,MLEg3)},lower=a,upper=b)$value
  
  
  tildec2_g1 = function(x){return (Kl1b[[2]](x)-inner_step3_1_1/(1-inner_step3_2_1)*(norm_score_g11(x)-Kl1b[[1]](x)))} # Same with below
  tildec2_g2 = function(x){return (Kl2b[[2]](x)-inner_step3_1_2/(1-inner_step3_2_2)*(norm_score_g21(x)-Kl2b[[1]](x)))} # Same with below
  tildec2_g3 = function(x){return (Kl3b[[2]](x)-inner_step3_1_3/(1-inner_step3_2_3)*(norm_score_g31(x)-Kl3b[[1]](x)))} # Same with below
  
  
  inner1_1h = inner1_2h = inner1_3h = inner_U1_h = inner_U2_h = inner_U3_h = c()
  U1_h = U2_h = U3_h = list()
  for (i in 1:(M1)){
    inner1_1h[i]=integrate(function(x) {(norm_score_g11(x)-Kl1b[[1]](x))*Kl1h[[i]](x)*densg1(x,MLEg1)},lower=a,upper=b)$value
    inner1_2h[i]=integrate(function(x) {(norm_score_g21(x)-Kl2b[[1]](x))*Kl2h[[i]](x)*densg2(x,MLEg2)},lower=a,upper=b)$value
    inner1_3h[i]=integrate(function(x) {(norm_score_g31(x)-Kl3b[[1]](x))*Kl3h[[i]](x)*densg3(x,MLEg3)},lower=a,upper=b)$value
  }
  for (i in 1:(M1)){
    call003 = paste("U1_h[[",i,"]]<-function(x) {Kl1h[[",i,"]](x)-inner1_1h[",i,"]/(1-inner_step3_2_1)*(norm_score_g11(x)-Kl1b[[1]](x))}",sep="")
    call004 = paste("U2_h[[",i,"]]<-function(x) {Kl2h[[",i,"]](x)-inner1_2h[",i,"]/(1-inner_step3_2_2)*(norm_score_g21(x)-Kl2b[[1]](x))}",sep="")
    call005 = paste("U3_h[[",i,"]]<-function(x) {Kl3h[[",i,"]](x)-inner1_3h[",i,"]/(1-inner_step3_2_3)*(norm_score_g31(x)-Kl3b[[1]](x))}",sep="")
    eval(parse(text=call003))
    eval(parse(text=call004))
    eval(parse(text=call005))
    inner_U1_h[i] = integrate(function(x) {(norm_score_g12(x)-tildec2_g1(x))*U1_h[[i]](x)*densg1(x,MLEg1)},lower=a,upper=b)$value
    inner_U2_h[i] = integrate(function(x) {(norm_score_g22(x)-tildec2_g2(x))*U2_h[[i]](x)*densg2(x,MLEg2)},lower=a,upper=b)$value
    inner_U3_h[i] = integrate(function(x) {(norm_score_g32(x)-tildec2_g3(x))*U3_h[[i]](x)*densg3(x,MLEg3)},lower=a,upper=b)$value
  }
  inner_U1=integrate(function(x) {norm_score_g12(x)*tildec2_g1(x)*densg1(x,MLEg1)},lower=a,upper=b)$value
  inner_U2=integrate(function(x) {norm_score_g22(x)*tildec2_g2(x)*densg2(x,MLEg2)},lower=a,upper=b)$value
  inner_U3=integrate(function(x) {norm_score_g32(x)*tildec2_g3(x)*densg3(x,MLEg3)},lower=a,upper=b)$value
  
  UKl1h=list();UKl2h=list();UKl3h=list()
  for (i in 1:(M1)){
    call01 = paste("UKl1h[[",i,"]]<-function(x) U1_h[[",i,"]](x)- inner_U1_h[",i,"]/(1-inner_U1)*(norm_score_g12(x) - tildec2_g1(x))",sep="")
    eval(parse(text=call01))
    call02 = paste("UKl2h[[",i,"]]<-function(x) U2_h[[",i,"]](x)- inner_U2_h[",i,"]/(1-inner_U2)*(norm_score_g22(x) - tildec2_g2(x))",sep="")
    eval(parse(text=call02))
    call03 = paste("UKl3h[[",i,"]]<-function(x) U3_h[[",i,"]](x)- inner_U3_h[",i,"]/(1-inner_U3)*(norm_score_g32(x) - tildec2_g3(x))",sep="")
    eval(parse(text=call03))
  }
  
  
  ##################################Calculating tilde h or tilde psi##################################
  Kvj1=list(); Kvj2=list();Kvj3=list()
  Kvj1[[1]] = function(x) return (1); Kvj2[[1]] = function(x) return (1); Kvj3[[1]] = function(x) return (1)
  for (i in 1:(M1)){
    call1 = paste("inner",i,"_11=integrate(function(x) {hj[[",i+1
                  ,"]](G_obsf(x))*norm_score_f1(x)*densf(x,MLEf) }, a, b)$value",sep="")
    call2 = paste("inner",i,"_12=integrate(function(x) {hj[[",i+1
                  ,"]](G_obsf(x))*norm_score_f2(x)*densf(x,MLEf) }, a, b)$value",sep="")
    call3 = paste("Kvj1[[",i+1,"]]= function(x) {UKl1h[[",i,"]](x) - norm_score_g11(x) * inner",i,
                  "_11 - norm_score_g12(x) * inner",i,"_12}",sep="")
    eval(parse(text=call1))
    eval(parse(text=call2))
    eval(parse(text=call3))
    
    call1 = paste("inner",i,"_21=integrate(function(x) {hj[[",i+1
                  ,"]](G_obsf(x))*norm_score_f1(x)*densf(x,MLEf) }, a, b)$value",sep="")
    call2 =paste("inner",i,"_22=integrate(function(x) {hj[[",i+1
                 ,"]](G_obsf(x))*norm_score_f2(x)*densf(x,MLEf) }, a, b)$value",sep="")
    call3 = paste("Kvj2[[",i+1,"]]= function(x) {UKl2h[[",i,"]](x) - norm_score_g21(x) * inner",i,
                  "_21 - norm_score_g22(x) * inner",i,"_22}",sep="")
    eval(parse(text=call1))
    eval(parse(text=call2))
    eval(parse(text=call3))
    
    call1 = paste("inner",i,"_31=integrate(function(x) {hj[[",i+1
                  ,"]](G_obsf(x))*norm_score_f1(x)*densf(x,MLEf) }, a, b)$value",sep="")
    call2 =paste("inner",i,"_32=integrate(function(x) {hj[[",i+1
                 ,"]](G_obsf(x))*norm_score_f2(x)*densf(x,MLEf) }, a, b)$value",sep="")
    call3 = paste("Kvj3[[",i+1,"]]= function(x) {UKl3h[[",i,"]](x) - norm_score_g31(x) * inner",i,
                  "_31 - norm_score_g32(x) * inner",i,"_32}",sep="")
    eval(parse(text=call1))
    eval(parse(text=call2))
    eval(parse(text=call3))
  }
  ########## Test statistics   ########## 
  for (j in 1:M1){
    vj_boot_f[j] <- sqrt(n)*mean(vj[[j+1]](xx_q))
    vj_boot_g1[j] <- sqrt(n)*mean(vj_g1[[j+1]](xx_q))
    vj_boot_g2[j] <- sqrt(n)*mean(vj_g2[[j+1]](xx_q))
    vj_boot_g3[j] <- sqrt(n)*mean(vj_g3[[j+1]](xx_q))
    vn_boot_K2_g1[j] <- sqrt(n)*mean(Kvj1[[j+1]](xx_q))
    vn_boot_K2_g2[j] <- sqrt(n)*mean(Kvj2[[j+1]](xx_q))
    vn_boot_K2_g3[j] <- sqrt(n)*mean(Kvj3[[j+1]](xx_q))
  }
  res_f <- max(cumsum((vj_boot_f)^2)/1:M1)
  res_order_f <- max(cumsum(sort((vj_boot_f)^2,decreasing=TRUE)/1:M1))
  
  res_boot_g1 <- max(cumsum((vj_boot_g1)^2)/1:M1)
  res_order_boot_g1 <- max(cumsum(sort((vj_boot_g1)^2, decreasing=TRUE)/1:M1))
  
  res_boot_g2 <- max(cumsum((vj_boot_g2)^2)/1:M1)
  res_order_boot_g2 <- max(cumsum(sort((vj_boot_g2)^2, decreasing=TRUE)/1:M1))
  
  res_boot_g3 <- max(cumsum((vj_boot_g3)^2)/1:M1)
  res_order_boot_g3 <- max(cumsum(sort((vj_boot_g3)^2, decreasing=TRUE)/1:M1))
  
  res_boot_K2_g1 <- max(cumsum((vn_boot_K2_g1)^2)/1:M1)
  res_order_boot_K2_g1 <- max(cumsum(sort((vn_boot_K2_g1)^2,decreasing=TRUE))/1:M1)
  
  res_boot_K2_g2 <- max(cumsum((vn_boot_K2_g2)^2)/1:M1)
  res_order_boot_K2_g2 <- max(cumsum(sort((vn_boot_K2_g2)^2, decreasing=TRUE)/1:M1))
  
  res_boot_K2_g3 <- max(cumsum((vn_boot_K2_g3)^2)/1:M1)
  res_order_boot_K2_g3 <- max(cumsum(sort((vn_boot_K2_g3)^2, decreasing=TRUE)/1:M1))
  
  power_f[o] = (sum(res_tn_boot_f>= res_f)+1) / (B+1)
  power_order_f[o] = (sum(res_tn_order_boot_f>= res_order_f)+1) / (B+1)
  
  power_g1[o] = (sum(res_tn_boot_g1>= res_boot_g1)+1) / (B+1)
  power_order_g1[o] = (sum(res_tn_order_boot_g1>= res_order_boot_g1)+1) / (B+1)
  
  power_g2[o] = (sum(res_tn_boot_g2>= res_boot_g2)+1) / (B+1)
  power_order_g2[o] = (sum(res_tn_order_boot_g2>= res_order_boot_g2)+1) / (B+1)
  
  power_g3[o] = (sum(res_tn_boot_g3>= res_boot_g3)+1) / (B+1)
  power_order_g3[o] = (sum(res_tn_order_boot_g3>= res_order_boot_g3)+1) / (B+1)
  
  power_K2_g1[o] = (sum(res_tn_boot_K2_g1>= res_boot_K2_g1)+1) / (B+1)
  power_order_K2_g1[o] = (sum(res_tn_order_boot_K2_g1>= res_order_boot_K2_g1)+1) / (B+1)
  
  power_K2_g2[o] = (sum(res_tn_boot_K2_g2>= res_boot_K2_g2)+1) / (B+1)
  power_order_K2_g2[o] = (sum(res_tn_order_boot_K2_g2>= res_order_boot_K2_g2)+1) / (B+1)
  
  power_K2_g3[o] = (sum(res_tn_boot_K2_g3>= res_boot_K2_g3)+1) / (B+1)
  power_order_K2_g3[o] = (sum(res_tn_order_boot_K2_g3>= res_order_boot_K2_g3)+1) / (B+1)
  
  print(o)
}

mean(power_f<0.001); mean(power_order_f<0.001);mean(power_f<0.05); mean(power_order_f<0.05); mean(power_f<0.1); mean(power_order_f<0.1)
mean(power_g1<0.001); mean(power_order_g1<0.001); mean(power_K2_g1<0.001); mean(power_order_K2_g1<0.001) 
mean(power_g1<0.05); mean(power_order_g1<0.05); mean(power_K2_g1<0.05); mean(power_order_K2_g1<0.05)
mean(power_g1<0.1); mean(power_order_g1<0.1); mean(power_K2_g1<0.1); mean(power_order_K2_g1<0.1)
mean(power_g2<0.001); mean(power_order_g2<0.001); mean(power_K2_g2<0.001); mean(power_order_K2_g2<0.001) 
mean(power_g2<0.05); mean(power_order_g2<0.05); mean(power_K2_g2<0.05); mean(power_order_K2_g2<0.05)
mean(power_g2<0.1); mean(power_order_g2<0.1); mean(power_K2_g2<0.1); mean(power_order_K2_g2<0.1)
mean(power_g3<0.001); mean(power_order_g3<0.001); mean(power_K2_g3<0.001); mean(power_order_K2_g3<0.001) 
mean(power_g3<0.05); mean(power_order_g3<0.05); mean(power_K2_g3<0.05); mean(power_order_K2_g3<0.05)
mean(power_g3<0.1); mean(power_order_g3<0.1); mean(power_K2_g3<0.1); mean(power_order_K2_g3<0.1)




