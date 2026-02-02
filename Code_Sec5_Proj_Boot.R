source("C:\\Users\\zhan6\\Downloads\\hjs.R")
library(univariateML)
library(nloptr)
library(LaplacesDemon)
library(expm)
library(truncnorm)
library(LPBkg)
library(orthopolynom)
library(pracma)
library(cascsim)
library(tolerance)
library(extraDistr)
bs1<-read.table("C:\\Users\\zhan6\\Downloads\\16688_src.txt")
bs2<-read.table("C:\\Users\\zhan6\\Downloads\\18710_src.txt")
bs_full<-c(bs1[,1],bs2[,1])
bs_abs<-abs(bs_full)
L1=1.65;U1=2.05;
bs11<-bs_abs[bs_abs>=L1&bs_abs<=U1]
n=length(bs11)

## Specify the postulated distribution G 
set.seed(1)
bkg1<-function(x) ifelse(x>=L1&x<=U1, 1/(U1-L1), 0)
pos = c(1.78499, 1.85247, 1.94365);sigma=0.0025
new_fe = function(x,c) ifelse(x>=L1+0.015-c & x<=U1-0.015-c,dlst(x, df=4, sigma=0.025)/(plst(U1-0.015-c,df=4,sigma=0.025)-plst(L1+0.015-c, df=4, sigma=0.025)),0)
line_int<-function(w0,x,c){dtruncnorm(w0,a=c-0.015,b=c+0.015,c,sigma)*new_fe(x-w0,c)}
fw<-function(x,c){ifelse(x>=L1&x<=U1,integrate(line_int,lower=c-0.015,upper=c+0.015,x=x,c=c,rel.tol=1e-8)$value,0)}
fw_norm<-Vectorize(fw)

new_fe_dist = function(x,c) (plst(x, df=4, sigma=0.025)-plst(L1+0.015-c,df=4, sigma=0.025))/(plst(U1-0.015-c, df=4, sigma=0.025)-plst(L1+0.015-c, df=4, sigma=0.025))
dist_int<-function(w0,x,c){dtruncnorm(w0,a=c-0.015,b=c+0.015,c,sigma)*new_fe_dist(x-w0,c)}
fw_dist<-function(x,c){integrate(dist_int,lower=c-0.015,upper=c+0.015,x=x,c=c,rel.tol = 1e-8)$value}
fw_dist<-Vectorize(fw_dist)


densg1 <- function(x,pars){
  (1-sum(pars[1:3]))*bkg1(x)+pars[1]*fw_norm(x,pos[1])+pars[2]*fw_norm(x,pos[2])+pars[3]*fw_norm(x,pos[3])
}
gradg1 <- function(pars){
  grad1 <- function(x) (-bkg1(x) + fw_norm(x,pos[1]))/densg1(x,pars)
  grad2 <- function(x) (-bkg1(x) + fw_norm(x,pos[2]))/densg1(x,pars)
  grad3 <- function(x) (-bkg1(x) + fw_norm(x,pos[3]))/densg1(x,pars)
  return(c(-sum(grad1(bs11)), -sum(grad2(bs11)), -sum(grad3(bs11))))
}

## Solve the MLE for G
start_time1 <- Sys.time()
ll_bs_new1<-function(pars){-sum(log(densg1(bs11,pars)))}
MLEg1 = constrOptim(theta=c(0.05,0.05,0.2), grad=gradg1, method ="BFGS",
                    f=ll_bs_new1, ui= rbind(diag(3),-diag(3),c(1,1,1),c(-1,-1,-1)),
                    ci=c(rep(0,3),rep(-1,3),0,-1))$par
end_time1 <- Sys.time()
end_time1 - start_time1

g_signal_nor1 <- function(x) ifelse(x>=1.65&x<=2.05, densg1(x,MLEg1),0)
integrate(g_signal_nor1,1.65,2.05)

## Specify the reference distribution F and solve for its MLE 
start_time2 <- Sys.time()
densf <- function(x,pars){
  (1-sum(pars[1:3]))*bkg1(x)+pars[1]*dtruncnorm(x,a=L1,b=U1,mean=pos[1],sd=0.05)+
    pars[2]*dtruncnorm(x,a=L1,b=U1,mean=pos[2],sd=0.05)+pars[3]*dtruncnorm(x,a=L1,b=U1,mean=pos[3],sd=0.05)
}
likef <- function(pars) -sum(log(densf(bs11,pars)))
MLEf = constrOptim(theta=c(0.05,0.05,0.2), grad=NULL,
                   f=likef, ui= rbind(diag(3),-diag(3),c(1,1,1),c(-1,-1,-1)),
                   ci=c(rep(0,3),rep(-1,3),0,-1))$par
end_time2 <- Sys.time()
end_time2-start_time2

f_signal <- function(x) ifelse(x>=1.65&x<=2.05, densf(x,MLEf),0)
f_signal_nor <- function(x) f_signal(x)/integrate(f_signal,L1,U1)$value
G_f<-function(x,pars){
  res = (1-sum(pars[1:3]))*(2.5*x-4.125)+pars[1]*ptruncnorm(x,L1,U1,pos[1],sd=0.05)+pars[2]*ptruncnorm(x,L1,U1,pos[2],sd=0.05)+pars[3]*ptruncnorm(x,L1,U1,pos[3],sd=0.05)
  return (res)
}

## Plot their densities
hist(bs11,prob=T,xlim=c(1.65,2.05),ylim=c(0,5),breaks =25)
curve(g_signal_nor1,add=T)
curve(f_signal_nor,add=T,col="red")

### Specifying the orthonormalzied score functions for G and F
norm_score_g1 <- function(pars){
  grad1 <- function(x) (-bkg1(x) + fw_norm(x,pos[1]))/densg1(x,pars)
  grad2 <- function(x) (-bkg1(x) + fw_norm(x,pos[2]))/densg1(x,pars)
  grad3 <- function(x) (-bkg1(x) + fw_norm(x,pos[3]))/densg1(x,pars)
  
  fish.exact <- matrix(NA, 3, 3)
  fish.exact[1, 1] <- integrate(function(x) {grad1(x)^2 *densg1(x,pars)}, L1, U1, rel.tol = 1e-8)$value
  fish.exact[2, 2] <- integrate(function(x) {grad2(x)^2 *densg1(x,pars)}, L1, U1, rel.tol = 1e-8)$value
  fish.exact[3, 3] <- integrate(function(x) {grad3(x)^2 *densg1(x,pars)}, L1, U1, rel.tol = 1e-8)$value
  
  fish.exact[1, 2] <- integrate(function(x) {grad1(x)*grad2(x) *densg1(x,pars)}, L1, U1, rel.tol = 1e-8)$value
  fish.exact[1, 3] <- integrate(function(x) {grad1(x)*grad3(x) *densg1(x,pars)}, L1, U1, rel.tol = 1e-8)$value
  fish.exact[2, 3] <- integrate(function(x) {grad2(x)*grad3(x) *densg1(x,pars)}, L1, U1, rel.tol = 1e-8)$value
  fish.exact[2, 1] <- fish.exact[1, 2]
  fish.exact[3, 1] <- fish.exact[1, 3]
  fish.exact[3, 2] <- fish.exact[2, 3]
  
  fisher_sqroot_inverse = (pinv(expm::sqrtm(fish.exact)))
  norm_score_func = function(x) fisher_sqroot_inverse %*% c(grad1(x),grad2(x),grad3(x)) 
  norm_score_func1 = function(x) as.numeric(fisher_sqroot_inverse[1,] %*% c(grad1(x),grad2(x),grad3(x)))
  norm_score_func2 = function(x) as.numeric(fisher_sqroot_inverse[2,] %*% c(grad1(x),grad2(x),grad3(x)))
  norm_score_func3 = function(x) as.numeric(fisher_sqroot_inverse[3,] %*% c(grad1(x),grad2(x),grad3(x)))
  norm_score_func = Vectorize(norm_score_func)
  norm_score_func1 = Vectorize(norm_score_func1)
  norm_score_func2 = Vectorize(norm_score_func2)
  norm_score_func3 = Vectorize(norm_score_func3)
  return (list(norm_score_func=norm_score_func,
               norm_score_func1=norm_score_func1,
               norm_score_func2=norm_score_func2,
               norm_score_func3=norm_score_func3))
}
norm_score_f <- function(pars){
  grad1 <- function(x) (-bkg1(x) + dtruncnorm(x,a=L1,b=U1,mean=pos[1],sd=0.05))/densf(x,pars)
  grad2 <- function(x) (-bkg1(x) + dtruncnorm(x,a=L1,b=U1,mean=pos[2],sd=0.05))/densf(x,pars)
  grad3 <- function(x) (-bkg1(x) + dtruncnorm(x,a=L1,b=U1,mean=pos[3],sd=0.05))/densf(x,pars)
  
  fish.exact <- matrix(NA, 3, 3)
  fish.exact[1, 1] <- integrate(function(x) {grad1(x)^2 *densf(x,pars)}, L1, U1, rel.tol = 1e-8)$value
  fish.exact[2, 2] <- integrate(function(x) {grad2(x)^2 *densf(x,pars)}, L1, U1, rel.tol = 1e-8)$value
  fish.exact[3, 3] <- integrate(function(x) {grad3(x)^2 *densf(x,pars)}, L1, U1, rel.tol = 1e-8)$value
  fish.exact[1, 2] <- integrate(function(x) {grad1(x)*grad2(x) *densf(x,pars)}, L1, U1, rel.tol = 1e-8)$value
  fish.exact[1, 3] <- integrate(function(x) {grad1(x)*grad3(x) *densf(x,pars)}, L1, U1, rel.tol = 1e-8)$value
  fish.exact[2, 3] <- integrate(function(x) {grad2(x)*grad3(x) *densf(x,pars)}, L1, U1, rel.tol = 1e-8)$value
  fish.exact[2, 1] <- fish.exact[1, 2]
  fish.exact[3, 1] <- fish.exact[1, 3]
  fish.exact[3, 2] <- fish.exact[2, 3]
  
  ## Sqaure root Inverse Fisher times score functions 
  fisher_sqroot_inverse <- pinv(expm::sqrtm(fish.exact))
  norm_score_func = function(x) fisher_sqroot_inverse %*% c(grad1(x),grad2(x),grad3(x)) 
  norm_score_func1 = function(x) as.numeric(fisher_sqroot_inverse[1,] %*% c(grad1(x),grad2(x),grad3(x)))
  norm_score_func2 = function(x) as.numeric(fisher_sqroot_inverse[2,] %*% c(grad1(x),grad2(x),grad3(x)))
  norm_score_func3 = function(x) as.numeric(fisher_sqroot_inverse[3,] %*% c(grad1(x),grad2(x),grad3(x)))
  norm_score_func = Vectorize(norm_score_func)
  norm_score_func1 = Vectorize(norm_score_func1)
  norm_score_func2 = Vectorize(norm_score_func2)
  norm_score_func3 = Vectorize(norm_score_func3)
  return (list(norm_score_func=norm_score_func,
               norm_score_func1=norm_score_func1,
               norm_score_func2=norm_score_func2,
               norm_score_func3=norm_score_func3))
}

score_g1 = norm_score_g1(MLEg1)
norm_score_g11 = score_g1$norm_score_func1
norm_score_g12 = score_g1$norm_score_func2
norm_score_g13 = score_g1$norm_score_func3

score_f = norm_score_f(MLEf)
norm_score_f1 = score_f$norm_score_func1
norm_score_f2 = score_f$norm_score_func2
norm_score_f3 = score_f$norm_score_func3

#########################   K2   ######################### 
### Specify isometry and apply it to the orthornomal basis functions h and orthonormalzied score functions b
M1=6
isometry=function(x) sqrt(densf(x,MLEf)/densg1(x,MLEg1))
lh=list();
for (i in 1:(M1)){
    call = paste("lh[[",i,"]]<-function(x) isometry(x) * hj[[",i+1,"]](G_f(x,MLEf))",sep="")
    eval(parse(text=call))
}
lb=list()
for (i in 1:3){
    call = paste("lb[[",i,"]]<-function(x) isometry(x) *norm_score_f",i,"(x)",sep="")
    eval(parse(text=call))

}

### Calculate necessary inner products for K
innerl = integrate(function(x) {sqrt(densf(x,MLEf)*densg1(x,MLEg1))},lower=L1,upper=U1,rel.tol=1e-8)$value
innerlh = innerlb = c()
for (i in 1:(M1)){
  innerlh[i]=integrate(function(x) {lh[[i]](x)*densg1(x,MLEg1)},lower=L1,upper=U1,rel.tol=1e-8)$value
}
for (i in 1:3){
  innerlb[i]=integrate(function(x) {lb[[i]](x)*densg1(x,MLEg1)},lower=L1,upper=U1,rel.tol=1e-8)$value
}

### Apply K to lh and lb 
Klh=list();Klb=list()
for (i in 1:(M1)){
    call = paste("Klh[[",i,"]]<-function(x) lh[[",i,"]](x)- (1-isometry(x))/(1-innerl)*innerlh[",i,"]",sep="")
    eval(parse(text=call))
}
for (i in 1:3){
    call = paste("Klb[[",i,"]]<-function(x) lb[[",i,"]](x)- (1-isometry(x))/(1-innerl)*innerlb[",i,"]",sep="")
    eval(parse(text=call))
}

### Calculate necessary inner products for U
inner_step3_1_1 = integrate(function(x) {(norm_score_g11(x)-Klb[[1]](x))*Klb[[2]](x)*densg1(x,MLEg1)},lower=L1, upper=U1,rel.tol=1e-8)$value
inner_step3_2_1 = integrate(function(x) {norm_score_g11(x)*Klb[[1]](x)*densg1(x,MLEg1)},lower=L1 ,upper=U1,rel.tol=1e-8)$value
tildec2_g1 = function(x){Klb[[2]](x)-inner_step3_1_1/(1-inner_step3_2_1)*(norm_score_g11(x)-Klb[[1]](x))} # Same with below

inner_step3_c3 = integrate(function(x) {(norm_score_g11(x)-Klb[[1]](x))*Klb[[3]](x)*densg1(x,MLEg1)},lower=L1, upper=U1,rel.tol=1e-8)$value
Ua1c1c3 = function(x) {Klb[[3]](x) - inner_step3_c3/(1-inner_step3_2_1)*(norm_score_g11(x)-Klb[[1]](x))}
inner_a2ctile2minus = integrate(function(x) {(norm_score_g12(x)-tildec2_g1(x))*Ua1c1c3(x)*densg1(x,MLEg1)},lower=L1,upper=U1,rel.tol=1e-8)$value
inner_a2ctile2 = integrate(function(x) {norm_score_g12(x)*tildec2_g1(x)*densg1(x,MLEg1)}, lower=L1 ,upper=U1,rel.tol=1e-8)$value
tildec3_g1 = function(x) {Ua1c1c3(x) - inner_a2ctile2minus/(1-inner_a2ctile2)*(norm_score_g12(x)-tildec2_g1(x))}

### Apply U to Klh and Klb 
inner1_h = inner_U1_h = inner_UU1_h=c()
U1_h = UU1_h = UKlh = list()
inner_Ua2c2 = integrate(function(x) {norm_score_g12(x)*tildec2_g1(x)*densg1(x,MLEg1)},lower=L1,upper=U1,rel.tol=1e-8)$value
inner_Ua3c3 = integrate(function(x) {norm_score_g13(x)*tildec3_g1(x)*densg1(x,MLEg1)},lower=L1,upper=U1,rel.tol=1e-8)$value

for (i in 1:(M1)){
  inner1_h[i]=integrate(function(x) {(norm_score_g11(x)-Klb[[1]](x))*Klh[[i]](x)*densg1(x,MLEg1)}, lower=L1, upper=U1, rel.tol=1e-6)$value
  call1 = paste("U1_h[[",i,"]]<-function(x) {Klh[[",i,"]](x)-inner1_h[",i,"]/(1-inner_step3_2_1)*(norm_score_g11(x)-Klb[[1]](x))}",sep="")
  eval(parse(text=call1))
  inner_U1_h[i] = integrate(function(x) {(norm_score_g12(x)-tildec2_g1(x))*U1_h[[i]](x)*densg1(x,MLEg1)}, lower=L1, upper=U1, rel.tol=1e-6)$value
  call2 = paste("UU1_h[[",i,"]]<-function(x) {U1_h[[",i,"]](x)-inner_U1_h[",i,"]/(1-inner_Ua2c2)*(norm_score_g12(x)-tildec2_g1(x))}",sep="")
  eval(parse(text=call2))
  inner_UU1_h[i] = integrate(function(x){(norm_score_g13(x)-tildec3_g1(x))*UU1_h[[i]](x)*densg1(x,MLEg1)},lower=L1,upper=U1,rel.tol=1e-8)$value
  call3 = paste("UKlh[[",i,"]]<-function(x) UU1_h[[",i,"]](x)- inner_UU1_h[",i,"]/(1-inner_Ua3c3)*(norm_score_g13(x) - tildec3_g1(x))",sep="")
  eval(parse(text=call3))
}


###Calculate tilde h for F based on the classical basis, that is, the composite of the shifted Legendre polynomials and CDF
vj[[1]] = function(x) return (1)
for (i in 1:(M1)){
  call01 = paste("inner",i,"_1=integrate(function(x) {hj[[",i+1
                 ,"]](G_f(x,MLEf))*norm_score_f1(x)*densf(x,MLEf) }, L1, U1, rel.tol = 1e-8)$value",sep="")
  call02 = paste("inner",i,"_2=integrate(function(x) {hj[[",i+1
                 ,"]](G_f(x,MLEf))*norm_score_f2(x)*densf(x,MLEf) }, L1, U1, rel.tol = 1e-8)$value",sep="")
  call03 = paste("inner",i,"_3=integrate(function(x) {hj[[",i+1
                 ,"]](G_f(x,MLEf))*norm_score_f3(x)*densf(x,MLEf) }, L1, U1, rel.tol = 1e-8)$value",sep="")
  call04 = paste("vj[[",i+1,"]]= function(x) {hj[[",i+1,"]](G_f(x,MLEf)) - norm_score_f1(x) * inner",i,
                 "_1 - norm_score_f2(x) * inner",i,"_2 - norm_score_f3(x) * inner",i,"_3}",sep="")
  eval(parse(text=call01))
  eval(parse(text=call02))
  eval(parse(text=call03))
  eval(parse(text=call04))
  print(i)
}

###Calculate tilde h for F based on the K2 basis
Kvj1=list()
Kvj1[[1]] = function(x) return (1)
for (i in 1:M1){
  call001 = paste("inner",i,"_11=integrate(function(x) {hj[[",i+1
                  ,"]](G_f(x,MLEf))*norm_score_f1(x)*densf(x,MLEf) }, L1, U1, rel.tol = 1e-8)$value",sep="")
  call002 = paste("inner",i,"_12=integrate(function(x) {hj[[",i+1
                  ,"]](G_f(x,MLEf))*norm_score_f2(x)*densf(x,MLEf) }, L1, U1, rel.tol = 1e-8)$value",sep="")
  call003 = paste("inner",i,"_13=integrate(function(x) {hj[[",i+1
                  ,"]](G_f(x,MLEf))*norm_score_f3(x)*densf(x,MLEf) }, L1, U1, rel.tol = 1e-8)$value",sep="")
  call004 = paste("Kvj1[[",i+1,"]] = function(x) {UKlh[[",i,"]](x) - norm_score_g11(x) * inner",i,
                  "_11 - norm_score_g12(x) * inner",i,"_12 - norm_score_g13(x) * inner",i,"_13 }",sep="")
  eval(parse(text=call001))
  eval(parse(text=call002))
  eval(parse(text=call003))
  eval(parse(text=call004))
  print(i)
}


################### Simulations ###########################
fw_dist=Vectorize(fw_dist)
G_g1<-function(x,pars){
  res = (1-sum(pars[1:3]))*(2.5*x-4.125)+pars[1]*fw_dist(x,pos[1])+pars[2]*fw_dist(x,pos[2])+pars[3]*fw_dist(x,pos[3])
  return (res)
}
###### From F ######
set.seed(1)
start_time1 <- Sys.time()
B <- 100000
res_tn_order_boot_f <- res_tn_boot_f  <- numeric(B); vj_boot_q <- matrix(NA,nrow=B,ncol=M1)
for (k in 1:B){
  num = rmultinom(n=1,size=n,prob=c(1-sum(MLEf),MLEf))
  xx_f = c(runif(num[1],L1,U1),rtruncnorm(num[2], a=L1, b=U1, mean = pos[1], sd = 0.05),
           rtruncnorm(num[3], a=L1, b=U1, mean = pos[2], sd = 0.05), 
           rtruncnorm(num[4], a=L1, b=U1, mean = pos[3], sd = 0.05))
  for (j in 1:M1){
    vj_boot_q[k,j] <- sqrt(n)*mean(vj[[j+1]](xx_f))}
  res_tn_boot_f[k] <- max(cumsum((vj_boot_q[k,])^2)/1:M1)
  res_tn_order_boot_f[k] <- max(cumsum(sort((vj_boot_q[k,])^2,decreasing=TRUE))/1:M1)
  print(k)
}
end_time1 <- Sys.time()
end_time1-start_time1


#################################### p-values #################################### 
vn_bs11 = c()
for (j in 1:M1){
  vn_bs11[j] <- sqrt(n)*mean(Kvj1[[j+1]](bs11))
}
res_tn_bs11 <- max(cumsum((vn_bs11)^2)/1:M1)
res_tn_order_bs11 <- max(cumsum(sort((vn_bs11)^2,decreasing=TRUE))/1:M1)

print((sum(res_tn_boot_f>=res_tn_bs11)+1)/(B+1))
print((sum(res_tn_order_boot_f>=res_tn_order_bs11)+1)/(B+1))


#################################### QQ-plots #################################### 
### This generates the Figure 4 in Section 5. Warning: this will be slow

sampler_c = function(n, cc){
  tt = rlst(2*n,df=4,sigma=0.025)
  res = (tt[tt>=L1+0.015-cc & tt<=U1-0.015-cc])[1:n]
  return(rtruncnorm(n,a=cc-0.015,b=cc+0.015,cc,sigma) + res)
}
sampler_x1 = function(n){
  indx = rmultinom(1,size=n,prob=c(1-sum(MLEg1),MLEg1))
  return (c(runif(indx[1],L1,U1), sampler_c(indx[2],pos[1]), sampler_c(indx[3],pos[2]),sampler_c(indx[4],pos[3])))
}

set.seed(1)
start_time2 <- Sys.time()
res_tn_order_boot_K2_g1 <- res_tn_boot_K2_g1 <- numeric(B);vn_boot_K2_g1 <- matrix(NA,nrow=B,ncol=M1)
for (k in 1:B){
  xx_k2 = sampler_x1(n)
  for (j in 1:M1){
    vn_boot_K2_g1[k,j] <- sqrt(n)*mean(Kvj1[[j+1]](xx_k2))
  }
  res_tn_boot_K2_g1[k] <- max(cumsum((vn_boot_K2_g1[k,])^2)/1:M1)
  res_tn_order_boot_K2_g1[k] <- max(cumsum(sort((vn_boot_K2_g1[k,])^2,decreasing=TRUE))/1:M1)
  print(k)
}
end_time2 <- Sys.time()
end_time2-start_time2


#################################### QQplots #################################### 
setEPS()
postscript("Real_plot11.eps",width = 9, height = 6)
par(mar = c(4, 5.5, 2, 2))
qqplot(res_tn_boot_K2_g1, res_tn_boot_f, ylab=expression("Quantiles of "*hat(T)[n]*" from "*G[beta]),
       xlab=expression("Quantiles of "*hat(T)[n]^K*" from "*F[gamma]*""),
       lwd = 2, cex.axis = 1.6, cex.lab = 1.8)
abline(0,1)
dev.off()


## Order selection statistics, ordered
setEPS()
postscript("Real_plot22.eps",width = 9, height = 6)
par(mar = c(4, 5.5, 2, 2))
qqplot(res_tn_order_boot_K2_g1, res_tn_order_boot_f, ylab=expression("Quantiles of "*tilde(T)[n]*" from "*G[beta]),
       xlab=expression("Quantiles of "*tilde(T)[n]^K*" from "*F[gamma]),
       lwd = 2, cex.axis = 1.6, cex.lab = 1.8)
abline(0,1)
dev.off()

