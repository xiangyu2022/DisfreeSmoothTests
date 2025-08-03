source("hjs.R")
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
bs1<-read.table("Downloads/16688_src.txt")
bs2<-read.table("Downloads/18710_src.txt")
bg1<-read.table("Downloads/16688_bkg.txt")
bg2<-read.table("Downloads/18710_bkg.txt")
bs_full<-c(bs1[,1],bs2[,1])
bs_abs<-abs(bs_full)
bg_full<-c(bg1[,1],bg2[,1])
bg_abs<-abs(bg_full)
L1=1.65;U1=2.05;
bs11<-bs_abs[bs_abs>=L1&bs_abs<=U1]
n=length(bs11)
pos = c(1.78499, 1.85247, 1.94365)
M1 = 6
bkg1<-function(x) ifelse(x>=L1&x<=U1, 1/(U1-L1), 0)
set.seed(1)
start_time <- Sys.time()
sigma=0.0025
new_fe = function(x,c) ifelse(x>=L1+0.015-c & x<=U1-0.015-c,dlst(x, df=4, sigma=0.025)/(plst(U1-0.015-c,df=4,sigma=0.025)-plst(L1+0.015-c, df=4, sigma=0.025)),0)
line_int<-function(w0,x,c){dtruncnorm(w0,a=c-0.015,b=c+0.015,c,sigma)*new_fe(x-w0,c)}
fw<-function(x,c){ifelse(x>=L1&x<=U1,integrate(line_int,lower=c-0.015,upper=c+0.015,x=x,c=c,rel.tol=1e-8)$value,0)}
fw_norm<-Vectorize(fw)

new_fe_dist = function(x,c) (plst(x, df=4, sigma=0.025)-plst(L1+0.015-c,df=4, sigma=0.025))/(plst(U1-0.015-c, df=4, sigma=0.025)-plst(L1+0.015-c, df=4, sigma=0.025))
dist_int<-function(w0,x,c){dtruncnorm(w0,a=c-0.015,b=c+0.015,c,sigma)*new_fe_dist(x-w0,c)}
fw_dist<-function(x,c){integrate(dist_int,lower=c-0.015,upper=c+0.015,x=x,c=c,rel.tol = 1e-8)$value}
fw_dist<-Vectorize(fw_dist)

densf1 <- function(x,pars){
  (1-sum(pars[1:3]))*bkg1(x)+pars[1]*fw_norm(x,pos[1])+pars[2]*fw_norm(x,pos[2])+pars[3]*fw_norm(x,pos[3])
}
gradf1 <- function(pars){
  grad1 <- function(x) (-bkg1(x) + fw_norm(x,pos[1]))/densf1(x,pars)
  grad2 <- function(x) (-bkg1(x) + fw_norm(x,pos[2]))/densf1(x,pars)
  grad3 <- function(x) (-bkg1(x) + fw_norm(x,pos[3]))/densf1(x,pars)
  return(c(-sum(grad1(bs11)), -sum(grad2(bs11)), -sum(grad3(bs11))))
}


start_time1 <- Sys.time()
ll_bs_new1<-function(pars){-sum(log(densf1(bs11,pars)))}
MLEf1 = constrOptim(theta=c(0.05,0.05,0.2), grad=gradf1, method ="BFGS",
                    f=ll_bs_new1, ui= rbind(diag(3),-diag(3),c(1,1,1),c(-1,-1,-1)),
                    ci=c(rep(0,3),rep(-1,3),0,-1))$par
end_time1 <- Sys.time()
end_time1 - start_time1

f_signal_nor1 <- function(x) ifelse(x>=1.65&x<=2.05, densf1(x,MLEf1),0)
integrate(f_signal_nor1,1.65,2.05)

start_time2 <- Sys.time()
densq <- function(x,pars){
  (1-sum(pars[1:3]))*bkg1(x)+pars[1]*dtruncnorm(x,a=L1,b=U1,mean=pos[1],sd=0.05)+
    pars[2]*dtruncnorm(x,a=L1,b=U1,mean=pos[2],sd=0.05)+pars[3]*dtruncnorm(x,a=L1,b=U1,mean=pos[3],sd=0.05)
}
likeq <- function(pars) -sum(log(densq(bs11,pars)))
MLEq = constrOptim(theta=c(0.05,0.05,0.2), grad=NULL,
                   f=likeq, ui= rbind(diag(3),-diag(3),c(1,1,1),c(-1,-1,-1)),
                   ci=c(rep(0,3),rep(-1,3),0,-1))$par
end_time2 <- Sys.time()
end_time2-start_time2
q_signal <- function(x) ifelse(x>=1.65&x<=2.05, densq(x,MLEq),0)
q_signal_nor <- function(x) q_signal(x)/integrate(q_signal,L1,U1)$value
G_q<-function(x,pars){
  res = (1-sum(pars[1:3]))*(2.5*x-4.125)+pars[1]*ptruncnorm(x,L1,U1,pos[1],sd=0.05)+pars[2]*ptruncnorm(x,L1,U1,pos[2],sd=0.05)+pars[3]*ptruncnorm(x,L1,U1,pos[3],sd=0.05)
  return (res)
}
## Plot their densities
hist(bs11,prob=T,xlim=c(1.65,2.05),ylim=c(0,5),breaks =25)
curve(f_signal_nor1,add=T)
curve(q_signal_nor,add=T,col="red")



#########################   K2   ######################### 
######################### Step 1 ######################### 
isometry1=function(x) sqrt(densq(x,MLEq)/densf1(x,MLEf1))
norm_score_q <- function(pars){
  grad1 <- function(x) (-bkg1(x) + dtruncnorm(x,a=L1,b=U1,mean=pos[1],sd=0.05))/densq(x,pars)
  grad2 <- function(x) (-bkg1(x) + dtruncnorm(x,a=L1,b=U1,mean=pos[2],sd=0.05))/densq(x,pars)
  grad3 <- function(x) (-bkg1(x) + dtruncnorm(x,a=L1,b=U1,mean=pos[3],sd=0.05))/densq(x,pars)
  
  fish.exact <- matrix(NA, 3, 3)
  fish.exact[1, 1] <- integrate(function(x) {grad1(x)^2 *densq(x,pars)}, L1, U1, rel.tol = 1e-8)$value
  fish.exact[2, 2] <- integrate(function(x) {grad2(x)^2 *densq(x,pars)}, L1, U1, rel.tol = 1e-8)$value
  fish.exact[3, 3] <- integrate(function(x) {grad3(x)^2 *densq(x,pars)}, L1, U1, rel.tol = 1e-8)$value
  fish.exact[1, 2] <- integrate(function(x) {grad1(x)*grad2(x) *densq(x,pars)}, L1, U1, rel.tol = 1e-8)$value
  fish.exact[1, 3] <- integrate(function(x) {grad1(x)*grad3(x) *densq(x,pars)}, L1, U1, rel.tol = 1e-8)$value
  fish.exact[2, 3] <- integrate(function(x) {grad2(x)*grad3(x) *densq(x,pars)}, L1, U1, rel.tol = 1e-8)$value
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
norm_score_f1 <- function(pars){
  grad1 <- function(x) (-bkg1(x) + fw_norm(x,pos[1]))/densf1(x,pars)
  grad2 <- function(x) (-bkg1(x) + fw_norm(x,pos[2]))/densf1(x,pars)
  grad3 <- function(x) (-bkg1(x) + fw_norm(x,pos[3]))/densf1(x,pars)
  
  fish.exact <- matrix(NA, 3, 3)
  fish.exact[1, 1] <- integrate(function(x) {grad1(x)^2 *densf1(x,pars)}, L1, U1, rel.tol = 1e-8)$value
  fish.exact[2, 2] <- integrate(function(x) {grad2(x)^2 *densf1(x,pars)}, L1, U1, rel.tol = 1e-8)$value
  fish.exact[3, 3] <- integrate(function(x) {grad3(x)^2 *densf1(x,pars)}, L1, U1, rel.tol = 1e-8)$value
  
  fish.exact[1, 2] <- integrate(function(x) {grad1(x)*grad2(x) *densf1(x,pars)}, L1, U1, rel.tol = 1e-8)$value
  fish.exact[1, 3] <- integrate(function(x) {grad1(x)*grad3(x) *densf1(x,pars)}, L1, U1, rel.tol = 1e-8)$value
  fish.exact[2, 3] <- integrate(function(x) {grad2(x)*grad3(x) *densf1(x,pars)}, L1, U1, rel.tol = 1e-8)$value
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

score_q = norm_score_q(MLEq)
norm_score_q1 = score_q$norm_score_func1
norm_score_q2 = score_q$norm_score_func2
norm_score_q3 = score_q$norm_score_func3

score_f1 = norm_score_f1(MLEf1)
norm_score_f11 = score_f1$norm_score_func1
norm_score_f12 = score_f1$norm_score_func2
norm_score_f13 = score_f1$norm_score_func3

l1h=list()
for (i in 1:(M1)){
  for (j in 1){
    call = paste("l",j,"h[[",i,"]]<-function(x) isometry",j,"(x) * hj[[",i+1,"]](G_q(x,MLEq))",sep="")
    eval(parse(text=call))
  }
}
l1b=list()
for (i in 1:3){
  for (j in 1){
    call = paste("l",j,"b[[",i,"]]<-function(x) isometry",j,"(x) *norm_score_q",i,"(x)",sep="")
    eval(parse(text=call))
  }
}
for (i in 1:M1){
  print(integrate(function(x)l1h[[i]](x)*l1h[[i]](x)*densf1(x,MLEf1),lower=L1,upper=U1))
  print(integrate(function(x)hj[[i]](G_q(x,MLEq))*hj[[i]](G_q(x,MLEq))*densq(x,MLEq),lower=L1,upper=U1))
}

######################### Step 2 ######################### 
innerl1 = integrate(function(x) {sqrt(densq(x,MLEq)*densf1(x,MLEf1))},lower=L1,upper=U1,rel.tol=1e-8)$value
inner1lh = inner1lb = c()
for (i in 1:(M1)){
  inner1lh[i]=integrate(function(x) {l1h[[i]](x)*densf1(x,MLEf1)},lower=L1,upper=U1,rel.tol=1e-8)$value
}
for (i in 1:3){
  inner1lb[i]=integrate(function(x) {l1b[[i]](x)*densf1(x,MLEf1)},lower=L1,upper=U1,rel.tol=1e-8)$value
}

Kl1h=list();Kl1b=list()
for (i in 1:(M1)){
  for (j in 1){
    call = paste("Kl",j,"h[[",i,"]]<-function(x) l",j,"h[[",i,"]](x)- (1-isometry",j,"(x))/(1-innerl",j,")*inner",j,"lh[",i,"]",sep="")
    eval(parse(text=call))
  }
}
for (i in 1:3){
  for (j in 1){
    call = paste("Kl",j,"b[[",i,"]]<-function(x) l",j,"b[[",i,"]](x)- (1-isometry",j,"(x))/(1-innerl",j,")*inner",j,"lb[",i,"]",sep="")
    eval(parse(text=call))
  }
}


##################################K2-STEP3##################################
### Tilde c2: U_a1c1 c2
inner_step3_1_1 = integrate(function(x) {(norm_score_f11(x)-Kl1b[[1]](x))*Kl1b[[2]](x)*densf1(x,MLEf1)},lower=L1, upper=U1,rel.tol=1e-8)$value
inner_step3_2_1 = integrate(function(x) {norm_score_f11(x)*Kl1b[[1]](x)*densf1(x,MLEf1)},lower=L1 ,upper=U1,rel.tol=1e-8)$value
tildec2_f1 = function(x){Kl1b[[2]](x)-inner_step3_1_1/(1-inner_step3_2_1)*(norm_score_f11(x)-Kl1b[[1]](x))} # Same with below
integrate(function(x)tildec2_f1(x)*norm_score_f11(x)*densf1(x,MLEf1),L1,U1)

### Tilde c3: U_a2tildec2 U_a1c1 c3 
inner_step3_c3 = integrate(function(x) {(norm_score_f11(x)-Kl1b[[1]](x))*Kl1b[[3]](x)*densf1(x,MLEf1)},lower=L1, upper=U1,rel.tol=1e-8)$value
Ua1c1c3 = function(x) {Kl1b[[3]](x) - inner_step3_c3/(1-inner_step3_2_1)*(norm_score_f11(x)-Kl1b[[1]](x))}
inner_a2ctile2minus = integrate(function(x) {(norm_score_f12(x)-tildec2_f1(x))*Ua1c1c3(x)*densf1(x,MLEf1)},lower=L1,upper=U1,rel.tol=1e-8)$value
inner_a2ctile2 = integrate(function(x) {norm_score_f12(x)*tildec2_f1(x)*densf1(x,MLEf1)}, lower=L1 ,upper=U1,rel.tol=1e-8)$value
tildec3_f1 = function(x) {Ua1c1c3(x) - inner_a2ctile2minus/(1-inner_a2ctile2)*(norm_score_f12(x)-tildec2_f1(x))}

##################################K2-STEP4##################################

### Ua3tildec3 Ua2tildec2 Ua1c1 Klh
inner1_h = inner_U1_h = inner_UU1_h=c()
U1_h = UU1_h = UKl1h = list()
inner_Ua2c2 = integrate(function(x) {norm_score_f12(x)*tildec2_f1(x)*densf1(x,MLEf1)},lower=L1,upper=U1,rel.tol=1e-8)$value
inner_Ua3c3 = integrate(function(x) {norm_score_f13(x)*tildec3_f1(x)*densf1(x,MLEf1)},lower=L1,upper=U1,rel.tol=1e-8)$value

for (i in 1:(M1)){
  inner1_h[i]=integrate(function(x) {(norm_score_f11(x)-Kl1b[[1]](x))*Kl1h[[i]](x)*densf1(x,MLEf1)}, lower=L1, upper=U1, rel.tol=1e-6)$value
  call1 = paste("U1_h[[",i,"]]<-function(x) {Kl1h[[",i,"]](x)-inner1_h[",i,"]/(1-inner_step3_2_1)*(norm_score_f11(x)-Kl1b[[1]](x))}",sep="")
  eval(parse(text=call1))
  inner_U1_h[i] = integrate(function(x) {(norm_score_f12(x)-tildec2_f1(x))*U1_h[[i]](x)*densf1(x,MLEf1)}, lower=L1, upper=U1, rel.tol=1e-6)$value
  call2 = paste("UU1_h[[",i,"]]<-function(x) {U1_h[[",i,"]](x)-inner_U1_h[",i,"]/(1-inner_Ua2c2)*(norm_score_f12(x)-tildec2_f1(x))}",sep="")
  eval(parse(text=call2))
  inner_UU1_h[i] = integrate(function(x){(norm_score_f13(x)-tildec3_f1(x))*UU1_h[[i]](x)*densf1(x,MLEf1)},lower=L1,upper=U1,rel.tol=1e-8)$value
  call3 = paste("UKl1h[[",i,"]]<-function(x) UU1_h[[",i,"]](x)- inner_UU1_h[",i,"]/(1-inner_Ua3c3)*(norm_score_f13(x) - tildec3_f1(x))",sep="")
  eval(parse(text=call3))
}


##################################Calculating tilde h or tilde psi##################################
vj1=list();vj=list()
vj[[1]] = function(x) return (1)
for (i in 1:(M1)){
  call01 = paste("inner",i,"_1=integrate(function(x) {hj[[",i+1
                 ,"]](G_q(x,MLEq))*norm_score_q1(x)*densq(x,MLEq) }, L1, U1, rel.tol = 1e-8)$value",sep="")
  call02 = paste("inner",i,"_2=integrate(function(x) {hj[[",i+1
                 ,"]](G_q(x,MLEq))*norm_score_q2(x)*densq(x,MLEq) }, L1, U1, rel.tol = 1e-8)$value",sep="")
  call03 = paste("inner",i,"_3=integrate(function(x) {hj[[",i+1
                 ,"]](G_q(x,MLEq))*norm_score_q3(x)*densq(x,MLEq) }, L1, U1, rel.tol = 1e-8)$value",sep="")
  call04 = paste("vj[[",i+1,"]]= function(x) {hj[[",i+1,"]](G_q(x,MLEq)) - norm_score_q1(x) * inner",i,
                 "_1 - norm_score_q2(x) * inner",i,"_2 - norm_score_q3(x) * inner",i,"_3}",sep="")
  eval(parse(text=call01))
  eval(parse(text=call02))
  eval(parse(text=call03))
  eval(parse(text=call04))
  print(i)
}

Kvj1=list()
Kvj1[[1]] = function(x) return (1)
for (i in 1:M1){
  call001 = paste("inner",i,"_11=integrate(function(x) {hj[[",i+1
                  ,"]](G_q(x,MLEq))*norm_score_q1(x)*densq(x,MLEq) }, L1, U1, rel.tol = 1e-8)$value",sep="")
  call002 = paste("inner",i,"_12=integrate(function(x) {hj[[",i+1
                  ,"]](G_q(x,MLEq))*norm_score_q2(x)*densq(x,MLEq) }, L1, U1, rel.tol = 1e-8)$value",sep="")
  call003 = paste("inner",i,"_13=integrate(function(x) {hj[[",i+1
                  ,"]](G_q(x,MLEq))*norm_score_q3(x)*densq(x,MLEq) }, L1, U1, rel.tol = 1e-8)$value",sep="")
  call004 = paste("Kvj1[[",i+1,"]] = function(x) {UKl1h[[",i,"]](x) - norm_score_f11(x) * inner",i,
                  "_11 - norm_score_f12(x) * inner",i,"_12 - norm_score_f13(x) * inner",i,"_13 }",sep="")
  eval(parse(text=call001))
  eval(parse(text=call002))
  eval(parse(text=call003))
  eval(parse(text=call004))
  print(i)
}

varcov_q=varcov_f1=matrix(nrow=M1+1,ncol=M1+1)
for (m in 1:(M1+1)){
  for (s in 1:m)
    varcov_q[m,s]= integrate(function(x) {vj[[m]](x)*vj[[s]](x)*densq(x,MLEq)}, L1, U1, rel.tol=1e-8)$value
}
for (m in 1:(M1+1)){
  for (s in 1:m){
    varcov_f1[m,s]= integrate(function(x) {Kvj1[[m]](x)*Kvj1[[s]](x)*densf1(x,MLEf1)}, L1, U1, rel.tol=1e-8)$value
  }
}

################### Simulations ###########################
fw_dist=Vectorize(fw_dist)
G_f1<-function(x,pars){
  res = (1-sum(pars[1:3]))*(2.5*x-4.125)+pars[1]*fw_dist(x,pos[1])+pars[2]*fw_dist(x,pos[2])+pars[3]*fw_dist(x,pos[3])
  return (res)
}
###### From Q ######
set.seed(1)
start_time1 <- Sys.time()
B <- 1000
res_tn_order_boot_q <- res_tn_boot_q <-res_gst_boot_q <- res_ks_boot_q <- numeric(B)
vj_boot_q <- matrix(NA,nrow=B,ncol=M1)
for (k in 1:B){
  num = rmultinom(n=1,size=n,prob=c(1-sum(MLEq),MLEq))
  xx_q = c(runif(num[1],L1,U1),rtruncnorm(num[2], a=L1, b=U1, mean = pos[1], sd = 0.05),
           rtruncnorm(num[3], a=L1, b=U1, mean = pos[2], sd = 0.05), 
           rtruncnorm(num[4], a=L1, b=U1, mean = pos[3], sd = 0.05))
  for (j in 1:M1){
    vj_boot_q[k,j] <- sqrt(n)*mean(vj[[j+1]](xx_q))}
  res_tn_boot_q[k] <- max(cumsum((vj_boot_q[k,])^2)/1:M1)
  res_tn_order_boot_q[k] <- max(cumsum(sort((vj_boot_q[k,])^2,decreasing=TRUE))/1:M1)
  print(k)
}
end_time1 <- Sys.time()
end_time1-start_time1


###### From F1 ######
sampler_c = function(n, cc){
  tt = rlst(2*n,df=4,sigma=0.025)
  res = (tt[tt>=L1+0.015-cc & tt<=U1-0.015-cc])[1:n]
  return(rtruncnorm(n,a=cc-0.015,b=cc+0.015,cc,sigma) + res)
}
sampler_x1 = function(n){
  indx = rmultinom(1,size=n,prob=c(1-sum(MLEf1),MLEf1))
  return (c(runif(indx[1],L1,U1), sampler_c(indx[2],pos[1]), sampler_c(indx[3],pos[2]),sampler_c(indx[4],pos[3])))
}
gradf1_func <- function(xx,pars){
  grad1 <- function(x) (-bkg1(x) + fw_norm(x,pos[1]))/densf1(x,pars)
  grad2 <- function(x) (-bkg1(x) + fw_norm(x,pos[2]))/densf1(x,pars)
  grad3 <- function(x) (-bkg1(x) + fw_norm(x,pos[3]))/densf1(x,pars)
  return(c(-sum(grad1(xx)), -sum(grad2(xx)), -sum(grad3(xx))))
}
set.seed(1)
start_time2 <- Sys.time()
res_ks_boot1 <- res_gst_boot1<- res_tn_boot1 <- res_tn_order_boot1 <- numeric(B)
S1_boot <- matrix(NA,n,M1)
vj_boot<- matrix(NA,nrow=B,ncol=M1)
for (k in 1:B){
  xx_boot = sampler_x1(n)
  for (j in 1:M1){
    vj_boot[k,j] <- sqrt(n)*mean(hj[[j+1]](G_f1(xx_boot,MLEf1)))
  }
  res_tn_boot1[k] <- max(cumsum((vj_boot[k,])^2)/1:M1)
  res_tn_order_boot1[k] <- max(cumsum(sort((vj_boot[k,])^2,decreasing=TRUE))/1:M1)
  print(k)
}
end_time2 <- Sys.time()
end_time2-start_time2

###### From K2F1 ######
set.seed(1)
start_time3 <- Sys.time()
res_tn_order_boot_K2_f1 <- res_tn_boot_K2_f1 <-res_gst_boot_K2_f1<-res_ks_boot_K2_f1 <- numeric(B)
vn_boot_K2_f1 <- matrix(NA,nrow=B,ncol=M1)
for (k in 1:B){
  xx_k2 = sampler_x1(n)
  for (j in 1:M1){
    vn_boot_K2_f1[k,j] <- sqrt(n)*mean(Kvj1[[j+1]](xx_k2))
  }
  res_ks_boot_K2_f1[k] <- max(abs(vn_boot_K2_f1[k,]))
  res_tn_boot_K2_f1[k] <- max(cumsum((vn_boot_K2_f1[k,])^2)/1:M1)
  res_tn_order_boot_K2_f1[k] <- max(cumsum(sort((vn_boot_K2_f1[k,])^2,decreasing=TRUE))/1:M1)
  print(k)
}
end_time3 <- Sys.time()
end_time3-start_time3


#################################### QQplots #################################### 
setEPS()
postscript("MNRAS1.eps",width = 9, height = 6)
par(mar = c(4, 5.5, 2, 2))
qqplot(res_tn_boot_f1,res_tn_boot_q, ylab=expression("Quantiles of "*hat(T)[n]*" from "*G[beta]),
       xlab=expression("Quantiles of "*hat(T)[n]*" from "*F[gamma]*""),
       lwd = 2, cex.axis = 1.6, cex.lab = 1.8)
abline(0,1)
dev.off()

setEPS()
postscript("Real_plot11.eps",width = 9, height = 6)
par(mar = c(4, 5.5, 2, 2))
qqplot(res_tn_boot_K2_f1, res_tn_boot_q, ylab=expression("Quantiles of "*hat(T)[n]*" from "*G[beta]),
       xlab=expression("Quantiles of "*hat(T)[n]*" from "*F[gamma]*""),
       lwd = 2, cex.axis = 1.6, cex.lab = 1.8)
abline(0,1)
dev.off()


## Order selection statistics, ordered
setEPS()
postscript("MNRAS3.eps",width = 9, height = 6)
par(mar = c(4, 5.5, 2, 2))
qqplot(res_tn_boot_f1, res_tn_order_boot_q, ylab=expression("Quantiles of "*tilde(T)[n]*" from "*G[beta]),
       xlab=expression("Quantiles of "*tilde(T)[n]*" from "*F[gamma]*" without K-2"),
       lwd = 2, cex.axis = 1.6, cex.lab = 1.8)
abline(0,1)
dev.off()

setEPS()
postscript("Real_plot22.eps",width = 9, height = 6)
par(mar = c(4, 5.5, 2, 2))
qqplot(res_tn_order_boot_K2_f1, res_tn_order_boot_q, ylab=expression("Quantiles of "*tilde(T)[n]*" from "*G[beta]),
       xlab=expression("Quantiles of "*tilde(T)[n]*" from "*F[gamma]),
       lwd = 2, cex.axis = 1.6, cex.lab = 1.8)
abline(0,1)
dev.off()




#################################### p-values #################################### 
vn_bs11 = c()
for (j in 1:M1){
  vn_bs11[j] <- sqrt(n)*mean(Kvj1[[j+1]](bs11))
}
res_tn_bs11 <- max(cumsum((vn_bs11)^2)/1:M1)
res_tn_order_bs11 <- max(cumsum(sort((vn_bs11)^2,decreasing=TRUE))/1:M1)

print((sum(res_tn_boot_q>=res_tn_bs11)+1)/(B+1))
print((sum(res_tn_order_boot_q>=res_tn_order_bs11)+1)/(B+1))


print((sum(res_tn_boot_K2_f1>=res_tn_bs11)+1)/(length(res_tn_boot_K2_f1)+1))
print((sum(res_tn_order_boot_K2_f1>=res_tn_order_bs11)+1)/(length(res_tn_order_boot_K2_f1)+1))

