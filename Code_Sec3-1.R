source("hjs.R")
library(nloptr)
library(LaplacesDemon)
library(expm)
library(pracma)
library(truncnorm)
library(LPBkg)
library(orthopolynom)

### Generate from the asymmetric Laplace distribution with location 0, scale 2, and kappa 0.1
n=100;location=0;scale=2;kappa=0.1
set.seed(1)
xx=ralaplace(n, location, scale, kappa=0.1)
alpha_func=function(x,theta) mean(sapply(x-theta, function(y) max(0,y)))
beta_func=function(x,theta) mean(sapply(x-theta, function(y) -min(0,y)))
func = function(k){
  1-2*k^2/(1+k^2)+sqrt(2)/scale*(beta_func(xx,location)/k-alpha_func(xx,location)*k)
}
### Obtain the MLE for the kappa 
MLE_ALP=uniroot(func,interval=c(.Machine$double.xmin^0.1,.Machine$double.xmax^0.1),tol=1e-10)$root

### Obtain the score functions to have 
dens <- function(x,pars) dalaplace(x, location, scale, kappa=MLE_ALP)
G_obs <- function(x) {palaplace(x, location, scale, kappa=MLE_ALP)}

start_time1 <- Sys.time()
norm_score <- function(pars){
  grad1 <- function(x) 1/pars - 2*pars/(1+pars^2) - ifelse(x>=location, sqrt(2)/scale*(x-location), sqrt(2)/(pars^2*scale)*(x-location))
  fish.exact <- matrix(NA, 1, 1)
  fish.exact[1, 1] <- 1/pars^2 + 4/(1+pars^2)^2
  
  ## Sqaure root Inverse Fisher times score functions 
  fisher_sqroot_inverse <- pinv(expm::sqrtm(fish.exact))
  norm_score_func = function(x) fisher_sqroot_inverse %*% c(grad1(x)) 
  norm_score_func1 = function(x) as.numeric(fisher_sqroot_inverse[1,] %*% c(grad1(x)))
  norm_score_func = Vectorize(norm_score_func)
  norm_score_func1 = Vectorize(norm_score_func1)
  return (list(norm_score_func=norm_score_func,
               norm_score_func1=norm_score_func1))
}

norm_score_function = norm_score(MLE_ALP)
norm_score_func1 = norm_score_function$norm_score_func1

## Calculating the tilde h 
M1=10
vj1=list();vj=list()
vj[[1]] = function(x) return (1)
for (i in 1:M1+1){
  call1 = paste("inner",i,"_1=integrate(function(x) {hj[[",i
                ,"]](G_obs(x))*norm_score_func1(x)*dens(x,MLEALP) }, -Inf, Inf)$value"
                ,sep="")
  call2 = paste("vj[[",i,"]]= function(x) {hj[[",i,"]](G_obs(x)) - norm_score_func1(x) * inner",i,
                "_1}",sep="")
  eval(parse(text=call1))
  print(paste(i,1))
  eval(parse(text=call2))
  print(paste(i,2))
}

B=100000
## Projected bootstrap
set.seed(1)
vn_boot <- numeric(M1)
res_tn_proj <- res_tn_order_proj <- numeric(B);
for (i in 1:B){
  xx_obs= ralaplace(n, location, scale, kappa=MLE_ALP)
  for (j in 1:M1){
    vn_boot[j] <- sqrt(n)*mean(vj[[j+1]](xx_obs))
  }
  res_tn_proj[i] <-  max(cumsum((vn_boot)^2)/1:M1)
  res_tn_order_proj[i] <- max(cumsum(sort((vn_boot)^2, decreasing=TRUE)/1:M1))
  print(i)
}
end_time1 <- Sys.time()
end_time1 - start_time1

## Classical parametric bootstrap  

start_time2 <- Sys.time()
set.seed(1)
res_tn_noproj  <- res_tn_order_noproj  <- numeric(B);vn_boot_noproj <- numeric(M1)
for (i in 1:B){
  xx_boot = ralaplace(n, location, scale, kappa=MLE_ALP)
  func_boot = function(k){1-2*k^2/(1+k^2)+sqrt(2)/scale*(beta_func(xx_boot,location)/k-alpha_func(xx_boot,location)*k)}
  MLE_ALP_boot = uniroot(func_boot,interval=c(.Machine$double.xmin^0.1,.Machine$double.xmax^0.1),tol=1e-10)$root
  G_boot <- function(x) palaplace(x, location, scale, kappa=MLE_ALP_boot)
  for (j in 1:M1){
    vn_boot_noproj[j] <- sqrt(n)*mean(hj[[j+1]](G_boot(xx_boot)))
    }
  res_tn_noproj[i] <-  max(cumsum((vn_boot_noproj)^2)/1:M1)
  res_tn_order_noproj[i] <- max(cumsum(sort((vn_boot_noproj)^2, decreasing=TRUE)/1:M1))
  print(i)
}
end_time2 <- Sys.time()
end_time2 - start_time2

## Monte Carlo
set.seed(1)
start_time3 <- Sys.time()
res_tn_noproj_MC  <- res_tn_order_noproj_MC  <- numeric(B);vn_boot_noproj_MC <- numeric(M1)

for (i in 1:B){
  xx_MC= ralaplace(n, location, scale, kappa=0.1)
  func_MC = function(k){
    1-2*k^2/(1+k^2)+sqrt(2)/scale*(beta_func(xx_MC,location)/k-alpha_func(xx_MC,location)*k)}
  MLEALP_MC = uniroot(func_MC,interval=c(.Machine$double.xmin^0.1,.Machine$double.xmax^0.1),tol=1e-10)$root
  G_MC <- function(x) palaplace(x,location, scale, kappa=MLEALP_MC)
  for (j in 1:M1){
    vn_boot_noproj_MC[j] <- sqrt(n)*mean(hj[[j+1]](G_MC(xx_MC)))
  }
  res_tn_noproj_MC[i] <- max(cumsum((vn_boot_noproj_MC)^2)/1:M1)
  res_tn_order_noproj_MC[i] <- max(cumsum(sort((vn_boot_noproj_MC)^2, decreasing=TRUE)/1:M1))
  
  print(i)
}
end_time3 <- Sys.time()
end_time3 - start_time3

# save.image("Sec3-1.RData")
# load(("Sec3-1.RData"))

cc=seq(0,4,length=1000)
prob_tn_proj <- prob_tn_noproj <- prob_tn_noproj_MC <- numeric(1000)
for (k in 1:1000){
  prob_tn_proj[k] <- sum(res_tn_proj<cc[k])/length(res_tn_proj)
  prob_tn_noproj[k] <- sum(res_tn_noproj<cc[k])/length(res_tn_noproj)
  prob_tn_noproj_MC[k] <- sum(res_tn_noproj_MC<cc[k])/length(res_tn_noproj_MC)
}

postscript("tnstat.eps", width = 8, height = 6, horizontal = FALSE, onefile = FALSE)
par(mar = c(5, 6.5, 2, 2),  oma = c(0,0,0,0) )
plot(cc,prob_tn_noproj_MC,type = "l",col="blue3",lty=1,cex=2.5,cex.axis=2.5,lwd=4,cex.lab=2.5,ylab=expression('P('*hat(T)["n"]<='c)'),xlab="c")
lines(cc,prob_tn_proj,type = "l", lty=2, col = "grey",cex=2,lwd=4)
lines(cc,prob_tn_noproj,type = "l",col="orange",lty=3,lwd=4)
legend("bottomright", legend=c(expression('Monte Carlo'),
                               expression('Projected Bootstrap'),
                               expression('Parametric Bootstrap')), lty=1:3, col = c("blue3","grey","orange"),lwd=4,cex=1.9,bty = "n")
dev.off()

cc=seq(0,15,length=1000)
prob_tn_order_proj <- prob_tn_order_noproj <- prob_tn_order_noproj_MC <- numeric(1000)
for (k in 1:1000){
  prob_tn_order_proj[k] <- sum(res_tn_order_proj<cc[k])/length(res_tn_order_proj)
  prob_tn_order_noproj[k] <- sum(res_tn_order_noproj<cc[k])/length(res_tn_order_noproj)
  prob_tn_order_noproj_MC[k] <- sum(res_tn_order_noproj_MC<cc[k])/length(res_tn_order_noproj_MC)
}

postscript("tnorderstat.eps", width = 8, height = 6, horizontal = FALSE, onefile = FALSE)
par(mar = c(5, 6.5, 2, 2),  oma = c(0,0,0,0))
plot(cc,prob_tn_order_noproj_MC,type = "l",col="blue3",lty=1,cex=2.5,cex.axis=2.5,lwd=4,cex.lab=2.5,ylab=expression('P('*tilde(T["n"])<='c)'),xlab="c")
lines(cc,prob_tn_order_proj,type = "l", lty=2,col="grey",lwd=4)
lines(cc, prob_tn_order_noproj,type = "l",col="orange",lty=3,lwd=4)
legend("bottomright", legend=c(expression('Monte Carlo'),
                               expression('Projected Bootstrap'),
                               expression('Parametric Bootstrap')),lty=1:3,col = c("blue3","grey","orange"),lwd=4,cex=1.9,bty = "n")
dev.off()


