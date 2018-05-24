#------------------------------------------------------------------------#
# Created on    :   January 1, 2018  
# AUTHOR        :   Hyokyoung G. Hong
# AFFILIATION   :   Michigan State University
# EMAIL         :   hhong@stt.msu.edu

## Conditional screening for survival data
## Hong, H. G., Kang, J., Li, Y. (2018) 
## Lifetime Data Analysis, 24, 45-71

## Inputs:
## time: time to death or time to censored
## delta: censoring indicator
## x: covariates
## Cset (conditioning set)
## Outputs: screening statistics
############################################################################################
library(glmnet)
library(mvtnorm)
library(survival)

##1) The partial likelihood approach
CS.plik<- function(x,delta,time, Cset){
p = ncol(x)
plik = numeric(p) ###log partial likelihood
fullset=1:p
scanset=fullset[-Cset]
one_model = function(j){
Xtemp=cbind(x[,Cset],x[,j]);
 fit1 = try(coxph(Surv(time,delta)~Xtemp))
      if(inherits(fit1,"try-error")){
          res= -100000
      }
      else res = fit1$loglik[2] 
}
plik[scanset] = sapply(scanset,one_model)
plik[Cset] = Inf
return(plik)
}

##2) Wald test approach
CS.wald<- function(x,delta,time, Cset){
p = ncol(x)
wald = numeric(p) 
fullset=1:p
scanset=fullset[-Cset]
one_model = function(j){
Xtemp=cbind(x[,Cset],x[,j]);
 fit1 = try(coxph(Surv(time,delta)~Xtemp))
      if(inherits(fit1,"try-error")){
          res= -100000
      }
      else res = abs(summary(fit1)$coef[length(Cset)+1,4]) 
}
wald[scanset]=sapply(scanset,one_model)
wald[Cset] = Inf
return(wald)
}

##3) The marginal magnitude approach
CS.mple <- function(x,delta,time, Cset){
  p = ncol(x)
  screen_stat = numeric(p) 
  fullset=1:p
  scanset=fullset[-Cset]
  one_model=function(j){
    Xtemp=cbind(x[,Cset],x[,j]);
    fit1 = try(coxph(Surv(time,delta)~Xtemp))
    if(inherits(fit1,"try-error")){
      res= -100000
    }
    else res = abs(summary(fit1)$coef[length(Cset)+1,1]) 
  }
  screen_stat[scanset] = sapply(scanset,one_model)
  screen_stat[Cset] = Inf
  return(screen_stat)
}

#######################################
##Example
simul_dat_example = function(N, p,seed,c_up_bound=2){
  if(!is.null(seed))
    set.seed(seed)
  active=c(1:6)
  x = matrix(rnorm(N*p,sd=sqrt(0.5)),nrow=N,ncol=p)
  z = rnorm(N,sd=sqrt(0.5))
  x = x + matrix(z,nrow=N,ncol=p)
  
  U=runif(N, 0,1)
  
  true_beta  = rep(0,length=p)
  true_beta[1:5] = 1
  true_beta[6] = -2.5
  dim(true_beta) = c(p,1)
  
  xbeta=x%*%true_beta
  pre_time=-log(U)/(1*exp(xbeta))
  
  pre_censoring=runif(N,1,30)
  pre_censoring=pre_censoring*(pre_censoring<c_up_bound)+c_up_bound*(pre_censoring>=c_up_bound)
  tcens=(pre_censoring<pre_time) # censoring indicator
  delta=1-tcens
  time=pre_time*(delta==1)+pre_censoring*(delta==0)
  
  delta = delta[order(time)]
  x = x[order(time),]
  time = time[order(time)]
  return(list(delta=delta,time=time,x=x,active=active))
}

dat=simul_dat_example(N=400,p=1000, seed=100)  
x=dat$x
time=dat$time
delta=dat$delta
active=dat$active


CS_plik_stat  = CS.plik(x=dat$x,delta=dat$delta,time=dat$time,Cset=c(1))
rank(-CS_plik_stat)[active]
CS_wald_stat = CS.wald(x=dat$x,delta=dat$delta,time=dat$time,Cset=c(1))
rank(-CS_wald_stat)[active]
CS_mple_stat  = CS.mple(x=dat$x,delta=dat$delta,time=dat$time,Cset=c(1))
rank(-CS_MPLE_stat )[active]
  
