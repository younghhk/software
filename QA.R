#------------------------------------------------------------------------#
# Created on    :   April 1, 2012  
# AUTHOR        :   Hyokyoung Grace Hong jointly with Xuming He and Lan Wang
# AFFILIATION   :   Michigan State University
# EMAIL         :   hhong@stt.msu.edu
#
# QUANTILE-ADAPTIVE MODEL-FREE VARIABLE SCREENING FOR HIGH-DIMENSIONAL HETEROGENEOUS DATA,
# 2013, The Annals of Statistics
# Input1:y (response), x (covariates), tau (quantiles of interest)
# Input2:time (failure times), delta (censoring indicator), x (covariates), tau (quantiles of interest)
# Output: Quantile adaptive statistics
# R: 2.15.1
#-------------------------------------------------------------------------#
### load packages
library(quantreg)
library(splines)
library(MASS)
library(survival)
#------------------------------------------------------------------------#
### Compute Quantile Adaptive sure Independence Screening (QaSIS)
## for uncensored dataset
QaSIS<-function(y,x,tau)
{
  p=dim(x)[2]
  n=length(y)
  fit<-numeric(p)
  y<-y-quantile(y,tau) #centered y
  x<-scale(x)
  for(j in 1:p){
    x0<-x[,j]
    knots=quantile(x0,c(1/3,2/3))
    a=bs(x0, knots=knots,degree=1)
    b=rq(y~a,tau=tau)
    fit[j] <-sum((b$fitted)^2)/n #avg of f^2 for each j
  }
  return(fit)
}


QaSIS.surv = function(x,time,delta,tau)
{
  N = length(time)
  p = ncol(x)
  survy=survfit(Surv(time,delta)~1)
  medy=survy$time[min(which(survy$surv<tau))]
  
  survc=survfit(Surv(time,1-delta)~1)
  wt=rep(0,N)
  for(i in 1:N)
  {
    k=max(which(round(survc$time,4)<=round(time[i],4)))
    wt[i]=delta[i]/survc$surv[k]
  }
  
fit=numeric(p)
  for(k in 1:p)
  {
    pix=bs(x[,k],df=3)
    betac=rq(time~pix, tau=tau, weight=wt)
    fit[k]=mean((predict(betac)-medy)^2)
  }
  return(fit)
}

#########################################
## Example for uncensored dataset
 simul_dat_example= function(N, p,seed=100){
  if(!is.null(seed))
  set.seed(seed)
  active=c(1:4)
  Sigma1=diag(p) ## covariance matrix
  for (i in 1:p){
    for (j in 1:p){
      if (i<j) {
        Sigma1[i,j]=(0.8)^abs(i-j)
        Sigma1[j,i]=Sigma1[i,j]
      }
    }
  }
  g1<-function(x) (x)
  g2<-function(x) ((2*x-1)^2)
  g3<-function(x) (sin(2*pi*x)/(2-sin(2*pi*x)))
  g4<-function(x)(0.1*sin(2*pi*x)+0.2*cos(2*pi*x)+0.3*sin(2*pi*x)^2+0.4*cos(2*pi*x)^3+0.5*sin(2*pi*x)^3)
  X=mvrnorm(N,mu=rep(0,p),Sigma=Sigma1)
  Y=5*g1(X[,1])+3*g2(X[,2])+4*g3(X[,3])+6*g4(X[,4])+sqrt(1.74)*rnorm(N,0,1) 
  return(list(x=X,y=Y,active=active))
}

dat=simul_dat_example(N=200,p=1000, seed=100) 
x=dat$x
y=dat$y
active=dat$active
## QaSIS(tau=.5)
out<-QaSIS(y=y,x=x,tau=.5)
rank(-out)[active] #rank of active variables
  
