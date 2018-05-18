#------------------------------------------------------------------------#
# Created on    :   May 18, 2018  
# AUTHOR        :   Hyokyoung G. Hong
# AFFILIATION   :   Michigan State University
# EMAIL         :   hhong@stt.msu.edu

##Lq-norm learning for ultrahigh-dimensional survival data: an integrative framework
##Statistica Sinica, under revision 2018

## Inputs:
## time: time to death or time to censored
## delta: censoring indicator
## x: covariates
## q: the power of interest 
## Outputs: screening statistics
############################################################################################
library(mvtnorm)
library(survival)
library(MASS)
library(prodlim)
### define some functions
pairwise.difference <- function(m){
npairs <- choose( ncol(m), 2 )
results <- matrix( NA, nc=npairs, nr=nrow(m) )
cnames <- rep(NA, npairs)
if(is.null(colnames(m))) colnames(m) <- paste("col", 1:ncol(m), sep="")
k <- 1
for(i in 1:ncol(m)){
for(j in 1:ncol(m)){
if(j <= i) next;
results[ ,k] <- m[ ,i] - m[ ,j]
cnames[k] <- paste(colnames(m)[ c(i, j) ], collapse=".vs.")
k <- k + 1
}
}
colnames(results) <- cnames
rownames(results) <- rownames(m)
return(results)
}

#For Continuous Xj
Cont.CMC= function(j,x,delta,time,q){
N=nrow(x)
Lambda=(3:round(log(N),0))
out=array(0,dim=c(length(Lambda),length(q)))
for( i in 1: length(Lambda)){##i
R=Lambda[i]
qt=c(quantile(x[,j],probs=c((1:(R-1))/R)))
#Slice Xj into 1,2,..,R
index=rep(1,nrow(x))
for (r in 1:R){##r
ind=which(x[,j]>qt[r])
index[ind]=r+1
}##r
#obtain dS
time.pool=numeric(R)
for(r in 1:R){##r
  time.pool[r]=max(time[index==r])
}
times=seq(min(time),min(time.pool),.1)
T=length(times)
km.all=survfit(Surv(time, delta)~ 1)
km.all.t=summary(km.all, times=times)
km.t<-km.all.t$surv
dS=km.t[-T]-km.t[-1]
#obtain S(t|R=r)
temp.km=array(0,dim=c(length(times),R))
for(r in 1:R){##r
sub.time=time[index==r]; sub.delta=delta[index==r]
km=survfit(Surv(sub.time, sub.delta)~ 1)
km.at.time.t=summary(km, times=times)
		temp.km[,r]<-km.at.time.t$surv
}##r
for (qq in 1:length(q)){ 
power=q[qq]	
metric=abs(pairwise.difference (temp.km))^power
out[i,qq]=max(apply(cbind(metric[-T,]*dS),2,sum)^(1/power))
}##q
}##i
Out=apply(out,2,sum)
return(Out)
}# end of cont.CMC function


### For discrete Xj
disc.CMC= function(j,x,delta,time,q){
N=nrow(x)
Out=numeric(length(q))
R=sort(unique(x[,j]))
##define the time interval
time.pool=numeric(length(R))
for(r in 1:length(R)){##r
  time.pool[r]=max(time[x[,j]==R[r]])
}
times=seq(min(time),min(time.pool),.1)
T=length(times)
##obtain dS=S(t)-S(t-1)
km.all=survfit(Surv(time, delta)~ 1)
km.all.t=summary(km.all, times=times)
km.t<-km.all.t$surv
dS=km.t[-T]-km.t[-1]
#get conditional S(t|xj)
temp.km=array(0,dim=c(length(times),length(R)))
for(r in 1:length(R)){##r
sub.time=time[x[,j]==R[r]]; sub.delta=delta[x[,j]==R[r]]
km=survfit(Surv(sub.time, sub.delta)~ 1)
km.at.time.t=summary(km, times=times)
temp.km[,r]<-km.at.time.t$surv
}
for (qq in 1:length(q)){
power=q[qq]	
metric=abs(pairwise.difference(temp.km))^power
Out[qq]=max(apply(cbind(metric[-T,]*dS),2,sum)^(1/power))
}
return(Out)
}#function

###Wrapper function
Lq.stat=function(x,delta,time,q){
x=scale(x)
p = ncol(x)
cep=numeric(p)
one_model=function(j){
	if( length(unique(x[,j]))<10) {res=disc.CMC(j,x,delta,time,q)
	} else{res=Cont.CMC(j,x,delta,time,q)
		}
		}
cep = sapply(1:p,one_model)
return(cep)
}

###############################################
###Example
###############################################
simul_dat_example= function(N, p,seed=100){
  if(!is.null(seed))
    set.seed(seed)
  active=c(1,2)
  x = matrix(rnorm(N*p,sd=1),nrow=N,ncol=p)
  z = rnorm(N,sd=sqrt(0.5))
  x = x + matrix(z,nrow=N,ncol=p)
  true_beta  = rep(0,length=p)
  true_beta[1:2]=1
  U=runif(N)
  pre_time=(-log(U)/(abs(x[,1])+abs(x[,2])))^(1/2)
  pre_censoring=runif(N, 0,3.8) 
  delta=(pre_censoring>pre_time) # censoring indicator
  time=pre_time*(delta==1)+pre_censoring*(delta==0)
  sum(1-delta)/N
  delta = delta[order(time)]
  x = x[order(time),]
  time = time[order(time)]
  return(list(x=x,time=time,delta=delta,active=active))
}

dat=simul_dat_example(N=200,p=1000, seed=100)
x=dat$x
time=dat$time
delta=dat$delta
active=dat$active
q=c(1,2,5,10,100) #note when q>100, use KS instead
#perform La-statistics; when x is continuous, it gives the fused Lq-norm results
res=Lq.stat(x,delta,time,q=q)
rownames(res)=as.character(q)
##The rank of active covariates
apply(-res,1,rank)[active,]
