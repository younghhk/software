#------------------------------------------------------------------------#
# Created on    :   July 19, 2017  
# AUTHOR        :   Hyokyoung G. Hong
# AFFILIATION   :   Michigan State University
# EMAIL         :   hhong@stt.msu.edu

# Integrated Powered Density: Screening Ultrahigh Dimensional
# Covariates with Survival Outcomes, Biometrics, 2018

## Inputs:
## time: time to death or time to censored
## delta: censoring indicator
## x: covariates
## gamma: the power of interest 
## note: bandwidth set as (90th percentile of time- min(time)) /10
## Outputs: screening statistics
# R: 3.1.1
#---------------------------------------------------------#
#load library
library(survPresmooth)
library(MASS)

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


IPOD.cont= function(j,x,delta,time,gamma){
  #tau=quantile(time, prob=.9)

  N=nrow(x); Lambda=(3:round(log(N),0))
  out=array(0,dim=c(length(gamma), length(Lambda)))
  
  for( i in 1: length(Lambda)){##i
    R=Lambda[i] #of slicings, R=3,4,5,6,...
    q=c(quantile(x[,j],probs=c((1:(R-1))/R),na.rm = TRUE))
    
    ##create index for subgroup
    index=rep(1,nrow(x))
    for (r in 1:length(q)){
      ind=which(x[,j]>=q[r])
      index[ind]=r+1
    }
    
    time.pool=numeric(R)
    for(r in 1:R){##r
      time.pool[r]=max(time[index==r])
    }
    t=seq(min(time),min(time.pool),.1)
    h=c(0,diff(range(t))/10)
    
    temp.a=array(0,dim=c(length(t)-1,length(gamma),R))
    R.result=numeric(length(gamma))
    for(r in 1:R){##r
      sub.time=time[index==r]; sub.delta=delta[index==r]
      f=presmooth(sub.time, sub.delta, x.est=t, estimand = "f", bw.selec="fixed",fixed.bw =h)
      
      for(a in 1:length(gamma)){ ##a
        g=f$estimate^gamma[a]
        rec=(t[-1]-t[-length(t)])*g[-length(t)]
        temp.a[,a,r]<- cumsum(rec[1:length(rec)])
      } ##a
    }##r
    
    for( a in 1:length(gamma)){#a
      R.result[a]=max(apply(abs(pairwise.difference (temp.a[,a,])),1,max))
    }##a
    out[,i]=R.result
  }##i
  out=apply(out,1,sum)
  return(out)
}#function

###discrete covariate
IPOD.disc= function(j,x,delta,time,gamma){
  R=sort(unique(x[,j]))
  time.pool=length(R)
  for(r in 1:length(R)){##r
    time.pool[r]=max(time[x[,j]==R[r]])
  }
  t=seq(min(time),min(time.pool),.1)
  h=c(0,diff(range(t))/10)
  out=numeric(length(gamma))
  temp.a=array(0,dim=c(length(t)-1,length(gamma),length(R)))
               for(r in 1:length(R)){##r
                 sub.time=time[x[,j]==R[r]]; sub.delta=delta[x[,j]==R[r]]
                               f=presmooth(sub.time, sub.delta, x.est=t, estimand = "f", bw.selec="fixed",fixed.bw =h)
                               
                               for(a in 1:length(gamma)){ ##a
                                 g=f$estimate^gamma[a]
                                 rec=(t[-1]-t[-length(t)])*g[-length(t)]
                                 temp.a[,a,r]<- cumsum(rec[1:length(rec)])
                               }##a
               }##r
               for( a in 1:length(gamma)){#a
                 out[a]=max(apply(abs(pairwise.difference (temp.a[,a,])),1,max))
               }##a
               return(out)
} 

### wrapper function for continuous and discrete covariates
IPOD=function(x,delta,time,gamma){
  p = ncol(x)
  cep=numeric(p)
  one_model=function(j){
    if( length(unique(x[,j]))<5) {res=IPOD.disc(j,x,delta,time,gamma)
    }else  {res=IPOD.cont(j,x,delta,time,gamma)
    }
  }
  cep = sapply(1:p,one_model)
  return(cep)
}

#======================================#
##Example from Case I of Li et al. (2016), Biometrics
simul_dat_example= function(N, p,rho, seed=100){
  if(!is.null(seed))
    set.seed(seed)
  mu0 = rep(0,p);
  sigij = function(i,j)
  { return(rho^abs(i-j)) }
  sigMa0 = matrix(0,p,p)
  
  for(i in 1:p)
  { for(j in 1:p)
  { sigMa0[i,j] = sigij(i,j)
  }
  }
  X = mvrnorm(N,mu0, sigMa0)
  eps = rnorm(N)
  # calculate the response Y
  ## define the calculation funtion of Y from x for Case 1
  g1 = function(x)
  { return(5*x) }
  
  g2 = function(x)
  { return(-4*x*(1-x)) }
  
  g3 = function(x)
  { return(10*(exp(-(3*x-1)^2)+exp(-4*(x-3)^2))-1.5) }
  
  g4 = function(x)
  { return(4*sin(2*pi*x)) }
  
  gY = function(x1,x2,x3,x4,eps)
  { return(g1(x1)+g2(x2)+g3(x3)+g4(x4)+eps) }
  
  Y = gY(X[,1],X[,2],X[,3],X[,4],eps)
  
  CeY = rnorm(N,0,4)-rnorm(N,5,1)+0.5*rnorm(N,25,1)
  
  active=1:4
  deltay = as.numeric(Y<=CeY); deltax = rep(1,N)
  obsY = apply(cbind(Y,CeY),1,min);
  
  delta = deltay[order(obsY)]
  x = X[order(obsY),]
  time = obsY[order(obsY)]
  return(list(x=x,time=time,delta=delta,active=active))
}
##==================================#
## Perform IPOD
p=1000; N=300 
dat=simul_dat_example(N, p,rho=0.8, seed=100)
x=dat$x;time=dat$time;delta=dat$delta;active=dat$active
gamma=c(1,1.2) #gamma=1 correspondes to K-S statistics
out=IPOD(x,delta,time,gamma=gamma)
rownames(out)=as.character(gamma)
apply(-out,1,rank)[active,] #the rank of active variables for gamma's
#out[1,] #The screening stat of p covariates for gamma=1
