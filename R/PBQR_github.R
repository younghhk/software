###########
#PBQR function for
#A QUANTILE REGRESSION DECOMPOSITION ESTIMATION OF DISPARITIES FOR COMPLEX SURVEY DATA"
#Refer to Hong et al. 2022+
#Inputs: y, X.d (race), X.s (age), M (majority), m (minority),
#a0, a1, (range of the age variable), x (covariates), wt (sampling weight), strat (strata), psu
#Example
library(NHANES)
dat=NHANESraw[NHANESraw$SurveyYr=="2009_10",] #2009-2010 NHANES dataset
gender2=ifelse(dat$Gender=="male",1,2) #male=1; female=2
race2=ifelse(dat$Race1=="White",1,2) #Whites=1; otherwise=2

y=dat$BMI #outcome
X.d=race2 #racial disparity 
X.s=dat$Age #age subgroup
x=data.frame(age=dat$Age, gender=gender2, poverty=dat$Poverty) #age should be included as covariate
wt=dat$WTMEC2YR #sampling weights
strat=dat$SDMVSTRA #strata
psu= dat$SDMVPSU #primary sampling unit
#K=no. of generated taus
#B= no. of permutation, B>=100 recommended

#run (slow, if B is large)
# To decompose  racial disparity (whites vs. others) with age subset (40 and 60yr old) given X=(gender,age, poverty)
out1=PBQR(y,x,X.d=X.d,M=1,m=2,X.s=X.s,a0=40,a1=60,B=100,K=500,wt=wt,psu=psu,strat=strat)
out1

# To decompose racial disparity (whites vs. others) given X=(gender,age, poverty)
out2=PBQR(y,x,X.d=X.d,M=1,m=2,X.s=X.s,a0=min(X.s),a1=max(X.s),B=2,K=500,wt=wt,psu=psu,strat=strat)
out2
###############################


PBQR=function(y, x, X.d, M, m, X.s, a0, a1, B, K, wt, psu, strat){
  library(stats)
  library(survey)
  library(quantreg)
  library(splines)
  library(Hmisc)
  
  newdat=data.frame(y=y,x,wt=wt,psu=psu,strat=strat,X.d=X.d,X.s=X.s)
  newdat=na.omit(newdat)
  taus=seq(.05,.95,.05)
  n.taus=length(taus)
  theta<-seq(1,K,1)/(K+1)
  
  qMM.out=matrix(0,nrow=B,ncol=n.taus) 
  qMm.out=matrix(0,nrow=B,ncol=n.taus) 
  qmm.out=matrix(0,nrow=B,ncol=n.taus) 
  
  temp.data<-newdat[which(newdat$X.d==M),]
  temp.data=temp.data[temp.data$X.s>=a0 &temp.data$X.s<=a1,]

  
  wt.normalized<-temp.data$wt/sum(temp.data$wt) 
  cumsum.wt.normalized<-cumsum(wt.normalized) 
  u.x<-runif(K); 
  temp.index<-NULL
  for (k in 1:K) { temp.index<-c(temp.index, sum(cumsum.wt.normalized<u.x[k])+1) }
  cat1.X<-temp.data[temp.index,]
  
  
  results<-NULL
  temp.data<-newdat[which(newdat$X.d==M),]
  for (k in 1:K) {  
    xx=as.matrix(temp.data[,c(colnames(x))])
    temp<-rq(y~xx, weights=wt, data=temp.data, tau=theta[k])
    results<-rbind(results, temp$coef)
  }
  
  cat1.results<-results
  for (v in 1:ncol(results)){
    ttemp<-lm(cat1.results[,v]~ns(theta, df=10, intercept=F)); 
    cat1.results[,v]<-ttemp$fitted
  }
  

  temp.data<-newdat[which(newdat$X.d==m),]

  
  wt.normalized<-temp.data$wt/sum(temp.data$wt) 
  cumsum.wt.normalized<-cumsum(wt.normalized) 
  u.x<-runif(K); 
  temp.index<-NULL
  for (k in 1:K) { temp.index<-c(temp.index, sum(cumsum.wt.normalized<u.x[k])+1) }
  cat2.X<-temp.data[temp.index,]
  
  
  results<-NULL
  for (k in 1:K) { 
    xx=as.matrix(temp.data[,c(colnames(x))])
    temp<-rq(y~xx, weights=wt, data=temp.data, tau=theta[k])
    results<-rbind(results, temp$coef)
  }
  cat2.results<-results
  for (v in 1:ncol(results)){
    ttemp<-lm(cat2.results[,v]~ns(theta, df=10, intercept=F)); 
    cat2.results[,v]<-ttemp$fitted
  }
  
  
  q11<-apply(as.matrix(cat1.X[,c(colnames(x))])*cat1.results[,-1], 1, sum)+cat1.results[,1]
  q22<-apply(as.matrix(cat2.X[,c(colnames(x))])*cat2.results[,-1], 1, sum)+cat2.results[,1]
  q12<-apply(as.matrix(cat2.X[,c(colnames(x))])*cat1.results[,-1], 1, sum)+cat1.results[,1]
 
  qMM=quantile(q11, prob=taus)
  qmm=quantile(q22, prob=taus)
  qMm=quantile(q12, prob=taus)
  
  D=qMM-qmm
  U=qMm-qmm
  
  res<-rbind(qMM,qmm,qMm,D=D, U=U)
  res=round(res,2)
  rownames(res)=c("q.majority","q.minority","q.counterfactual","overall disparity","unexplained dispairty")

  
  
  ## Permutation-based variance estimation
  for(b in 1:B){
    #set.seed(b)
    R.hj=rep(0, nrow(newdat))
    S=unique(newdat$strat)
    J=unique(newdat$psu)
    for (h in 1:length(S)){
      for (j in 1:length(J)){
        index=which(newdat$strat==S[h]&newdat$psu==J[j])
        R.hj[index]<-rexp(n=1,rate=1)
    
      }
    }
    
    newdat$R.hj<-R.hj
    new.wt=newdat$R.hj*newdat$wt
    newdat$new.wt=new.wt
    
    #=permute qMM
    temp.data<-newdat[which(newdat$X.d==M),]
    xx=as.matrix(temp.data[,c(colnames(x))])
    results<-NULL
    for (k in 1:K) { 
      temp<-rq(y~xx,weights=new.wt, data=temp.data, tau=theta[k])
      results<-rbind(results, temp$coef)
    }
    cat1.results<-results;
    for (v in 1:ncol(results)){
      ttemp<-lm(cat1.results[,v]~ns(theta, df=10, intercept=F)); 
      cat1.results[,v]<-ttemp$fitted
    }
    
    temp.data=temp.data[temp.data$X.s>=a0 &temp.data$X.s<=a1,]
    wt.normalized<-temp.data$wt/sum(temp.data$wt)
    cumsum.wt.normalized<-cumsum(wt.normalized);
    u.x<-runif(K); 
    temp.index<-NULL
    for (k in 1:K) { temp.index<-c(temp.index, sum(cumsum.wt.normalized<u.x[k])+1) }
    cat1.X<-temp.data[temp.index,]
    
    q11<-apply(as.matrix(cat1.X[,c(colnames(x))])*cat1.results[,-1], 1, sum)+cat1.results[,1]
    qMM.out[b,]<-quantile(q11, prob=taus)
    
    ##permute qMm 
    temp.data<-newdat[which(newdat$X.d==m),]
    temp.data=temp.data[temp.data$X.s>=a0 &temp.data$X.s<=a1,]
    wt.normalized<-temp.data$wt/sum(temp.data$wt)
    cumsum.wt.normalized<-cumsum(wt.normalized)
    u.x<-runif(K); 
    temp.index<-NULL
    for (k in 1:K) { temp.index<-c(temp.index, sum(cumsum.wt.normalized<u.x[k])+1) }
    cat2.X<-temp.data[temp.index,]
    
    q12<-apply(as.matrix(cat2.X[,c(colnames(x))])*cat1.results[,-1], 1, sum)+cat1.results[,1]
    qMm.out[b,]<-quantile(q12, prob=taus)
    
    ##permute qmm 
    results<-NULL
    temp.data<-newdat[which(newdat$X.d==m),]
    xx=as.matrix(temp.data[,c(colnames(x))])
    for (k in 1:K) {  
      temp<-rq(y~xx,weights=new.wt, data=temp.data, tau=theta[k])
      results<-rbind(results, temp$coef)
    }
    cat2.results<-results;
    for (v in 1:ncol(results)){
      ttemp<-lm(cat2.results[,v]~ns(theta, df=10, intercept=F)); 
      cat2.results[,v]<-ttemp$fitted
    }
    
    
    q22<-apply(as.matrix(cat2.X[,c(colnames(x))])*cat2.results[,-1], 1, sum)+cat2.results[,1]
    qmm.out[b,]<-quantile(q22, prob=taus)
  cat("B=",b)
  }
  
  vU=round(apply(qMm.out-qmm.out,2,var),3)
  out=list(res=res,  var.unexplained.disparity=vU)
  out
}

  
  