#CODE
#####################################
# Note (last update: 7/28/2020)
# tvpSIR.R
# run the time-varing Poisson SIR model (tvpSIR.R).
# the goal is to obtain predicted  R0(t)

#Example
out=tvpSIR(Co="Korea",esp=.5,P=0,p=8,q=5)
T=nrow(out$all)
data.frame(out$all[,1], out$all$R0)[(T-14):T,]
########################################
#Required the following
library(splines)
library("numDeriv")
#Co:select country
#P:prediction interval
#p,q: the number of knots
#esp: regression penalty

tvpSIR=function(Co,esp,P,p,q){
  data=choose.country(Co)
  dat=data$k
  date=dat[,1]
  N=data$N
  head(dat)
  a=which(dat[,3]!=0)
  b=which(dat[,2]!=0)
  first.date=intersect(a,b)[1]  ;first.date 
  #first.date=1
  rr=dat[first.date:nrow(dat),3]
  ii=dat[first.date:nrow(dat),2]
  r0=rr[1]/N;r0
  i0=ii[1]/N;i0
  Z.r=rr/N
  Z.i=ii/N

  T=length(Z.r)-P #duration for estimation;remaining for validation
  times=1:T
  spdeg=3
  


  knots1 <- seq(from=min(times),to=max(times),length=p+2)[2:(p+1)]
  knots2<-seq(from=min(times),to=max(times),length=q+2)[2:(q+1)]
  result0=nlm(NegLogcst, c(0,0), Z.r[1:T], Z.i[1:T],r0,i0, hessian=F, times,N, iterlim=1000)
  initial= c( rep(result0$estimate[1], length(knots1)+spdeg+1), rep(result0$estimate[2], length(knots2)+spdeg+1)) 
  
  
  result=nlm(NegLoglik, p=initial, Z.r[1:T], Z.i[1:T],r0,i0, hessian=T, times,N, eps, knots1=knots1,knots2=knots2, 
             spdeg=spdeg, steptol=1e-7, iterlim=1000)
  
  x=result$estimate 
  
  
  est=SIR(times, x, r0,i0, knots1,knots2, spdeg)
  

  ##estimate CI of I and R
  remove=numeric(T)
  infect=numeric(T)
  lo.infect=numeric(T)
  up.infect=numeric(T)
  lo.remove=numeric(T)
  up.remove=numeric(T)
  
  COV = solve(result$hessian)
  #print(COV)
  for (i in 1:T) 
  {
    jacob = jacobian(RI, x=x, r0=r0, i0=i0, time=i, knots1=knots1,
                     knots2=knots2,spdeg=spdeg,times=times)
    
    
    # se.remove= sqrt( t(jacob[1,])%*%diag(COV)%*%(jacob[1,]) )
    se.remove= sqrt(  sum(  jacob[1,]^2 * diag(COV)) )
    #  se.infect= sqrt( t(jacob[2,])%*%diag(COV)%*%(jacob[2,]) )
    se.infect= sqrt(  sum(jacob[2,]^2 * diag(COV)) )
    
    RIcall = RI(x, r0, i0, i,knots1=knots1,
                knots2=knots2, spdeg, times)
    remove[i]= RIcall[1]
    infect[i]= RIcall[2]
    #cat(remove[i],se.remove,  infect[i], se.infect, "\n")
    lo.infect[i]= (infect[i] -1.96*se.infect)
    up.infect[i]= (infect[i] +1.96*se.infect)
    lo.remove[i]= (remove[i] -1.96*se.remove)
    up.remove[i]= (remove[i] +1.96*se.remove)
  }
  
  
 
  
  newtime=date[first.date:(nrow(dat)-P)]
  all0= N*cbind( Z.r[1:T], est$r,  Z.i[1:T], est$i, lo.infect, up.infect,lo.remove, up.remove)
  all=data.frame(newtime=newtime, est$beta,est$gamma,est$R0, all0)
  colnames(all)=c("newtime", "beta","gamma","R0","obs.r","est.r","obs.i","est.i","lo.i","up.i","lo.r","up.r")
  
  
  
 
  func.beta = splinefun(x=times, y=est$beta, method="monoH.FC",  ties = mean)
  func.gamma = splinefun(x=times, y=est$gamma, method="monoH.FC",  ties = mean)
  func.r = splinefun(x=times, y=est$r, method="monoH.FC",  ties = mean)
  func.i = splinefun(x=times, y=est$i, method="monoH.FC",  ties = mean)
  ftimes=(T+1):(T+P)
  pred.r=func.r(ftimes)
  pred.i=func.i(ftimes)
  pred.beta=func.beta(ftimes)
  pred.gamma=func.gamma(ftimes)
  pred.R0=pred.beta/pred.gamma
  
  pred=data.frame(newtime=date[(nrow(dat)-P+1):nrow(dat)],pred.beta,pred.gamma, pred.R0, Z.r[(T+1):(T+P)], pred.r,  Z.i[(T+1):(T+P)], pred.i, lo.infect=NA, up.infect=NA,lo.remove=NA, 
                  up.remove=NA)
  colnames(pred)=colnames(all)
  list(all=all,pred=pred )
}

SIR<- function(times, x, r0,i0, knots1,knots2,spdeg){
  T=length(times)
  i=numeric(T)
  r=numeric(T)
  i[1]=i0
  r[1]=r0

  
 
  t0.beta<-bs(times,  knots =knots1,  degree = spdeg, intercept =TRUE,Boundary.knots=c(-20,T+3))
  t0.gamma<-bs(times,  knots =knots2,  degree = spdeg, intercept =TRUE,Boundary.knots=c(-20,T+3))
  k=dim(t0.beta)[2]

  
  beta = exp( t0.beta%*% x[1:k]) 
  gamma = exp(t0.gamma%*% x[-(1:k)]) 
  
  for(j in 1: (T-1)){
    r[j+1]=r[j]+gamma[j]*i[j]
    i[j+1]=i[j]+beta[j]*(1-r[j]-i[j])*i[j]-gamma[j]*i[j] 
  }
  R0=beta/gamma
  all = cbind(times, r,i,beta,gamma,R0)
  colnames(all)=c("times", "r","i","beta","gamma","R0") #revised to add times: helps for the interpolation
  all <- as.data.frame(all)
  return(all)
}


##### constant SIR
SIR0<- function(times, x, r0,i0){
  T=length(times)
  i=numeric(T)
  r=numeric(T)
  i[1]=i0
  r[1]=r0
  
  beta = exp(x[1]) 
  gamma = exp(x[2]) 
  
  for(j in 1: (T-1)){
    r[j+1]=r[j]+gamma*i[j]
    i[j+1]=i[j]+beta*(1-r[j]-i[j])*i[j]-gamma*i[j] 
  }
  
  all = cbind(r,i)
  colnames(all)=c("r","i")
  all <- as.data.frame(all)
  return(all)
  
}


################## likelihood with constant 
NegLogcst=function(x,Z.r,Z.i,r0,i0,times,N){
  z=SIR0(times, x, r0, i0)
  R=z$r
  I=z$i
  I=z$i
  Loglik=N*sum(-R+Z.r*log(R) -I+Z.i*log(I) )   
  return(-Loglik)
}



###########################################
NegLoglik=function(x,Z.r,Z.i,r0,i0,times,N, knots1,knots2,spdeg,eps){
  z=SIR(times, x, r0, i0,knots1,knots2, spdeg)
  R=z$r
  I=z$i
  Loglik=N*sum(-R+Z.r*log(R) -I+Z.i*log(I) )   
  T=length(times)
  eps=  eps*T  ## important to add this term to ensure the gamma will not be estimated to be 0
    return(-Loglik + eps*t(x)%*%x)
}

#####################################
RI=function(x, r0,i0, time, knots1,knots2, spdeg, times){
  z=SIR(times, x, r0,i0,knots1,knots2, spdeg)
  remove=z$r[time]
  infect=z$i[time]
  return(c(remove,infect))
}

###########################
choose.country=function(Country){
  if(Country=="Korea"){
    data=read.csv("https://raw.githubusercontent.com/ulklc/covid19-timeseries/master/countryReport/country/KR.csv")
    k=make.dat(data)
    N=	51269185
  }
  if(Country=="Sweden"){
    data=read.csv("https://raw.githubusercontent.com/ulklc/covid19-timeseries/master/countryReport/country/SE.csv")
    k=make.dat(data)
    N=10099265
  }
 

  if(Country=="US"){
    data=read.csv("https://raw.githubusercontent.com/ulklc/covid19-timeseries/master/countryReport/country/US.csv")
    k=make.dat(data)
    N=331002651
    }
  list(k=k,N=N)
}

make.dat=function(dat){
  case=as.numeric(dat[,7])
  recovered=as.numeric(dat[,8])
  death=as.numeric(dat[,9])
  case=case-recovered-death
  removal=death+recovered
  date=dat[,1]
  k=data.frame(date=date,case=case,removal=removal)
  return(k)
}
