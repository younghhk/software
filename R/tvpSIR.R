#CODE

#tvpSIR.R

#Here is an example code that you can copy/paste into your R envriromentto 
#run the time-varing Poisson SIR model (tvpSIR.R).
# the goal is to obtain predicted I(t), R(t) and R0(t)

#required library
library(splines)
library("numDeriv")

#import the following sources
#source("SIR0.r")
#source("SIR.r")
#source("NegLogcst.r")
#source("NegLoglik.r")
#source("RI.r")


#####################################
# Note (last update: 6/3/2020)

#Inputs: Data, number.knots, N

##Data have the following time-series dataframe format consistting of [time,i,r], which is a Tx3 matrix,
##where 
#time =1,...,T
##r is the fraction of removed
##i is the fraction of infectious
##number.knots is the number of inner knots for cubic B-spline
##N is the population size

#Outputs:
#estimated beta(t), gamma(t), R_0(t), I(t), R(t), and lower and upper bounds of I(t) and R(t).
#####################################################################
#beginning of the code
tvpSIR=function(Data, number.knots, N){
  r=dat$r
  i=dat$i
  r0=r[1]/N #initial condition
  i0=i[1]/N #initial condition
  
  ##Note: for stable results, set t=0 as the date when both r0 and i0 are nonzero 
  Z.r=r/N 
  Z.i=i/N
  
  T=length(Z.r)
  times=1:T
  spdegree=3
  
  
  ##esimtate the initial parameters for the likelihood
  knots <- seq(from=min(times),to=max(times),length=number.knots+2)[2:(number.knots+1)]
  result0=nlm(NegLogcst, c(0,0), Z.r, Z.i,r0,i0, hessian=F, times,N, iterlim=1000)
  initial= c( rep(result0$estimate[1], length(knots)+spdegree+1), rep(result0$estimate[2], length(knots)+spdegree+1)) 
  
  ##obtain the B-spline coefficients using N-R
  result=nlm(NegLoglik, initial, Z.r, Z.i,r0,i0, hessian=T, times,N,  knots=knots, spdeg=spdegree, steptol=1e-7, iterlim=1000)
  x=result$estimate 
  
  ##estimate  beta,gamma, R0, I(t),R(t)
  est.dat=SIR(times, x, r0,i0, knots, spdegree)
  
  ##construct CI of I(t) and R(t)
  remove=numeric(T)
  infect=numeric(T)
  lo.infect=numeric(T)
  up.infect=numeric(T)
  lo.remove=numeric(T)
  up.remove=numeric(T)
  
  COV = solve(result$hessian)
  
  for (i in 1:T) 
  {
    jacob = jacobian(RI, x=x, r0=r0, i0=i0, time=i, knots=knots,spdeg=spdegree,times=times)
    se.remove= sqrt(sum(jacob[1,]^2 * diag(COV)) )
    se.infect= sqrt(sum(jacob[2,]^2 * diag(COV)) )
    
    RIcall = RI(x, r0, i0, i, knots, spdegree, times)
    remove[i]= RIcall[1]
    infect[i]= RIcall[2]
    
    lo.infect[i]= (infect[i] -1.96*se.infect)
    up.infect[i]= (infect[i] +1.96*se.infect)
    lo.remove[i]= (remove[i] -1.96*se.remove)
    up.remove[i]= (remove[i] +1.96*se.remove)
  }
  
  
  out0= N*cbind(Z.r, est.dat$r,  Z.i, est.dat$i, lo.infect, up.infect,lo.remove, up.remove)
  out=data.frame(times=times, est.dat$beta,est.dat$gamma,est.dat$R0, out0)
  colnames(out)=c("times", "est.beta","est.gamma","est.R0","obs r","est r","obs i","est i","lo.i","up.i","lo.r","up.r")
  
  list(out=out)
}

#end of the code
##########################################
#Following functions need to be sourced
##################################################
## constant SIR
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


##############################################
##time-varying SIR
SIR<- function(times, x, r0,i0, knots, spdeg){
  T=length(times)
  #beta=numeric(T)
  #gamma=numeric(T)
  i=numeric(T)
  r=numeric(T)
  i[1]=i0
  r[1]=r0
  
  
  t0<-bs(times,  knots =knots, degree = spdeg, intercept =TRUE,Boundary.knots=c(-3,T+3))
  k=dim(t0)[2]
  
  
  beta = exp( t0%*% x[1:k]) 
  gamma = exp(t0%*% x[-(1:k)]) 
  
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

################## 
#Likelihood for constant model
NegLogcst=function(x,Z.r,Z.i,r0,i0,times,N){
  z=SIR0(times, x, r0, i0)
  R=z$r
  I=z$i
  I=z$i
  Loglik=N*sum(-R+Z.r*log(R) -I+Z.i*log(I) )   
  return(-Loglik)
}



###########################################
#Likelihood for time-varying model
NegLoglik=function(x,Z.r,Z.i,r0,i0,times,N, knots, spdeg){
  
  z=SIR(times, x, r0, i0,knots, spdeg)
  R=z$r
  I=z$i
  Loglik=N*sum(-R+Z.r*log(R) -I+Z.i*log(I) )   
  T=length(times)
  eps=  0.5*T  
  return(-Loglik + eps*t(x)%*%x)
}

#################################################
RI=function(x, r0,i0, time, knots, spdeg, times){
  z=SIR(times, x, r0,i0,knots,spdeg)
  remove=z$r[time]
  infect=z$i[time]
  return(c(remove,infect))
}
