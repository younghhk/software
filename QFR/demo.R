library(np)
library(MASS)
library(Rglpk)
library(quantreg)
library("parallel")
library("magrittr")


library(quantreg)
library(splines)
library(MASS)
library(survival)
library(tictoc)

# quantile check function
check <- function(x,tau)
{
  return((tau - as.numeric(x<=0))* x)
}

denzero<-function(vec)
{
  return(sum(vec!=0))
}

wherezero<-function(vec)
{
  return(which(vec!=0))
}

# derivative of SCAD function
Pprime<-function(theta,lambda)
{
  a<-3.7
  y=(theta <= lambda)*lambda + (a*lambda - theta)/(a-1) *((theta>lambda)&(theta<a*lambda)) + 0*(theta>a*lambda)
  return(y)
}

# Quantile regression solver
qreg<-function(Ymat,Xmat,tau) # no intercept
{
  p2<-dim(Xmat)[2]
  n2<-length(Ymat)
  obj <- c(rep(tau,n2),rep((1-tau),n2),rep(0,2*p2))
  mat <- cbind(diag(n2),-diag(n2),Xmat,-Xmat)
  dir <- rep("==",n2)
  rhs <- Ymat
  max <- FALSE
  
  opti<-Rglpk_solve_LP(obj,mat,dir,rhs,max=max)$solution
  betaqr = opti[(2*n2+1):(2*n2+p2)]-opti[(2*n2+p2+1):(2*n2+2*p2)]
  return(betaqr)
} 


# BIC selection procedure in the survival setting
selection_bic_surv<-function(Y,X,l_max,l_width,Cn=log(p),delta,T,tau){
  
  # here Y= min(Y,T)
  # choice for an initial estimator
  n=dim(X)[1]
  p=dim(X)[2]
  
  survc=survfit(Surv(T,1-delta)~1)
  wt=rep(0,n)
  for(i in 1:n)
  {
    k=max(which(round(survc$time,4)<=round(T[i],4)))
    wt[i]=delta[i]/survc$surv[k]
  }
  wt_index = which(as.matrix(wt > 0.0001) + 0==1)
  
  beta_ini<-rep(0,p) # zero initial
  IterN<-5 # how many updates? - zero initial
  
  #beta_ini<-qreg(Y,X,tau) # one step
  #IterN<-1 # how many updates? - one step 
  
  # SCAD solution path
  lam_grid<-seq(0,l_max,l_width)
  beta_iscad<-matrix(0,length(beta_ini),length(lam_grid))
  
  ##############################################################################
  # SCAD                                                                       #
  ##############################################################################   
  
  ################################
  # path generation              #
  ################################
  
  # path - beta_si
  
  iter_func = function(i,iter,beta_ini,beta_iscad,lam_grid,tau,n,p,X,Y){  
    if (iter == 1) {beta_i<-beta_ini}
    if (iter !=1) {beta_i<-beta_iscad[,i]}
    
    lambda<-lam_grid[i]
    lamvec= n*Pprime(abs(beta_i),lambda);
    
    obj <- c(rep(tau,n),rep((1-tau),n),lamvec,lamvec)
    mat <- cbind(diag(n),-diag(n),X,-X)
    dir <- rep("==",n)
    rhs <- Y
    max <- FALSE
    
    opt<-Rglpk_solve_LP(obj,mat,dir,rhs,max=max)$solution
    betaq <- opt[(2*n+1):(2*n+p)]-opt[(2*n+p+1):(2*n+2*p)]
    return(betaq)}
  
  for (i in 1:n){
    X[i,]=X[i,]*wt[i]
  }
  X=X[wt_index,]
  Y=Y*as.matrix(wt)
  Y=Y[wt_index]
  n=length(Y)
  
  for ( iter in 1:IterN)
  {
    cat("ITER #",iter,"\n")
    
    results = mclapply(1:length(lam_grid),iter_func, iter,beta_ini,beta_iscad,lam_grid,tau,n,p,X,Y)
    #results = mclapply(1:2,iter_func, iter,beta_ini,beta_iscad,lam_grid,tau,n,p,X,Y)
    beta_iscad = matrix(unlist(results), byrow=F, nrow=p)
    # beta_iscad[,i]<-betaq  
    #tl<-paste("Update stage=",as.character(iter),"  Sample",as.character(Siter))
    #matplot(lam_grid,t(beta_iscad),type='l',main=tl)
  }
  
  
  si_ind<-which(apply(beta_iscad,2,denzero)<=50)
  
  beta_si<-beta_iscad[,si_ind]
  if ( sum(apply(abs(beta_si),2,sum)==0)>0)
  {
    beta_si<-beta_si[,-which(apply(abs(beta_si),2,sum)==0)] # remove the null fit
  }
  
  ################ CASE 2 (Use the path as it is) ############
  
  selm<-beta_si!=0
  selm = as.matrix(selm)
  
  
  for ( i in 1:dim(selm)[2])
  {
    beta_si[selm[,i],i]<-qreg(Y,as.matrix(X[,selm[,i]]),tau) # fit an unpenalized estimator for each selected model in the path
  }
  
  
  BIC_L_scad<-c() # ordinary BIC
  BIC_H_scad<-c() # H-BIC
  
  for ( i in 1:dim(selm)[2])
  {
    BIC_L_scad<-c(BIC_L_scad,log(sum(check(Y - X %*% beta_si[,i],tau))/n) + log(n)*denzero(beta_si[,i])/(2*n))
    BIC_H_scad<-c(BIC_H_scad,log(sum(check(Y - X %*% beta_si[,i],tau))/n) + Cn*log(n)*denzero(beta_si[,i])/(2*n))
  }
  
  sel_L_scad<-which(beta_si[,which.min(BIC_L_scad)]!=0)
  sel_H_scad<-which(beta_si[,which.min(BIC_H_scad)]!=0)
  
  
  return(list(sel_O=sel_L_scad, sel_H=sel_H_scad))  
}



# BIC selection procedure

selection_bic<-function(Y,X,l_max,l_width,Cn=log(p)){
  # choice for an initial estimator
  n=dim(X)[1]
  p=dim(X)[2]
  
  beta_ini<-rep(0,p) # zero initial
  IterN<-5 # how many updates? - zero initial
  
  #beta_ini<-qreg(Y,X,tau) # one step
  #IterN<-1 # how many updates? - one step 
  
  # SCAD solution path
  lam_grid<-seq(0,l_max,l_width)
  beta_iscad<-matrix(0,length(beta_ini),length(lam_grid))
  
  ##############################################################################
  # SCAD                                                                       #
  ##############################################################################   
  
  ################################
  # path generation              #
  ################################
  
  # path - beta_si
  
iter_func = function(i,iter,beta_ini,beta_iscad,lam_grid,tau,n,p,X,Y){  
  if (iter == 1) {beta_i<-beta_ini}
  if (iter !=1) {beta_i<-beta_iscad[,i]}
  
  lambda<-lam_grid[i]
  lamvec= n*Pprime(abs(beta_i),lambda);
  
  obj <- c(rep(tau,n),rep((1-tau),n),lamvec,lamvec)
  mat <- cbind(diag(n),-diag(n),X,-X)
  dir <- rep("==",n)
  rhs <- Y
  max <- FALSE
  
  opt<-Rglpk_solve_LP(obj,mat,dir,rhs,max=max)$solution
  betaq <- opt[(2*n+1):(2*n+p)]-opt[(2*n+p+1):(2*n+2*p)]
  return(betaq)}
  
  
  for ( iter in 1:IterN)
  {
    cat("ITER #",iter,"\n")
    
    results = mclapply(1:length(lam_grid),iter_func, iter,beta_ini,beta_iscad,lam_grid,tau,n,p,X,Y)
    #results = mclapply(1:2,iter_func, iter,beta_ini,beta_iscad,lam_grid,tau,n,p,X,Y)
    beta_iscad = matrix(unlist(results), byrow=F, nrow=p)
      # beta_iscad[,i]<-betaq  
    
    #tl<-paste("Update stage=",as.character(iter),"  Sample",as.character(Siter))
    #matplot(lam_grid,t(beta_iscad),type='l',main=tl)
    
  }
  
  
  si_ind<-which(apply(beta_iscad,2,denzero)<=50)
  
  beta_si<-beta_iscad[,si_ind]
  if ( sum(apply(abs(beta_si),2,sum)==0)>0)
  {
    beta_si<-beta_si[,-which(apply(abs(beta_si),2,sum)==0)] # remove the null fit
  }
  
  ################ CASE 2 (Use the path as it is) ############
  
  selm<-beta_si!=0
  
  for ( i in 1:dim(selm)[2])
  {
    beta_si[selm[,i],i]<-qreg(Y,as.matrix(X[,selm[,i]]),tau) # fit an unpenalized estimator for each selected model in the path
  }
  
  
  BIC_L_scad<-c() # ordinary BIC
  BIC_H_scad<-c() # H-BIC
  
  for ( i in 1:dim(selm)[2])
  {
    BIC_L_scad<-c(BIC_L_scad,log(sum(check(Y - X %*% beta_si[,i],tau))/n) + log(n)*denzero(beta_si[,i])/(2*n))
    BIC_H_scad<-c(BIC_H_scad,log(sum(check(Y - X %*% beta_si[,i],tau))/n) + Cn*log(n)*denzero(beta_si[,i])/(2*n))
  }
  
  sel_L_scad<-which(beta_si[,which.min(BIC_L_scad)]!=0)
  sel_H_scad<-which(beta_si[,which.min(BIC_H_scad)]!=0)
  
  
  return(list(sel_O=sel_L_scad, sel_H=sel_H_scad))  
}



# Stepwise selection procedure in the survival setting


STEPWISE_surv <- function(x, y, eta1 = 0, eta2 = 1, initvars = NULL, tau=0.5, delta, T){
  # here Y= min(Y,T)
  
  doParallel::registerDoParallel(5)
  `%dopar%` <- foreach::`%dopar%`
  `%>%` <- magrittr::`%>%`
  N <- x %>% nrow(); p <- x %>% ncol(); eta1 <- eta1; eta2 <- eta2
  
  survc=survfit(Surv(T,1-delta)~1)
  wt=rep(0,N)
  for(i in 1:N)
  {
    k=max(which(round(survc$time,4)<=round(T[i],4)))
    wt[i]=delta[i]/survc$surv[k]
  }
  
  for (i in 1:N){
    x[i,]=x[i,]*wt[i]
  }
  y=y*wt
  wt_index = which(wt>0.00001)
  x=x[wt_index,]
  y=y[wt_index]
  
  
  
  EBIC <- min(N,p) %>% numeric(); S <- min(N,p) %>% numeric(); x <-  as.matrix(x)
  
  M <- foreach::foreach (i = c(1:dim(x)[2]), .combine = c, .packages = c("MASS", "magrittr","quantreg","Rglpk")) %dopar% {
     
      check <- function(x,tau){return((tau - as.numeric(x<=0))* x)}
                                                                           
      qreg<-function(Ymat,Xmat,tau) # no intercept
        {
          p2<-dim(Xmat)[2]
          n2<-length(Ymat)
          obj <- c(rep(tau,n2),rep((1-tau),n2),rep(0,2*p2))
          mat <- cbind(diag(n2),-diag(n2),Xmat,-Xmat)
          dir <- rep("==",n2)
          rhs <- Ymat
          max <- FALSE
                                                                             
          opti<-Rglpk_solve_LP(obj,mat,dir,rhs,max=max)$solution
          betaqr = opti[(2*n2+1):(2*n2+p2)]-opti[(2*n2+p2+1):(2*n2+2*p2)]
          return(betaqr)} 
                                                                           
        log(sum(check(y - matrix(x[,i],ncol=1) %*% matrix(qreg(y, matrix(x[, i],ncol=1),tau), ncol=1), tau)))}
  
  S[1] <- which.min(M); newset <- S[1]
  
  fit1 = log(sum(check(y - matrix(x[,S[1]],ncol=1) %*% matrix(qreg(y, matrix(x[, S[1]], ncol=1),tau), ncol=1), tau)))
  
  MM <- fit1 
  
  EBIC[1] <- 2*MM + log(p) * (log(N))/N 
  
  for(k in 2:min(N,p)) {
    
    M <- foreach::foreach (i = c(1:dim(x)[2]), .combine = c, .packages = c("MASS", "magrittr","quantreg","Rglpk")) %dopar% {
                                                                             
    check <- function(x,tau)
    {
      return((tau - as.numeric(x<=0))* x)
     }
                                                                             
     qreg<-function(Ymat,Xmat,tau) # no intercept
     {
       p2<-dim(Xmat)[2]
      n2<-length(Ymat)
      obj <- c(rep(tau,n2),rep((1-tau),n2),rep(0,2*p2))
      mat <- cbind(diag(n2),-diag(n2),Xmat,-Xmat)
      dir <- rep("==",n2)
      rhs <- Ymat
      max <- FALSE
      opti<-Rglpk_solve_LP(obj,mat,dir,rhs,max=max)$solution
      betaqr = opti[(2*n2+1):(2*n2+p2)]-opti[(2*n2+p2+1):(2*n2+2*p2)]
      return(betaqr)                                                                           } 
                                                                             
                                                                             
                                                                             
      if (length(unique(c(newset,i)))==1){
      log(sum(check(y - matrix(x[,unique(c(newset,i))],ncol=1) %*% matrix(qreg(y,matrix(x[, unique(c(newset,i))],ncol=1),tau),ncol=1), tau)))
      } else{
        log(sum(check(y - x[,unique(c(newset,i))] %*% matrix(qreg(y,(x[, unique(c(newset,i))]),tau),ncol=1), tau)))
            }                                                                       
                                                                             
     }
    M[newset] <- 10^10; S[k] <- which.min(M); newset <- S[1:k]
    
    
    if (length(newset)==1){
      fit1= log(sum(check(y - matrix(x[,newset],ncol=1) %*% matrix(qreg(y,matrix(x[, newset],ncol=1),tau),ncol=1), tau)))
    } else{
      fit1= log(sum(check(y - x[,newset] %*% matrix(qreg(y,(x[, newset]),tau),ncol=1), tau)))
    }
    
    MM <- fit1 
    
    EBIC[k] <- 2*MM + k*log(p)*(log(N))/N 
    
    if(EBIC[k]>EBIC[k-1] | k == min(N,p))  break 
  }
  
  finalset <- S[1:(k-1)]
  
  if (length(finalset) == 1){
    fit1=   matrix(qreg(y,matrix(x[, finalset],ncol=1),tau),ncol=1)
  } else{
    fit1=   matrix(qreg(y,(x[, finalset]),tau),ncol=1)
  }
  
  list(finalset = finalset,  fit1 = fit1)
}




# Stepwise selection procedure 

STEPWISE <- function(x, y, eta1 = 0, eta2 = 1, initvars = NULL, tau=0.5){
  
doParallel::registerDoParallel(5)
`%dopar%` <- foreach::`%dopar%`
`%>%` <- magrittr::`%>%`
N <- x %>% nrow(); p <- x %>% ncol(); eta1 <- eta1; eta2 <- eta2
EBIC <- min(N,p) %>% numeric(); S <- min(N,p) %>% numeric(); x <-  as.matrix(x)

M <- foreach::foreach (i = c(1:dim(x)[2]), .combine = c, .packages = c("MASS", "magrittr","quantreg",
                                                                           "Rglpk")) %dopar% {
                                                                             
check <- function(x,tau)
{
return((tau - as.numeric(x<=0))* x)}
      
qreg<-function(Ymat,Xmat,tau) # no intercept
{
p2<-dim(Xmat)[2]
n2<-length(Ymat)
obj <- c(rep(tau,n2),rep((1-tau),n2),rep(0,2*p2))
mat <- cbind(diag(n2),-diag(n2),Xmat,-Xmat)
dir <- rep("==",n2)
rhs <- Ymat
max <- FALSE
                                                                               
opti<-Rglpk_solve_LP(obj,mat,dir,rhs,max=max)$solution
betaqr = opti[(2*n2+1):(2*n2+p2)]-opti[(2*n2+p2+1):(2*n2+2*p2)]
return(betaqr)
} 

log(sum(check(y - matrix(x[,i],ncol=1) %*% matrix(qreg(y, matrix(x[, i],ncol=1),tau), ncol=1), tau)))
}
    
S[1] <- which.min(M); newset <- S[1]
    
fit1 = log(sum(check(y - matrix(x[,S[1]],ncol=1) %*% matrix(qreg(y, matrix(x[, S[1]], ncol=1),tau), ncol=1), tau)))
    
MM <- fit1 
    
EBIC[1] <- 2*MM + log(p) * (log(N))/N 
    
for(k in 2:min(N,p)) {
      
M <- foreach::foreach (i = c(1:dim(x)[2]), .combine = c, .packages = c("MASS", "magrittr","quantreg",
                                                                             "Rglpk")) %dopar% {
                                                                               
check <- function(x,tau)
{
return((tau - as.numeric(x<=0))* x)
}
                                                                               
qreg<-function(Ymat,Xmat,tau) # no intercept
{
p2<-dim(Xmat)[2]
n2<-length(Ymat)
obj <- c(rep(tau,n2),rep((1-tau),n2),rep(0,2*p2))
mat <- cbind(diag(n2),-diag(n2),Xmat,-Xmat)
dir <- rep("==",n2)
rhs <- Ymat
max <- FALSE
                                                                                 
opti<-Rglpk_solve_LP(obj,mat,dir,rhs,max=max)$solution
betaqr = opti[(2*n2+1):(2*n2+p2)]-opti[(2*n2+p2+1):(2*n2+2*p2)]
return(betaqr)
} 
                                                                               
        
if (length(unique(c(newset,i)))==1){
log(sum(check(y - matrix(x[,unique(c(newset,i))],ncol=1) %*% matrix(qreg(y,matrix(x[, unique(c(newset,i))],ncol=1),tau),ncol=1), tau)))
} else{
log(sum(check(y - x[,unique(c(newset,i))] %*% matrix(qreg(y,(x[, unique(c(newset,i))]),tau),ncol=1), tau)))
}                                                                       
}
M[newset] <- 10^10; S[k] <- which.min(M); newset <- S[1:k]
      
if (length(newset)==1){
fit1= log(sum(check(y - matrix(x[,newset],ncol=1) %*% matrix(qreg(y,matrix(x[, newset],ncol=1),tau),ncol=1), tau)))
} else{
fit1= log(sum(check(y - x[,newset] %*% matrix(qreg(y,(x[, newset]),tau),ncol=1), tau)))
}
   
MM <- fit1 
      
EBIC[k] <- 2*MM + k*log(p)*(log(N))/N 
      
if(EBIC[k]>EBIC[k-1] | k == min(N,p))  break 
}
    
finalset <- S[1:(k-1)]
    
if (length(finalset) == 1){
fit1=   matrix(qreg(y,matrix(x[, finalset],ncol=1),tau),ncol=1)
} else{
fit1=   matrix(qreg(y,(x[, finalset]),tau),ncol=1)
}
    
list(finalset = finalset,  fit1 = fit1)
}


# quantile adaptive screening procedure
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


# quantile adaptive screening procedure in the survival setting
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


## Simulation analyses

library(mvtnorm)
library(parallel)
sdd=1


## Example: nonsurvival simulation setting

example1_QA = function(repnumb,n,p,sig=0.5,tau){
  sigma = array(0,c(p-1,p-1))
  for (ii in 1:(p-1)){
    for (jj in 1:(p-1)){
      sigma[ii,jj] = sig^(abs(ii-jj))
    }
  }
  set.seed(repnumb)
  X = rmvnorm(n, rep(0,p-1), sigma)
  X=cbind(X[,1:7], matrix(0.5*runif(n)+0.5,ncol=1), X[,8:(p-1)])
  beta = c(3,-1.5,0,0,2, rep(0,p-5)); beta = matrix(beta, ncol=1)
  Y = X %*% beta +  2*X[,8]*matrix(rnorm(n),ncol=1)    
  
  out<-QaSIS(y=Y,x=X,tau=tau)
  return(rank(-out)[active]) #rank of active variables
}

example1_BIC = function(repnumb,n,p,sig=0.5,tau){
  sigma = array(0,c(p-1,p-1))
  for (ii in 1:(p-1)){
    for (jj in 1:(p-1)){
      sigma[ii,jj] = sig^(abs(ii-jj))
    }
  }
  set.seed(repnumb)
  X = rmvnorm(n, rep(0,p-1), sigma)
  X=cbind(X[,1:7], matrix(0.5*runif(n)+0.5,ncol=1), X[,8:(p-1)])
  beta = c(3,-1.5,0,0,2, rep(0,p-5)); beta = matrix(beta, ncol=1)
  Y = X %*% beta +  2*X[,8]*matrix(rnorm(n),ncol=1)     
  
  l_width<-0.05
  l_max<-0.35
  result = selection_bic(Y,X,l_max,l_width,Cn=log(p))
  return(result) 
}

example1_seq = function(repnumb,n,p,sig=0.5,tau){
  sigma = array(0,c(p-1,p-1))
  for (ii in 1:(p-1)){
    for (jj in 1:(p-1)){
      sigma[ii,jj] = sig^(abs(ii-jj))
    }
  }
  set.seed(repnumb)
  X = rmvnorm(n, rep(0,p-1), sigma)
  X=cbind(X[,1:7], matrix(0.5*runif(n)+0.5,ncol=1), X[,8:(p-1)])
  beta = c(3,-1.5,0,0,2, rep(0,p-5)); beta = matrix(beta, ncol=1)
  Y = X %*% beta +  2*X[,8]*matrix(rnorm(n),ncol=1)    
  result = STEPWISE(X, Y, eta1 = 0, eta2 = 1, initvars = NULL, tau=tau)
  return(result) 
}




## Example: survival simulation setting

example1_QA_surv = function(repnumb,n,p,sig=0.5,tau){
  sigma = array(0,c(p-1,p-1))
  for (ii in 1:(p-1)){
    for (jj in 1:(p-1)){
      sigma[ii,jj] = sig^(abs(ii-jj))
    }
  }
  set.seed(repnumb)
  X = rmvnorm(n, rep(0,p-1), sigma)
  X=cbind(X[,1:7], matrix(0.5*runif(n)+0.5,ncol=1), X[,8:(p-1)])
  beta = c(3,-1.5,0,0,2, rep(0,p-5)); beta = matrix(beta, ncol=1)
  Y = X %*% beta +  2*X[,8]*matrix(rnorm(n),ncol=1)    
  ind_mix = rmultinom(n,1,c(0.5,0.25,0.25))
  T = matrix(ind_mix[1,]*rnorm(n,-5,4) + ind_mix[2,]*rnorm(n,5,1) + ind_mix[3,]*rnorm(n,25,1), c(n,1))
  T = matrix(rep(100,n), ncol=1)
  delta = (Y <= T)+0
  YY=apply(cbind(T,Y),1,min)
  out = QaSIS.surv(X,YY,delta,tau)
  active=c(1,2,3,4,5,6,7,8)
  return(rank(-out)[active]) #rank of active variables
}


example1_BIC_surv = function(repnumb,n,p,sig=0.5,tau){
  sigma = array(0,c(p-1,p-1))
  for (ii in 1:(p-1)){
    for (jj in 1:(p-1)){
      sigma[ii,jj] = sig^(abs(ii-jj))
    }
  }
  set.seed(repnumb)
  X = rmvnorm(n, rep(0,p-1), sigma)
  X=cbind(X[,1:7], matrix(0.5*runif(n)+0.5,ncol=1), X[,8:(p-1)])
  beta = c(3,-1.5,0,0,2, rep(0,p-5)); beta = matrix(beta, ncol=1)
  Y = X %*% beta +  2*X[,8]*matrix(rnorm(n),ncol=1)    
  ind_mix = rmultinom(n,1,c(0.5,0.25,0.25))
  T = matrix(ind_mix[1,]*rnorm(n,-5,4) + ind_mix[2,]*rnorm(n,5,1) + ind_mix[3,]*rnorm(n,25,1), c(n,1))
  T = matrix(rep(100,n), ncol=1)
  delta = (Y <= T)+0
  YY=apply(cbind(T,Y),1,min)
  l_width<-0.05
  l_max<-0.35
  result = selection_bic_surv(YY,X,l_max,l_width,Cn=log(p),delta,T,tau)
  return(result)  
}


example1_seq_surv = function(repnumb,n,p,sig=0.5,tau){
  sigma = array(0,c(p-1,p-1))
  for (ii in 1:(p-1)){
    for (jj in 1:(p-1)){
      sigma[ii,jj] = sig^(abs(ii-jj))
    }
  }
  set.seed(repnumb)
  X = rmvnorm(n, rep(0,p-1), sigma)
  X=cbind(X[,1:7], matrix(0.5*runif(n)+0.5,ncol=1), X[,8:(p-1)])
  beta = c(3,-1.5,0,0,2, rep(0,p-5)); beta = matrix(beta, ncol=1)
  Y = X %*% beta +  2*X[,8]*matrix(rnorm(n),ncol=1)   
  ind_mix = rmultinom(n,1,c(0.5,0.25,0.25))
  T = matrix(ind_mix[1,]*rnorm(n,-5,4) + ind_mix[2,]*rnorm(n,5,1) + ind_mix[3,]*rnorm(n,25,1), c(n,1))
  T = matrix(rep(100,n), ncol=1)
  delta = (Y <= T)+0
  YY=apply(cbind(T,Y),1,min)
  result = STEPWISE_surv(X, YY, eta1 = 0, eta2 = 1, initvars = NULL, tau=tau, delta, T)
  return(result)  
}



#Example 1 (tau=0.5)
NNN=100;
tic()
ex1_05_QA = mclapply(1:NNN, example1_QA_surv, n=300, p=1000, sig=0.5,tau=0.5, mc.cores = 4)
ex1_05_QA = matrix(unlist(ex1_05_QA), byrow=T, ncol=8)
toc_QA=toc()
tic()
ex1_05_BIC = mclapply(1:NNN, example1_BIC_surv, n=300, p=1000, sig=0.5, tau=0.5,mc.cores = 4)
toc_BIC=toc()
tic()
ex1_05_seq = mclapply(1:NNN, example1_seq_surv, n=300, p=1000, sig=0.5, tau=0.5,mc.cores = 5)
toc_seq=toc()


#Investigating model selection property
ex1_05_BIC_neat=NULL; ex1_05_seq_neat=NULL; fp_BIC_neat=NULL; fp_seq_neat=NULL; 
for (ii in 1:NNN){
  fp_BIC_neat = c(fp_BIC_neat, length(setdiff(ex1_05_BIC[[ii]]$sel_H, c(1,2,5))))
  fp_seq_neat = c(fp_seq_neat, length(setdiff(ex1_05_seq[[ii]]$sel_H, c(1,2,5))))
  ex1_05_BIC_n = NULL; ex1_05_seq_n=NULL
  for (kk in 1:8){
    ex1_05_BIC_n=c(ex1_05_BIC_n,is.element(kk,ex1_05_BIC[[ii]]$sel_H)+0)
    ex1_05_seq_n=c(ex1_05_seq_n,is.element(kk,ex1_05_seq[[ii]]$finalset)+0)
  }
  ex1_05_BIC_neat=rbind(ex1_05_BIC_neat,ex1_05_BIC_n)
  ex1_05_seq_neat=rbind(ex1_05_seq_neat,ex1_05_seq_n)
}




##Example 1 (tau=0.3)
NNN=100;
tic()
ex1_05_QA = mclapply(1:NNN, example1_QA_surv, n=300, p=1000, sig=0.5,tau=0.3, mc.cores = 4)
ex1_05_QA = matrix(unlist(ex1_05_QA), byrow=T, ncol=8)
toc_QA=toc()
tic()
ex1_05_BIC = mclapply(1:NNN, example1_BIC_surv, n=300, p=1000, sig=0.5, tau=0.3,mc.cores = 4)
toc_BIC=toc()
tic()
ex1_05_seq = mclapply(1:NNN, example1_seq_surv, n=300, p=1000, sig=0.5, tau=0.3,mc.cores = 5)
toc_seq=toc()


#Investigating model selection property
ex1_05_BIC_neat=NULL; ex1_05_seq_neat=NULL; fp_BIC_neat=NULL; fp_seq_neat=NULL; 
for (ii in 1:NNN){
  fp_BIC_neat = c(fp_BIC_neat, length(setdiff(ex1_05_BIC[[ii]]$sel_H, c(1,2,5))))
  fp_seq_neat = c(fp_seq_neat, length(setdiff(ex1_05_seq[[ii]]$sel_H, c(1,2,5))))
  ex1_05_BIC_n = NULL; ex1_05_seq_n=NULL
  for (kk in 1:8){
    ex1_05_BIC_n=c(ex1_05_BIC_n,is.element(kk,ex1_05_BIC[[ii]]$sel_H)+0)
    ex1_05_seq_n=c(ex1_05_seq_n,is.element(kk,ex1_05_seq[[ii]]$finalset)+0)
  }
  ex1_05_BIC_neat=rbind(ex1_05_BIC_neat,ex1_05_BIC_n)
  ex1_05_seq_neat=rbind(ex1_05_seq_neat,ex1_05_seq_n)
}




##Example 1 (tau=0.7)
NNN=100;
tic()
ex1_05_QA = mclapply(1:NNN, example1_QA_surv, n=300, p=1000, sig=0.5,tau=0.7, mc.cores = 4)
ex1_05_QA = matrix(unlist(ex1_05_QA), byrow=T, ncol=8)
toc_QA=toc()
tic()
ex1_05_BIC = mclapply(1:NNN, example1_BIC_surv, n=300, p=1000, sig=0.5, tau=0.7,mc.cores = 4)
toc_BIC=toc()
tic()
ex1_05_seq = mclapply(1:NNN, example1_seq_surv, n=300, p=1000, sig=0.5, tau=0.7,mc.cores = 5)
toc_seq=toc()


#Investigating model selection property
ex1_05_BIC_neat=NULL; ex1_05_seq_neat=NULL; fp_BIC_neat=NULL; fp_seq_neat=NULL; 
for (ii in 1:NNN){
  fp_BIC_neat = c(fp_BIC_neat, length(setdiff(ex1_05_BIC[[ii]]$sel_H, c(1,2,5))))
  fp_seq_neat = c(fp_seq_neat, length(setdiff(ex1_05_seq[[ii]]$sel_H, c(1,2,5))))
  ex1_05_BIC_n = NULL; ex1_05_seq_n=NULL
  for (kk in 1:8){
    ex1_05_BIC_n=c(ex1_05_BIC_n,is.element(kk,ex1_05_BIC[[ii]]$sel_H)+0)
    ex1_05_seq_n=c(ex1_05_seq_n,is.element(kk,ex1_05_seq[[ii]]$finalset)+0)
  }
  ex1_05_BIC_neat=rbind(ex1_05_BIC_neat,ex1_05_BIC_n)
  ex1_05_seq_neat=rbind(ex1_05_seq_neat,ex1_05_seq_n)
}




