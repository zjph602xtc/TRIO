###################################################################################################
###### DIrect Sparse regression procedure using COvariance from Multi-modality data (DISCOM) ######
###################################################################################################
######                                      Guan Yu                                          ######
######                                  Date: 12/30/2018                                     ######
###################################################################################################
# R codes for the following methods:
# 1.Lasso for the complete data
# 2.Lasso for the imputed data
# 3.Ridge regression for the complete data
# 4.Ridge regression for the imputed data
# 5.iMSF
# 6.DISCOM
# 7.DISCOM-Huber
# 8.Fast-DISCOM
# 9.Fast-DISCOM-Huber 
###################################################################################################


###################################################################################################
# The following functions are used to compute the initial covariance matrices 
###################################################################################################
compute.xtx<-function(x,robust=0,k_value=1) 
  # robust=1 for huber estimate, k_value is used in huber function
{
  p=ncol(x)
  cov.matrix=matrix(NA,p,p)
  if(robust==0){cov.matrix=cov(x,use="pairwise.complete.obs")}
  else{
    for(i in 1:p)
    {
      for(j in 1:p){
        index=which((!is.na(x[,i]))&(!is.na(x[,j])))
        x1=x[index,i]
        x2=x[index,j]
        cov.matrix[i,j]=huberM((x1-mean(x1))*(x2-mean(x2)),k=k_value*sqrt(length(index)/log(p)))$mu
      }
    }
  }
  cov.matrix  
}

compute.xty<-function(x,y,robust=0,k_value=1)
  # robust=1 for huber estimate, k_value is used in huber function
{
  p=ncol(x)
  cov.vector=rep(NA,p)
  if(robust==0){cov.vector=cov(x,y,use="pairwise.complete.obs")}
  else{
    for(i in 1:p){
      index=which((!is.na(x[,i])))
      x1=x[index,i]
      x2=y[index]
      cov.vector[i]=huberM((x1-mean(x1))*(x2-mean(x2)),k=k_value*sqrt(length(index)/log(p)))$mu
    }
  }
  cov.vector
}

###################################################################################################
# Load packages
###################################################################################################
library(MASS)
library(glmnet)
library(scout)
library(Matrix)
library(softImpute)
library(grplasso)
library(robustbase)

###################################################################################################
# Set parameters in a simulation example (Example 1 in the paper) with three modalities
###################################################################################################
nsim=30             # number of simulations
n1=100              # number of samples with complete observations
n2=100              # number of samples with observations from modality 1 and modality 2
n3=100              # number of samples with observations from modality 1 and modality 3
n4=100              # number of samples with observations only from modality 1 

n.tun=200           # number of tuning samples
n.test=400          # number of test samples
p1=100              # number of predictors from modality 1
p2=100              # number of predictors from modality 2
p3=100              # number of predictors from modality 3
p=p1+p2+p3          # the total number of predictors

###################################################################################################
# Set the true regression coefficient vector in the linear model
###################################################################################################
beta.true=c(rep(0.5,3),rep(0,p1-3),rep(0.5,3),rep(0,p2-3),rep(0.5,3),rep(0,p3-3))

###################################################################################################
# Set the true covariance matrix of the predictor vector X
###################################################################################################
cov.true=matrix(NA,p,p)   
for(i in 1:p){for(j in 1:p){cov.true[i,j]=0.6^(abs(i-j))}} 

###################################################################################################
# Use the following five measures to compare different methods
###################################################################################################
est_results=matrix(NA,nsim,9)      # estimation error
mse_results=matrix(NA,nsim,9)      # mean squared error 
fpr_results=matrix(NA,nsim,9)      # false positive rate
fnr_results=matrix(NA,nsim,9)      # false negative rate
time_results=matrix(NA,nsim,10)    # the elapsed time in seconds using R

for(sim in 1:nsim){
  ###################################################################################################
  # Generate the training data, imputed training data, tuning data, and testing data
  ###################################################################################################
  set.seed(sim)
  X.all=mvrnorm(n1+n2+n3+n4,rep(0,p),cov.true)
  Y=X.all%*%beta.true+rnorm(n1+n2+n3+n4)
  X.all[(n1+1):(n1+n2),(p1+p2+1):(p1+p2+p3)]=NA
  X.all[(n1+n2+1):(n1+n2+n3),(p1+1):(p1+p2)]=NA
  X.all[(n1+n2+n3+1):(n1+n2+n3+n4),(p1+1):(p1+p2+p3)]=NA  # X.all and Y will be the training data
  
  start=proc.time()
  fit1=softImpute(X.all,rank=50,lambda=30)
  X.imputed=complete(X.all,fit1)    # Use the Soft-thresholded SVD method to impute the missing data
  time=proc.time()-start
  time_results[sim,1]=as.vector(time[3])
  
  X.tun=mvrnorm(n.tun,rep(0,p),cov.true)
  Y.tun=X.tun%*%beta.true+rnorm(n.tun)                    # X.tun and Y.tun will be the tuning data
  
  X.test=mvrnorm(n.test,rep(0,p),cov.true)
  Y.test=X.test%*%beta.true+rnorm(n.test)                 # X.test and Y.test will be the testing data
  
  ###################################################################################################
  # Lasso method for complete data 
  ###################################################################################################
  {
    start=proc.time()
    nlambda=30         
    lasso=glmnet(X.all[1:n1,],Y[1:n1],nlambda=nlambda,intercept=F,alpha=1)
    lambda.all=as.vector(lasso$lambda)
    nlambda=length(lambda.all)
    lasso.tun.error=rep(NA,nlambda)
    for(i in 1:nlambda)
    {
      beta.lasso=as.vector(lasso$beta[,i])
      lasso.tun.values=as.vector(X.tun%*%beta.lasso)
      lasso.tun.error[i]=mean((Y.tun-lasso.tun.values)^2)
    }
    opt.lambda=lambda.all[which.min(lasso.tun.error)]
    beta.lasso=as.vector(glmnet(X.all[1:n1,],Y[1:n1],lambda=opt.lambda,intercept=F,alpha=1)$beta)
    lasso.test.values=as.vector(X.test%*%beta.lasso)
    lasso.test.error=mean((Y.test-lasso.test.values)^2)
    lasso.fpr=sum((beta.true==0)&(beta.lasso!=0))/sum(beta.true==0)
    lasso.fnr=sum((beta.true!=0)&(beta.lasso==0))/sum(beta.true!=0)
    lasso.est.error=sqrt(sum((beta.lasso-beta.true)^2))
    time=proc.time()-start
    time_results[sim,2]=as.vector(time[3])
  } 
  
  ###################################################################################################
  # Lasso method for the imputed data 
  ###################################################################################################
  {
    start=proc.time()
    nlambda=30         
    lasso=glmnet(X.imputed,Y,nlambda=nlambda,intercept=F,alpha=1)
    lambda.all=as.vector(lasso$lambda)
    nlambda=length(lambda.all)
    lasso.tun.error=rep(NA,nlambda)
    for(i in 1:nlambda)
    {
      beta.lasso=as.vector(lasso$beta[,i])
      lasso.tun.values=as.vector(X.tun%*%beta.lasso)
      lasso.tun.error[i]=mean((Y.tun-lasso.tun.values)^2)
    }
    opt.lambda=lambda.all[which.min(lasso.tun.error)]
    beta.lasso=as.vector(glmnet(X.imputed,Y,lambda=opt.lambda,intercept=F,alpha=1)$beta)
    lasso.test.values=as.vector(X.test%*%beta.lasso)
    lasso.imputed.test.error=mean((Y.test-lasso.test.values)^2)
    lasso.imputed.fpr=sum((beta.true==0)&(beta.lasso!=0))/sum(beta.true==0)
    lasso.imputed.fnr=sum((beta.true!=0)&(beta.lasso==0))/sum(beta.true!=0)
    lasso.imputed.est.error=sqrt(sum((beta.lasso-beta.true)^2))
    time=proc.time()-start
    time_results[sim,3]=as.vector(time[3])+time_results[sim,1]
  } 
  
  ###################################################################################################
  # Ridge regression for complete data
  ###################################################################################################
  {
    start=proc.time()
    nlambda=30        
    ridge=glmnet(X.all[1:n1,],Y[1:n1],nlambda=nlambda,intercept=F,alpha=0)
    lambda.all=as.vector(ridge$lambda)
    nlambda=length(lambda.all)
    ridge.tun.error=rep(NA,nlambda)
    for(i in 1:nlambda)
    {
      beta.ridge=as.vector(ridge$beta[,i])
      ridge.tun.values=as.vector(X.tun%*%beta.ridge)
      ridge.tun.error[i]=mean((Y.tun-ridge.tun.values)^2)
    }
    opt.lambda=lambda.all[which.min(ridge.tun.error)]
    beta.ridge=as.vector(glmnet(X.all[1:n1,],Y[1:n1],lambda=opt.lambda,intercept=F,alpha=0)$beta)
    ridge.test.values=as.vector(X.test%*%beta.ridge)
    ridge.test.error=mean((Y.test-ridge.test.values)^2)
    ridge.fpr=sum((beta.true==0)&(beta.ridge!=0))/sum(beta.true==0)
    ridge.fnr=sum((beta.true!=0)&(beta.ridge==0))/sum(beta.true!=0)
    ridge.est.error=sqrt(sum((beta.ridge-beta.true)^2))
    time=proc.time()-start
    time_results[sim,4]=as.vector(time[3])
  } 
  
  ###################################################################################################
  # Ridge regression for the imputed data
  ###################################################################################################
  {
    start=proc.time()
    nlambda=30         
    ridge=glmnet(X.imputed,Y,nlambda=nlambda,intercept=F,alpha=0)
    lambda.all=as.vector(ridge$lambda)
    nlambda=length(lambda.all)
    ridge.tun.error=rep(NA,nlambda)
    for(i in 1:nlambda)
    {
      beta.ridge=as.vector(ridge$beta[,i])
      ridge.tun.values=as.vector(X.tun%*%beta.ridge)
      ridge.tun.error[i]=mean((Y.tun-ridge.tun.values)^2)
    }
    opt.lambda=lambda.all[which.min(ridge.tun.error)]
    beta.ridge=as.vector(glmnet(X.imputed,Y,lambda=opt.lambda,intercept=F,alpha=0)$beta)
    ridge.test.values=as.vector(X.test%*%beta.ridge)
    ridge.imputed.test.error=mean((Y.test-ridge.test.values)^2)
    ridge.imputed.fpr=sum((beta.true==0)&(beta.ridge!=0))/sum(beta.true==0)
    ridge.imputed.fnr=sum((beta.true!=0)&(beta.ridge==0))/sum(beta.true!=0)
    ridge.imputed.est.error=sqrt(sum((beta.ridge-beta.true)^2))
    time=proc.time()-start
    time_results[sim,5]=as.vector(time[3])+time_results[sim,1]
  }  
  
  ###################################################################################################
  # iMSF method
  ###################################################################################################
  {
    start=proc.time()
    Y.train.new=c(Y[1:n1]/sqrt(n1),Y[(1+n1):(n1+n2)]/sqrt(n2),Y[(1+n1+n2):(n1+n2+n3)]/sqrt(n3),
                  Y[(1+n1+n2+n3):(n1+n2+n3+n4)]/sqrt(n4))
    X.train.new=as.matrix(bdiag(X.all[1:n1,]/sqrt(n1),X.all[(1+n1):(n1+n2),1:(p1+p2)]/sqrt(n2),
                                X.all[(1+n1+n2):(n1+n2+n3),c(1:p1,(p1+p2+1):(p1+p2+p3))]/sqrt(n3),
                                X.all[(1+n1+n2+n3):(n1+n2+n3+n4),1:p1]/sqrt(n4)))
    group.index=c(1:(p1+p2+p3),1:(p1+p2),1:p1,(p1+p2+1):(p1+p2+p3),1:p1)
    nlambda=30
    lambda.all=seq(1,0.05,length=nlambda)
    imsf.tun.error=rep(NA,nlambda)
    
    for(i in 1:nlambda){
      imsf=grplasso(X.train.new,Y.train.new,group.index,lambda=lambda.all[i],center=F,
                    standardize=F,model=LinReg(),
                    control = grpl.control(trace = 0))
      beta.imsf=as.vector(imsf$coefficients[1:(p1+p2+p3)])
      imsf.tun.values=as.vector(X.tun%*%beta.imsf)
      imsf.tun.error[i]=mean((Y.tun-imsf.tun.values)^2)
    }
    opt.lambda=lambda.all[which.min(imsf.tun.error)]
    imsf.fit=grplasso(X.train.new,Y.train.new,group.index,lambda=opt.lambda,center=F,
                      standardize=F,model=LinReg(),
                      control = grpl.control(trace = 0))
    beta.imsf=as.vector(as.vector(imsf.fit$coefficients[1:(p1+p2+p3)]))
    imsf.test.values=as.vector(X.test%*%beta.imsf)
    imsf.test.error=mean((Y.test-imsf.test.values)^2)
    imsf.fpr=sum((beta.true==0)&(beta.imsf!=0))/sum(beta.true==0)
    imsf.fnr=sum((beta.true!=0)&(beta.imsf==0))/sum(beta.true!=0)
    imsf.est.error=sqrt(sum((beta.imsf-beta.true)^2))
    time=proc.time()-start
    time_results[sim,6]=as.vector(time[3])
  }
  
  ###################################################################################################
  # DISCOM 
  ###################################################################################################
  {
    start=proc.time()
    nalpha=10
    nlambda=30    
    alpha.all=10^seq(0,-3,length=nalpha)
    lambda.all=10^(seq(log10(2),-3,length=nlambda))
    
    DISCOM.tun.error=array(NA,dim=c(nalpha,nalpha,nlambda))
    
    X.all <- read.csv('sig_p_0.40_cor_0.30_missing_0.40.csv', header = FALSE)
    Y <- read.csv('sig_p_0.40_cor_0.30_missing_0.40_y.csv', header = FALSE)
    
    xtx.raw=compute.xtx(X.all)
    xty=compute.xty(X.all,Y)
    xtx.raw.I=as.matrix(bdiag(xtx.raw[1:p1,1:p1],xtx.raw[(p1+1):(p1+p2),(p1+1):(p1+p2)],
                              xtx.raw[(p1+p2+1):(p1+p2+p3),(p1+p2+1):(p1+p2+p3)]))
    xtx.raw.C=xtx.raw-xtx.raw.I
    shrink.target=sum(diag(xtx.raw))/p
    
    for(i in 1:nalpha){  
      for(j in 1:nalpha){
        alpha1=alpha.all[i]
        alpha2=alpha.all[j]
        xtx=alpha1*xtx.raw.I+alpha2*xtx.raw.C+(1-alpha1)*shrink.target*diag(p)
        if(min(eigen(xtx)$values)<0){DISCOM.tun.error[i,j,]=10^8}
        else{
          beta.initial=rep(0,p)
          for(k in 1:nlambda){
            beta.cov=as.vector(crossProdLasso(xtx,xty,lambda.all[k],beta.init=beta.initial)$beta)
            betal.initial=beta.cov
            DISCOM.tun.values=as.vector(as.matrix(X.tun)%*%beta.cov)
            DISCOM.tun.error[i,j,k]=mean((Y.tun-DISCOM.tun.values)^2)
          }
        }
      }
    }
    
    opt.index=as.vector(which(DISCOM.tun.error==min(DISCOM.tun.error),arr.ind=T)[1,])
    opt.alpha1=alpha.all[opt.index[1]]
    opt.alpha2=alpha.all[opt.index[2]]
    opt.lambda=lambda.all[opt.index[3]]
    xtx=opt.alpha1*xtx.raw.I+opt.alpha2*xtx.raw.C+(1-opt.alpha1)*shrink.target*diag(p)
    beta.cov=as.vector(crossProdLasso(xtx,xty,opt.lambda)$beta)
    predict.test.values=as.vector(X.test%*%beta.cov)
    
    DISCOM.test.error=mean((Y.test-predict.test.values)^2)
    DISCOM.fpr=sum((beta.true==0)&(beta.cov!=0))/sum(beta.true==0)
    DISCOM.fnr=sum((beta.true!=0)&(beta.cov==0))/sum(beta.true!=0)
    DISCOM.est.error=sqrt(sum((beta.cov-beta.true)^2))
    time=proc.time()-start
    time_results[sim,7]=as.vector(time[3])
  }
  
  ###################################################################################################
  # DISCOM-Huber
  ###################################################################################################
  {
    start=proc.time()
    nalpha=10
    nlambda=30    
    alpha.all=10^seq(0,-3,length=nalpha)
    lambda.all=10^(seq(log10(2),-3,length=nlambda))
    
    DISCOM.Huber.tun.error=array(NA,dim=c(nalpha,nalpha,nlambda))
    
    xtx.raw=compute.xtx(X.all,1,1)
    xty=compute.xty(X.all,Y,1,1) 
    xtx.raw.I=as.matrix(bdiag(xtx.raw[1:p1,1:p1],xtx.raw[(p1+1):(p1+p2),(p1+1):(p1+p2)],
                              xtx.raw[(p1+p2+1):(p1+p2+p3),(p1+p2+1):(p1+p2+p3)]))
    xtx.raw.C=xtx.raw-xtx.raw.I
    shrink.target=sum(diag(xtx.raw))/p
    
    for(i in 1:nalpha){  
      for(j in 1:nalpha){
        alpha1=alpha.all[i]
        alpha2=alpha.all[j]
        xtx=alpha1*xtx.raw.I+alpha2*xtx.raw.C+(1-alpha1)*shrink.target*diag(p)
        if(min(eigen(xtx)$values)<0){DISCOM.Huber.tun.error[i,j,]=10^8}
        else{
          beta.initial=rep(0,p)
          for(k in 1:nlambda){
            beta.cov=as.vector(crossProdLasso(xtx,xty,lambda.all[k],beta.init=beta.initial)$beta)
            betal.initial=beta.cov
            DISCOM.Huber.tun.values=as.vector(as.matrix(X.tun)%*%beta.cov)
            DISCOM.Huber.tun.error[i,j,k]=mean((Y.tun-DISCOM.Huber.tun.values)^2)
          }
        }
      }
    }
    
    opt.index=as.vector(which(DISCOM.Huber.tun.error==min(DISCOM.Huber.tun.error),arr.ind=T)[1,])
    opt.alpha1=alpha.all[opt.index[1]]
    opt.alpha2=alpha.all[opt.index[2]]
    opt.lambda=lambda.all[opt.index[3]]
    xtx=opt.alpha1*xtx.raw.I+opt.alpha2*xtx.raw.C+(1-opt.alpha1)*shrink.target*diag(p)
    beta.cov=as.vector(crossProdLasso(xtx,xty,opt.lambda)$beta)
    predict.test.values=as.vector(X.test%*%beta.cov)
    
    DISCOM.Huber.test.error=mean((Y.test-predict.test.values)^2)
    DISCOM.Huber.fpr=sum((beta.true==0)&(beta.cov!=0))/sum(beta.true==0)
    DISCOM.Huber.fnr=sum((beta.true!=0)&(beta.cov==0))/sum(beta.true!=0)
    DISCOM.Huber.est.error=sqrt(sum((beta.cov-beta.true)^2))
    time=proc.time()-start
    time_results[sim,8]=as.vector(time[3])
  }
  
  ###################################################################################################
  # Fast-DISCOM
  ###################################################################################################
  {
    start=proc.time()
    
    xtx.raw=compute.xtx(X.all)
    xty=compute.xty(X.all,Y)
    xtx.raw.I=as.matrix(bdiag(xtx.raw[1:p1,1:p1],xtx.raw[(p1+1):(p1+p2),(p1+1):(p1+p2)],
                              xtx.raw[(p1+p2+1):(p1+p2+p3),(p1+p2+1):(p1+p2+p3)]))
    xtx.raw.C=xtx.raw-xtx.raw.I
    
    shrink.target=sum(diag(xtx.raw))/p
    c1=sqrt(log(p)/200)
    c2=sqrt(log(p)/100)
    A.matrix=(c2-c1)*xtx.raw.I+c1*shrink.target*diag(p)
    K.min=max(0,min(eigen(xtx.raw)$values)/(c2*min(eigen(xtx.raw)$values)-min(eigen(A.matrix)$values)))
    K.max=1/c2
    n.K=30  
    nlambda=30
    K.all=seq(K.min,K.max,length=n.K)
    lambda.all=10^(seq(log10(2),-3,length=nlambda))
    DISCOM.fast.tun.error=matrix(NA,n.K,nlambda)
    
    for(i in 1:n.K){
      alpha1=1-K.all[i]*c1
      alpha2=1-K.all[i]*c2
      xtx=alpha1*xtx.raw.I+alpha2*xtx.raw.C+(1-alpha1)*shrink.target*diag(p)
      beta.initial=rep(0,p)
      for(k in 1:nlambda){
        beta.cov=as.vector(crossProdLasso(xtx,xty,lambda.all[k],beta.init=beta.initial)$beta)
        betal.initial=beta.cov
        DISCOM.fast.tun.values=as.vector(as.matrix(X.tun)%*%beta.cov)
        DISCOM.fast.tun.error[i,k]=mean((Y.tun-DISCOM.fast.tun.values)^2)
      }
    }
    
    opt.index=as.vector(which(DISCOM.fast.tun.error==min(DISCOM.fast.tun.error),arr.ind=T)[1,])
    opt.alpha1=1-K.all[opt.index[1]]*c1
    opt.alpha2=1-K.all[opt.index[1]]*c2
    opt.lambda=lambda.all[opt.index[2]]
    xtx=opt.alpha1*xtx.raw.I+opt.alpha2*xtx.raw.C+(1-opt.alpha1)*shrink.target*diag(p)
    beta.cov=as.vector(crossProdLasso(xtx,xty,opt.lambda)$beta)
    predict.test.values=as.vector(X.test%*%beta.cov)
    
    DISCOM.fast.test.error=mean((Y.test-predict.test.values)^2)
    DISCOM.fast.fpr=sum((beta.true==0)&(beta.cov!=0))/sum(beta.true==0)
    DISCOM.fast.fnr=sum((beta.true!=0)&(beta.cov==0))/sum(beta.true!=0)
    DISCOM.fast.est.error=sqrt(sum((beta.cov-beta.true)^2))
    time=proc.time()-start
    time_results[sim,9]=as.vector(time[3])
  }
  
  ###################################################################################################
  # Fast-DISCOM-Huber
  ###################################################################################################
  {
    start=proc.time()
    
    xtx.raw=compute.xtx(X.all,1,1)
    xty=compute.xty(X.all,Y,1,1) 
    xtx.raw.I=as.matrix(bdiag(xtx.raw[1:p1,1:p1],xtx.raw[(p1+1):(p1+p2),(p1+1):(p1+p2)],
                              xtx.raw[(p1+p2+1):(p1+p2+p3),(p1+p2+1):(p1+p2+p3)]))
    xtx.raw.C=xtx.raw-xtx.raw.I
    
    shrink.target=sum(diag(xtx.raw))/p
    c1=sqrt(log(p)/200)
    c2=sqrt(log(p)/100)
    A.matrix=(c2-c1)*xtx.raw.I+c1*shrink.target*diag(p)
    K.min=max(0,min(eigen(xtx.raw)$values)/(c2*min(eigen(xtx.raw)$values)-min(eigen(A.matrix)$values)))
    K.max=1/c2
    n.K=30  
    nlambda=30
    K.all=seq(K.min,K.max,length=n.K)
    lambda.all=10^(seq(log10(2),-3,length=nlambda))
    
    DISCOM.Huber.fast.tun.error=matrix(NA,n.K,nlambda)
    
    for(i in 1:n.K){
      alpha1=1-K.all[i]*c1
      alpha2=1-K.all[i]*c2
      xtx=alpha1*xtx.raw.I+alpha2*xtx.raw.C+(1-alpha1)*shrink.target*diag(p)
      beta.initial=rep(0,p)
      for(k in 1:nlambda){
        beta.cov=as.vector(crossProdLasso(xtx,xty,lambda.all[k],beta.init=beta.initial)$beta)
        betal.initial=beta.cov
        DISCOM.Huber.fast.tun.values=as.vector(as.matrix(X.tun)%*%beta.cov)
        DISCOM.Huber.fast.tun.error[i,k]=mean((Y.tun-DISCOM.Huber.fast.tun.values)^2)
      }
    }
    
    opt.index=as.vector(which(DISCOM.Huber.fast.tun.error==min(DISCOM.Huber.fast.tun.error),
                              arr.ind=T)[1,])
    opt.alpha1=1-K.all[opt.index[1]]*c1
    opt.alpha2=1-K.all[opt.index[1]]*c2
    opt.lambda=lambda.all[opt.index[2]]
    xtx=opt.alpha1*xtx.raw.I+opt.alpha2*xtx.raw.C+(1-opt.alpha1)*shrink.target*diag(p)
    beta.cov=as.vector(crossProdLasso(xtx,xty,opt.lambda)$beta)
    predict.test.values=as.vector(X.test%*%beta.cov)
    
    DISCOM.Huber.fast.test.error=mean((Y.test-predict.test.values)^2)
    DISCOM.Huber.fast.fpr=sum((beta.true==0)&(beta.cov!=0))/sum(beta.true==0)
    DISCOM.Huber.fast.fnr=sum((beta.true!=0)&(beta.cov==0))/sum(beta.true!=0)
    DISCOM.Huber.fast.est.error=sqrt(sum((beta.cov-beta.true)^2))
    time=proc.time()-start
    time_results[sim,10]=as.vector(time[3])
  }
  
  mse_results[sim,]=c(lasso.test.error,lasso.imputed.test.error,ridge.test.error,
                      ridge.imputed.test.error,
                      imsf.test.error,DISCOM.test.error,DISCOM.Huber.test.error,
                      DISCOM.fast.test.error,DISCOM.Huber.fast.test.error)
  
  fpr_results[sim,]=c(lasso.fpr,lasso.imputed.fpr,ridge.fpr,ridge.imputed.fpr,
                      imsf.fpr,DISCOM.fpr,DISCOM.Huber.fpr,DISCOM.fast.fpr,DISCOM.Huber.fast.fpr)
  
  fnr_results[sim,]=c(lasso.fnr,lasso.imputed.fnr,ridge.fnr,ridge.imputed.fnr,
                      imsf.fnr,DISCOM.fnr,DISCOM.Huber.fnr,DISCOM.fast.fnr,DISCOM.Huber.fast.fnr)
  
  est_results[sim,]=c(lasso.est.error,lasso.imputed.est.error,ridge.est.error,ridge.imputed.est.error,
                      imsf.est.error,DISCOM.est.error,DISCOM.Huber.est.error,
                      DISCOM.fast.est.error,DISCOM.Huber.fast.est.error)
}

###################################################################################################
# Report the results
###################################################################################################
mse.mean=apply(mse_results,2,mean)
mse.sd=apply(mse_results,2,sd)/sqrt(nsim)
fpr.mean=apply(fpr_results,2,mean)
fpr.sd=apply(fpr_results,2,sd)/sqrt(nsim)
fnr.mean=apply(fnr_results,2,mean)
fnr.sd=apply(fnr_results,2,sd)/sqrt(nsim)
est.mean=apply(est_results,2,mean)
est.sd=apply(est_results,2,sd)/sqrt(nsim)
time.mean=apply(time_results[,2:10],2,mean)
time.sd=apply(time_results[,2:10],2,sd)/sqrt(nsim)

finalresult=data.frame(Method=c("Lasso","ImputeLasso","Ridge","ImputeRidge","Imsf","DISCOM",
                                "DISCOM_Huber","Fast_DISCOM","Fast_DISCOM_Huber"),
                       est.mean,mse.mean,fpr.mean,fnr.mean,time.mean) 
finalresult
