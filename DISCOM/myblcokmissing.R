library(MASS)
library(glmnet)
library(scout)
library(Matrix)
library(softImpute)
library(grplasso)
library(robustbase)
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

cleanx <- function(x, nsub){
  for (i in 1:(nrow(x)/nsub)){
    if (i==1){
      X <- x[((i-1)*nsub+1):(i*nsub),]
    }else{
      X <- cbind(X, x[((i-1)*nsub+1):(i*nsub),])
    }
  }
  return(X)
}


# DISCOM.test.error < - NULL
# DISCOM.test.std <- NULL

err <- NULL
estd <- NULL
allres <- NULL
p <- 100
for (r in c(1,3,5,7,10)){
  for (smooth_p in c(1)){
    for (missing_p in c(0)){
      print(smooth_p)
      # sig_prop <- 1
      # smooth_p <- 1
      # missing_p <- 0
      X.test <- read.csv(paste0('G:/tensor_data_rank/',sprintf('t_sig_p_%.2f_Vcor_%d_missing_%.2f_r_%d_%d.csv',sig_prop, smooth_p, missing_p,r, 1)), header = F)
      Y.test <- read.csv(paste0('G:/tensor_data_rank/',sprintf('t_sig_p_%.2f_Vcor_%d_missing_%.2f_y_r_%d_%d.csv',sig_prop, smooth_p, missing_p,r, 1)), header = F)
      X.test <- cleanx(X.test, NROW(Y.test))
      
      myerror <- NULL
      for (iter in 1:100){
        # print(iter)
        x <- read.csv(paste0('G:/tensor_data_rank/',sprintf('sig_p_%.2f_Vcor_%d_missing_%.2f_r_%d_%d.csv',sig_prop, smooth_p, missing_p,r, iter)), header = F)
        y <- read.csv(paste0('G:/tensor_data_rank/',sprintf('sig_p_%.2f_Vcor_%d_missing_%.2f_y_r_%d_%d.csv',sig_prop, smooth_p, missing_p,r, iter)), header = F)
        x <- cleanx(x, 500)
        
        X.tun <- x
        for (j in 1:ncol(x)){
          X.tun[is.na(X.tun[,j]),j] <- mean(X.tun[,j], na.rm=T)
        }
        Y.tun <- y
        Y.tun <- as.matrix(Y.tun)
        
        
        nalpha=10
        nlambda=30    
        alpha.all=10^seq(0,-3,length=nalpha)
        lambda.all=10^(seq(log10(2),-3,length=nlambda))
        
        DISCOM.tun.error=array(NA,dim=c(nalpha,nalpha,nlambda))
        
        xtx.raw=compute.xtx(x)
        xty=compute.xty(x,y)
        
        blockxtx <- list()
        for (blk in 1:(nrow(xtx.raw)/p)){
          blockxtx[[blk]] <- xtx.raw[((blk-1)*p+1):(blk*p),((blk-1)*p+1):(blk*p)]
        }
        xtx.raw.I <- as.matrix(do.call(bdiag, blockxtx))
        xtx.raw.C <- xtx.raw-xtx.raw.I
        
        pp <- nrow(xtx.raw)
        shrink.target <- sum(diag(xtx.raw))/pp
        
        if (iter==1){
          for(i in 1:nalpha){
            for(j in 1:nalpha){
              alpha1=alpha.all[i]
              alpha2=alpha.all[j]
              xtx=alpha1*xtx.raw.I+alpha2*xtx.raw.C+(1-alpha1)*shrink.target*diag(pp)
              if(min(eigen(xtx)$values)<0){DISCOM.tun.error[i,j,]=10^8}
              else{
                beta.initial=rep(0,pp)
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
        }
        
        opt.alpha1=alpha.all[opt.index[1]]
        opt.alpha2=alpha.all[opt.index[2]]
        opt.lambda=lambda.all[opt.index[3]]
        xtx=opt.alpha1*xtx.raw.I+opt.alpha2*xtx.raw.C+(1-opt.alpha1)*shrink.target*diag(pp)
        beta.cov=as.vector(crossProdLasso(xtx,xty,opt.lambda)$beta)
        
        predict.test.values=as.vector(as.matrix(X.test)%*%beta.cov)
        
        myerror <- c(myerror, 
                     sqrt(sum((as.matrix(Y.test)-predict.test.values)^2))/sqrt(NROW(Y.test)))
        # DISCOM.fpr=sum((beta.true==0)&(beta.cov!=0))/sum(beta.true==0)
        # DISCOM.fnr=sum((beta.true!=0)&(beta.cov==0))/sum(beta.true!=0)
        # DISCOM.est.error=sqrt(sum((beta.cov-beta.true)^2))
        # time=proc.time()-start
        # time_results[sim,7]=as.vector(time[3])
        
      }
      err <- c(err, median(myerror))
      estd <- c(estd, IQR(myerror))
      allres <- rbind(allres, myerror)
    }
  }
}
write.csv(allres,'rank_sig1.csv',row.names = F)

# beep()
# mean(DISCOM.test.error)

# varying p ---------------------------------------------------------------
err <- NULL
estd <- NULL
allres <- NULL
r <- 3
smooth_p <- 1
missing_p <- 0
sig_prop <- 0.1
for (p in seq(100,500,100)){
  print(p)
  X.test <- read.csv(paste0('G:/tensor_data/fig2_n/',sprintf('t_sig_p_%.2f_Vcor_%d_missing_%.2f_p%d_%d.csv',sig_prop, smooth_p, missing_p,p,1)), header = F)
  Y.test <- read.csv(paste0('G:/tensor_data/fig2_n/',sprintf('t_sig_p_%.2f_Vcor_%d_missing_%.2f_p%d_y_%d.csv',sig_prop, smooth_p, missing_p,p,1)), header = F)
  X.test <- cleanx(X.test, NROW(Y.test))
  
  myerror <- NULL
  for (iter in c(1,1:100)){
    # print(iter)
    x <- read.csv(paste0('G:/tensor_data/fig2_n/',sprintf('sig_p_%.2f_Vcor_%d_missing_%.2f_p%d_%d.csv',sig_prop, smooth_p, missing_p,p,iter)), header = F)
    y <- read.csv(paste0('G:/tensor_data/fig2_n/',sprintf('sig_p_%.2f_Vcor_%d_missing_%.2f_p%d_y_%d.csv',sig_prop, smooth_p, missing_p,p,iter)), header = F)
    x <- cleanx(x, 500)
    
    X.tun <- x
    for (j in 1:ncol(x)){
      X.tun[is.na(X.tun[,j]),j] <- mean(X.tun[,j], na.rm=T)
    }
    Y.tun <- y
    Y.tun <- as.matrix(Y.tun)
    
    
    nalpha=10
    nlambda=30    
    alpha.all=10^seq(0,-3,length=nalpha)
    lambda.all=10^(seq(log10(2),-3,length=nlambda))
    
    DISCOM.tun.error=array(NA,dim=c(nalpha,nalpha,nlambda))
    
    xtx.raw=compute.xtx(x)
    xty=compute.xty(x,y)

    blockxtx <- list()
    for (blk in 1:(nrow(xtx.raw)/p)){
      blockxtx[[blk]] <- xtx.raw[((blk-1)*p+1):(blk*p),((blk-1)*p+1):(blk*p)]
    }
    xtx.raw.I <- as.matrix(do.call(bdiag, blockxtx))
    xtx.raw.C <- xtx.raw-xtx.raw.I
    
    pp <- nrow(xtx.raw)
    shrink.target <- sum(diag(xtx.raw))/pp
    
    if (TRUE){
      for(i in 1:nalpha){
        for(j in 1:nalpha){
          print(j)
          alpha1=alpha.all[i]
          alpha2=alpha.all[j]
          xtx=alpha1*xtx.raw.I+alpha2*xtx.raw.C+(1-alpha1)*shrink.target*diag(pp)
          if(min(eigen(xtx)$values)<0){DISCOM.tun.error[i,j,]=10^8}
          else{
            beta.initial=rep(0,pp)
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
    }
    
    opt.alpha1=alpha.all[opt.index[1]]
    opt.alpha2=alpha.all[opt.index[2]]
    opt.lambda=lambda.all[opt.index[3]]
    xtx=opt.alpha1*xtx.raw.I+opt.alpha2*xtx.raw.C+(1-opt.alpha1)*shrink.target*diag(pp)
    beta.cov=as.vector(crossProdLasso(xtx,xty,opt.lambda)$beta)
    
    predict.test.values=as.vector(as.matrix(X.test)%*%beta.cov)
    
    myerror <- c(myerror, 
                 sqrt(sum((as.matrix(Y.test)-predict.test.values)^2))/sqrt(NROW(Y.test)))
    # DISCOM.fpr=sum((beta.true==0)&(beta.cov!=0))/sum(beta.true==0)
    # DISCOM.fnr=sum((beta.true!=0)&(beta.cov==0))/sum(beta.true!=0)
    # DISCOM.est.error=sqrt(sum((beta.cov-beta.true)^2))
    # time=proc.time()-start
    # time_results[sim,7]=as.vector(time[3])
    
  }
  err <- c(err, median(myerror))
  estd <- c(estd, IQR(myerror))
  allres <- rbind(allres, myerror)
  
}
write.csv(allres,'varying_p.csv',row.names = F)

