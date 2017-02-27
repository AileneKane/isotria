#This function runs FLM based on data that is formatted properly and nUnits (which determines the number of knots)
#This script returns both coefficients and probed values
FuncFLM<-function(data,nUnits,WEEK,SENS=TRUE){
  
  doFit <- function(data) {
    return(gam(logarea.t1 ~ logarea.t0 + s(lags, k=knots, by=tcovar, bs="cs") 
               + s(lags, k=knots, by=pcovar, bs="cs"),
               data=data, method="GCV.Cp",gamma=1.2))   
  }
  
  datag_all <- data#$datag; #if by chunks, datag
  
  knots=36  #this can be user defined.
    
  gam.ridge <- doFit(datag_all) 
  #summary(gam.ridge); 
  
  getBeta.data <- data.frame(lags=sort(unique(datag_all$lags)), pcovar=1, tcovar=0,Group=1,logarea.t0=0,W=0)
  terms.data <- predict(gam.ridge, newdata=getBeta.data, type="terms",se=TRUE)
  betas_mirrored <- cbind(lag=getBeta.data$lags,
                          beta=terms.data$fit[,"s(lags):pcovar"],
                          se=terms.data$se[,"s(lags):pcovar"])
  lags=betas_mirrored[,1]; beta_precip = betas_mirrored[,2]; se_precip = qnorm(0.975)*betas_mirrored[,3];
  
  getBeta.data <- data.frame(lags=sort(unique(datag_all$lags)), pcovar=0, tcovar=1, Group=1,logarea.t0=0,W=0)
  terms.data <- predict(gam.ridge, newdata=getBeta.data, type="terms",se=TRUE)
  betas_mirrored <- cbind(lag=getBeta.data$lags,
                          beta=terms.data$fit[,"s(lags):tcovar"],
                          se=terms.data$se[,"s(lags):tcovar"])
  lags=betas_mirrored[,1]; beta_temp = betas_mirrored[,2]; se_temp = qnorm(0.975)*betas_mirrored[,3];
  
  FLMlist<-list(beta_precip, se_precip, beta_temp, se_temp)
  names(FLMlist)<-c("beta_precip", "se_precip", "beta_temp", "se_temp")
  
  if(SENS){
    
    SensP=numeric(ncol(datag_all$pcovar)); newDat<-datag_all;
    p.orig<-as.vector(predict(object=gam.ridge, newdata=newDat,type="response"))
    
    for (ii in 1:ncol(datag_all$pcovar)){
      newDat<-datag_all;
      eps<-.05*mean(apply(datag_all$pcovar,2,sd))
      newDat$pcovar[,ii]<-datag_all$pcovar[,ii]+eps; 
      p.up<-as.vector(predict(object=gam.ridge, newdata=newDat,type="response"))
      SensP[ii]=(mean(p.up)-mean(p.orig))/eps; 
    }
    
    SensT=numeric(ncol(datag_all$tcovar)); 
    for (ii in 1:ncol(datag_all$tcovar)){
      newDat<-datag_all;
      eps<-.05*mean(apply(datag_all$tcovar,2,sd))
      newDat$tcovar[,ii]<-datag_all$tcovar[,ii]+eps; 
      p.up<-as.vector(predict(object=gam.ridge, newdata=newDat,type="response"))
      SensT[ii]=(mean(p.up)-mean(p.orig))/eps;
    }
    
  FLMlist<-list(beta_precip, se_precip, beta_temp, se_temp, SensP, SensT)
  names(FLMlist)<-c("beta_precip", "se_precip", "beta_temp", "se_temp", "SensP", "SensT")
  }
  
  return(FLMlist)}
