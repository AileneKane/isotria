#########################################################################
### Analysis of authetic climate/demographic data 
### Methods: FLM (splines), RF and Lasso
#########################################################################
rm(list=ls(all=TRUE)); graphics.off();
require(mgcv); require(randomForest); library(glmnet);

## Provide directory information
mainDir<-"C:/repos/drivers/Manuscripts/Sheffield/Rcode"
utilDir<-paste(mainDir,"/Climate/AuthenticAnalysis/Utilities", sep="")
dataDir<-paste(mainDir,"/Data", sep="")

## Set working directory to AuthenticAnalysis folder
setwd(mainDir)

## Import Demographic Data
sppList=sort(c("ARTR","HECO"))
doSpp=sppList[1]

#Aggregate by months or years
WEEKset=FALSE
if(WEEKset==TRUE){nUnits<-156}else{nUnits<-36} #how far back to you want to go (weeks/months)?

source(paste(utilDir,"/fetchDemoData.R", sep="")) #Synthesize Demog and Climate
data<-fetchDemoData(doSpp=doSpp, nUnits=nUnits, WEEK=WEEKset)

fit.noClimate <- lm(logarea.t1 ~ W + factor(Group), data=data); 
data$logarea.t1 <- data$logarea.t1 - predict(fit.noClimate,data=data);

#############################################
## Flm Routine
#############################################
source(paste(utilDir,"/FuncFLM.R", sep=""))
FLMRes<-FuncFLM(data=data,nUnits=nUnits, WEEK=WEEKset)

#############################################
## RF Routine
#############################################
y <- data$logarea.t1; 
evars <- which(substr(names(data),2,2)==".")
xvars <- c(evars); 
x <- as.data.frame(data[,xvars]);
x$la<-data$logarea.t0

source(paste(utilDir,"/FuncRF.R", sep=""))
RFRes<-FuncRF(x=x,y=y)

#############################################
# Lasso routine
#############################################
y <- data$logarea.t1
evars <- which(substr(names(data),2,2)==".")
zvars <- which(substr(names(data),8,10)=="t0"); 
xvars <- c(zvars,evars); 
x <- as.data.frame(data[,xvars]); 

source(paste(utilDir,"/FuncLas.R", sep=""))
LARes<-FuncLas(data=data,x=x,y=y)

#############################################
# Plotting routine
#############################################
colCF<-c("#FF0000", "#BF3EFF", "#009E73");
nUnits<-length(FLMRes$beta_precip)

par(mfrow=c(2,2), mar=c(2,2,1,1))
##Ppt Coefficients for FLM and LASSO
pymax<-max(FLMRes$beta_precip+FLMRes$se_precip);pymin<-min(FLMRes$beta_precip-FLMRes$se_precip)
plot(beta_precip~c(0:(nUnits-1)), data=FLMRes, type="l", col=colCF[1], ylim=c(pymin,pymax), 
     xlim=c(0,nUnits))
lines(beta_precip-se_precip~c(0:(nUnits-1)), data=FLMRes, type="l", lty=2, col=colCF[1])
lines(beta_precip+se_precip~c(0:(nUnits-1)), data=FLMRes, type="l", lty=2, col=colCF[1])
lines(LARes$pptCoef~c(0:(nUnits-1)), type="l", lty=1, col=colCF[3])
abline(h=0, lty=2)

##Temp Coefficients for FLM and LASSO
tymax<-max(FLMRes$beta_temp+FLMRes$se_temp);tymin<-min(FLMRes$beta_temp-FLMRes$se_temp)
plot(beta_temp~c(0:(nUnits-1)), data=FLMRes, type="l", ylim=c(tymin,tymax), 
     xlim=c(0,nUnits), col=colCF[1])
lines(beta_temp-se_temp~c(0:(nUnits-1)), data=FLMRes, type="l", lty=2, col=colCF[1])
lines(beta_temp+se_temp~c(0:(nUnits-1)), data=FLMRes, type="l", lty=2, col=colCF[1])
lines(LARes$tmpCoef~c(0:(nUnits-1)),type="l", lty=1, col=colCF[3])
abline(h=0, lty=2)

##Ppt Sensitivities for FLM, RF, and LASSO
pymax<-max(FLMRes$SensP);
plot(SensP~c(0:(nUnits-1)), data=FLMRes, type="l", col=colCF[1], #ylim=c(0,pymax), 
     xlim=c(0,nUnits))
lines(RFRes$SensP~c(0:(nUnits-1)), type="l", lty=1, col=colCF[2])
lines(LARes$SensP~c(1:(nUnits-1)), type="l", lty=1, col=colCF[3])
abline(h=0, lty=2)

pymax<-max(FLMRes$SensT);
plot(SensT~c(0:(nUnits-1)), data=FLMRes, type="l", col=colCF[1], #ylim=c(0,.01), 
     xlim=c(0,nUnits))
lines(RFRes$SensT~c(0:(nUnits-1)), type="l", lty=1, col=colCF[2])
lines(LARes$SensT~c(0:(nUnits-1)), type="l", lty=1, col=colCF[3])
abline(h=0, lty=2)

