rm(list=ls(all=TRUE)); #graphics.off();
require(mgcv); require(randomForest); library(glmnet); require(MASS);

#Directories 
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/repos"); # modify as needed
workDir=paste(root,"/drivers/Manuscripts/Sheffield/Rcode",sep=""); 
dataDir=paste(workDir,"/Data",sep="");
utilDir=paste(workDir,"/Climate/ArtificialAnalysis/Utilities",sep="")
setwd(paste(workDir,"/Climate/ArtificialAnalysis", sep=""))

set.seed(111)
##############################################################
#Write-in constants
##############################################################
#Simulation Parameters 
monthStart=6; nYrs=3; nMonths=nYrs*12; WEEK=FALSE;
nFyears=50; nIndiv=300 #3000/nFyears;

#Climate Data
DETREND=FALSE; SCALE=FALSE; 
noiseSD_logP=0.6; noiseSD_T= 2.1; 
tauHalf.T = 60; tauHalf.P = 15; rPT=0.5;

#Demographic Data
sppList=sort(c("ARTR","HECO"))
doSpp<-sppList[1]
alphaList <- c(0.2208319, 0.5672076, 0.5160774, 0.3499988);
midRings <- c(seq(1,19,2),seq(22.5,47.5,5),seq(55,145,10));

## Constants that affect simulated response variables
nEffect=24; sigmaFrac=1; fPrecip=.5; varFrac=.5;

paramList<-list(monthStart=monthStart,nYrs=nYrs,nMonths=nMonths,WEEK=WEEK,
             nFyears=nFyears,nIndiv=nIndiv,DETREND=DETREND,SCALE=SCALE,
             noiseSD_logP=noiseSD_logP,noiseSD_T=noiseSD_T,
             tauHalf.T=tauHalf.T,tauHalf.P=tauHalf.P,rPT=rPT,
             sppList=sppList,doSpp=doSpp,alphaList=alphaList,
             midRings=midRings,nEffect=nEffect,sigmaFrac=sigmaFrac,
             fPrecip=fPrecip,varFrac=varFrac)

#####################################
#Import artficial data
#####################################
source(paste(utilDir,"/FetchData.R", sep=""))
data<-FetchData(paramList=paramList)

par(mfrow=c(2,1))
plot(data$p.beta~c(1:length(data$p.beta)), type="o")
plot(data$t.beta~c(1:length(data$t.beta)), type="o")

#############################################
## Flm Routine
#############################################
source(paste(utilDir,"/FuncFLM.R", sep=""))

nUnits = ifelse(WEEK,nMonths*4,nMonths); 
FLMRes<-FuncFLM(data=data$datag,nUnits,WEEK=WEEK,SENS=TRUE)

#############################################
## RF Routine
#############################################
y <- data$datag$logarea.t1; 
evars <- which(substr(names(data$datag),2,2)==".")
xvars <- c(evars); 
x <- as.data.frame(data$datag[,xvars]);
x$la<-data$datag$logarea.t0

source(paste(utilDir,"/FuncRF.R", sep=""))
RFRes<-FuncRF(x=x,y=y)

#############################################
# Lasso routine
#############################################
y <- data$datag$logarea.t1
evars <- which(substr(names(data$datag),2,2)==".")
zvars <- which(names(data$datag)=="logarea.t0"); 
xvars <- c(zvars,evars); 
x <- as.data.frame(data$datag[,xvars]); 

source(paste(utilDir,"/FuncLas.R", sep=""))
LARes<-FuncLas(data=data$datag,x=x,y=y, SENS=TRUE)

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
lines(data$p.beta~c(1:length(data$p.beta)))

##Temp Coefficients for FLM and LASSO
tymax<-max(FLMRes$beta_temp+FLMRes$se_temp);tymin<-min(FLMRes$beta_temp-FLMRes$se_temp)
plot(beta_temp~c(0:(nUnits-1)), data=FLMRes, type="l", ylim=c(tymin,tymax), 
     xlim=c(0,nUnits), col=colCF[1])
lines(beta_temp-se_temp~c(0:(nUnits-1)), data=FLMRes, type="l", lty=2, col=colCF[1])
lines(beta_temp+se_temp~c(0:(nUnits-1)), data=FLMRes, type="l", lty=2, col=colCF[1])
lines(LARes$tmpCoef~c(0:(nUnits-1)),type="l", lty=1, col=colCF[3])
abline(h=0, lty=2)
lines(data$t.beta~c(1:length(data$t.beta)))

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

