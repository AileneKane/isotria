rm(list=ls(all=TRUE)); #graphics.off();
require(mgcv); require(randomForest); library(glmnet); require(MASS);

#Directories 
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/repos"); # modify as needed
workDir=paste(root,"/drivers/Manuscripts/Sheffield/Rcode",sep=""); 
dataDir=paste(workDir,"/Data",sep="");
utilDir=paste(workDir,"/Climate/ArtificialAnalysis/Utilities/",sep="")
simDatDir=paste(workDir,"/Climate/ArtificialAnalysis/SimData/",sep="")

setwd(paste(workDir,"/Climate/ArtificialAnalysis", sep=""))

tips<-10 #number of iterations of each set of params
ListID<-"MS_test" #For naming files later

##############################################################
#Write-in constants
##############################################################
#Simulation Parameters 
monthStart=6; nYrs=3; nMonths=nYrs*12; WEEK=c(FALSE);
#nFyears=30; nIndiv=3000/nFyears;

#Climate Data
DETREND=FALSE; SCALE=FALSE; 
noiseSD_logP=0.6; noiseSD_T= 2.1; 
tauHalf.T = 60; tauHalf.P = 60; rPT=c(0,0.25,.5);

#Demographic Data
sppList=sort(c("ARTR","HECO"))
doSpp<-sppList[1]
alphaList <- c(0.2208319, 0.5672076, 0.5160774, 0.3499988);
midRings <- c(seq(1,19,2),seq(22.5,47.5,5),seq(55,145,10));

## Constants that affect simulated response variables
nEffect=24; sigmaFrac=1; fPrecip=.5; 

varFrac=c(0.05,0.5);

#n=3000
#nIndiv=c(n/10,n/15,n/20,n/30,n/50)#individuals per year 
#nFyears=n/nIndiv+3

nIndiv=rep(100,5)
nFyears=c(10,15,20,30,50)+3
combo<-rep(seq(1:length(nFyears)), times=tips)

###########################################################
##Expanded grid for moving parts of big simulation
###########################################################
params<-expand.grid(doSpp,nMonths,combo, tauHalf.T, tauHalf.P, noiseSD_logP, noiseSD_T, nEffect, varFrac, sigmaFrac, fPrecip, rPT, WEEK)
names(params)<-c("doSpp","nMonths","combo","tauHalf.T", "tauHalf.P", "noiseSD_logP","noiseSD_T","nEffect", "varFrac","sigmaFrac", "fPrecip", "rPT", "WEEK")

#Make sure the number of individuals are held constant
for (i in 1:length(combo)){params$nFyears[params$combo==i]<-nFyears[i]}
for (i in 1:length(combo)){params$nIndiv[params$combo==i]<-nIndiv[i]}

params<-params[,-which(names(params)=="combo")]

write.table(params,paste(simDatDir,ListID,"_params.txt",sep=""),
            col.names=T, row.names=F)

nrow(params);((2000*nrow(params)/60)/60)/24;

###########################################################
##Big loop for big simulation
###########################################################
strt<-proc.time()
source(paste(utilDir,"AnalysisFun.R", sep="")) 
consList<-list(params,ListID,strt)

#it=1
#To run chunks STOP HERE and open "utiltites/AnalysisFun.R"

#####################################
#Loop for singular, or a few runs
#####################################

#system.time(for (it in 1){
#AnalysisFun(it=it,consList=consList)})

#####################################
#Parallel structure for many runs
#####################################
library(foreach)
library(doSNOW)
cl <- makeCluster(7, type = "SOCK", outfile="")
registerDoSNOW(cl)

listPak<-c('mgcv','randomForest','glmnet','MASS');
foreach(it=c(1:nrow(params)), .packages=listPak, .inorder=F) %dopar% AnalysisFun(it=it,consList=consList)

stopCluster(cl)

