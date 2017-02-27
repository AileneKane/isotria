################################################ 
# Model fitting to real data, unknown kernels
################################################
rm(list=ls(all=TRUE))
require(fda);  require(minqa); require(MASS); require(speedglm); 

## R's working directory must be the folder this file is in. 
## Subfolders named ARTR and PSSP (with data sets) must exist. 
## Data to analyze are chosen by the doSpp variable on line 20(ish). 

# setwd("C:/Repos/drivers/Manuscripts/Sheffield/Rcode/Competition") 

root=getwd(); 
source("LoadData.R"); 
source("DistanceSplineSubs.R"); 

doYear=FALSE; 
if(doYear) source("DistanceSplineSubs_T.R"); 

################### Specify which data set will be analyzed 
doSpp="ARTR"; datName=doSpp; dataDir=paste(root,"/",sep=""); 
description=paste("Fitting real data,",doSpp,sep="");
analysisNum="21AugAM"; #For identifying multiple runs with the same data

lambdaList=10^seq(-7,6,by=1); #roughness penalties to test
lambdaList=10^seq(-3,3,by=1); #short list to try it out quickly

midRings=c(seq(1,19,2),seq(22.5,47.5,5),seq(55,145,10))
mk_breaks=function(n,xmin=1,xmax=150) {
  seq(xmin,xmax,length=n)
}

# breaks for spline basis, used in create.bspline.basis in the fda package 
breaks=c(mk_breaks(n=36, xmin=min(midRings),xmax=13), 16, 18, 20, 25, 30, 40, 60, 80, 100, 120, max(midRings)); 

savName=paste("fit",datName,sep=""); 

######################################################
# Read in the data #
######################################################
outS=loadSurvivalData(dataDir,doSpp); 
D.s=outS$D; # data on survival and covariates 
sppData.s=outS$sppData; # annulus data only 
nonCompLength.s=outS$nonCompLength.s

outG=loadGrowthData(dataDir,doSpp); 
D.g=outG$D; # data on growth and covariates 
sppData.g=outG$sppData; # annulus data only 
nonCompLength.g=outG$nonCompLength.g

####################################################
# Making a place to save our figures and results   #
####################################################
#NOTE: THIS CURRENTLY OVERWRITES PREVIOUS WORK WITH SAME DATA AND "ANALYSISNUM" IDENTIFIER
dirName=paste(doSpp,analysisNum,sep=''); # move to the dataset directory
unlink(dirName,recursive=TRUE)
dir.create(dirName)
setwd(dirName)

#############################
# CALL NEW FITTING FUNCTION #
#############################
report=TRUE; 
res <- fit_spline_aic(sppData.s=sppData.s,sppData.g=sppData.g,
                    nonCompLength.s=nonCompLength.s, nonCompLength.g=nonCompLength.g,
                    lambdaList=lambdaList,breaks=breaks,savName=savName,
                    description=description,midRings=midRings)

Y=res$Y; Pars=res$Pars; aic=res$aic; Ymid=res$Ymid; dof=res$dof; dfit2=res$dfit2; emptyCols=res$emptyCols;

save(list=c("Y","Ymid","Pars","lambdaList","dfit2","aic","dof","emptyCols","warnings","breaks","description","dirName",
            "midRings","doSpp","analysisNum"), file=paste(savName,".Rdata",sep="")); 

##########################################################
# Everything past this point is just for making figures.  
# If you used just a few lambda values, it might not work. 
##########################################################

closeDist=15 #maximum distance to plot for the "close" plots
colfun=colorRampPalette(c("blue","red"))
colors=colfun(length(lambdaList))

graphics.off(); 
dev.nu <- function(height,width) {
    dev.new(height=height,width=width); 
    par(mgp=c(2.5,1,0),cex.axis=1.6,cex.lab=1.6,xaxs="i",bty="n",mar=c(5,5,2,2))
}

#### Make a clean plot for all distances
dev.nu(height=6,width=9); 
 
matplot(midRings,Ymid,type='l',lty=1,col=colors,main=paste("Fits for",doSpp),
        xlab='Distance',ylab='Relative interaction strength')
legend("topright",inset=.05,legend=as.character(rev(lambdaList)),pch=19,col=rev(colors),title="Lambda",cex=.8)
abline(h=0,lty=2)
dev.copy2pdf(file=sprintf("SplineDat%s%sActual_all.pdf",doSpp,analysisNum)); 

#####Clean plot, only the near distances
graphics.off(); dev.nu(height=6,width=9); 
matplot(midRings,Ymid,type='l',lty=1,col=colors,main=paste("Fits for",doSpp),
        xlab='Distance',ylab='Relative interaction strength',xlim=c(0,closeDist))
legend("topright",inset=.05,legend=as.character(lambdaList),pch=19,col=colors,title="Lambda",cex=.8)
abline(h=0,lty=2)
dev.copy2pdf(file=sprintf("SplineDat%s%sActual_close.pdf",doSpp,analysisNum)); 

#### Plot of aic vs lambda, etc. 
graphics.off(); dev.nu(height=9,width=9); 
par(mfrow=c(2,2)); 
plot(log10(lambdaList),aic,type='o',ylab='Approximate AIC',xlab='Lambda (log10 scale)',main='AIC vs. lambda')
plot(log10(lambdaList),dof,type='o',ylab='Approximate DOF',xlab='Lambda (log10 scale)',main='DOF vs. lambda')
plot(dof,aic,type='o',ylab='Approximate AIC',xlab='Approximate DOF',main='AIC vs. DOF')
dev.copy2pdf(file=sprintf("aic%s%s.pdf",doSpp,analysisNum));

#### Plot best N number of fits at near distances 
## smoothest peak of AIC as a function of dof 
findPeak=function(aic,dof) {   
    e = order(dof); 
    sortAIC= aic[order(dof)]; 
    aicPeaks = which(peaks(-sortAIC)==TRUE); 
    aicBest=sortAIC[aicPeaks[1]];     
    j = which(aic==aicBest); 
    return(j);
}  

# choose the best fits based on AIC 
e=order(dof); 
dofSave=dof; aicSave=aic; 
dof=dof[e]; aic=aic[e]; Ymid=Ymid[,e]; 
aicPeak=findPeak(aic,dof); 

# choose the best fits based on AIC 
aicPeak=findPeak(aic,dof); 
span = seq(aicPeak-2,aicPeak+2,by=1); 
span= span[span>0]; span=span[span<=length(aic)]; 
N=min(3,length(span)); 
#bestAICPos=match(sort(aic[span])[1:N],aic[span]) 
# bestAICPos=sort(span[bestAICPos]); 

aicmin=aic[aicPeak];  
bestAICPos=c(aicPeak, which((aic>aicmin)&(aic<(4+aicmin))));

e=order(aic[bestAICPos]); bestAICPos=bestAICPos[e]; 
YBest=Ymid[,bestAICPos]; aicBest=aic[bestAICPos];
YMean=apply(YBest,1,mean); 

graphics.off(); dev.nu(height=6,width=9); 
colors=c("black",colfun(length(bestAICPos)-1));  
matplot(midRings,YBest,type='b',pch=1,col=colors, lty=c(1:20), cex=1.3,
        xlab='Distance (cm)',ylab='Relative interaction strength',xlim=c(0,closeDist+0.1),
        cex.lab=1.75,cex.main=1.75,lwd=c(3,rep(2,20)))
matpoints(midRings,YBest[,1],type='l',lwd=3)

abline(h=0,lty=2)
legend("topright",inset=.15,legend=as.character(round(aicBest,digits=1)),pch=1,col=colors,title=doSpp,cex=1.65,
    lty=c(1:20),lwd=2,bty="n")
dev.copy2pdf(file=sprintf("SplineDat%s%sBest%s_close.pdf",doSpp,analysisNum,N)); 

dof=dofSave; aic=aicSave; 

# save.image(file="allInfo.Rdata"); # optional, to save all the results 