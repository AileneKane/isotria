# Test model fitting on aritifical data with known kernel. 
# For use with data created using CreateTestData.R
# Select data files with the "fileNum" variable found on line 35(ish)

### The R working directory must be set to the folder this file is in. 
### Subfolders named ARTR and PSSP (with data sets) must exist, and  
### subfolders (created by CreateTestData.R) with the artificial data

rm(list=ls(all=TRUE))

require(fda); require(minqa); require(pracma); 
source("DistanceKernelSubs.R"); 
source("DistanceSplineSubs.R"); 

midRings=c(seq(1,19,2),seq(22.5,47.5,5),seq(55,145,10))
mk_breaks=function(n,xmin=1,xmax=150) {
  seq(xmin,xmax,length=n)
}

## define breakpoints for the spline basis, used in create.bspline.basis in fda package. 
breaks=c(mk_breaks(n=30, xmin=min(midRings),xmax=30), 40, 60, 80, 100, 120, max(midRings)); 

###### Specify which artificial data set will be analyzed #############
#  1xx = ARTR, 2xxx = PSSP (real data that the fake data are based on) 
#  x0x = Exp, x1x=Gaussian, x2x=3par, x3x=Sigmoid  (kernel used in fake data) 
##############################################################
analysisNum="Aug22PM"; #For identifying multiple runs with the same data

description="ARTR data Sigmoid" 
typeOfCurve="Sigmoid"

# lambdaList=10^seq(-6,8,by=1); # grid of smoothing parameters to try 
lambdaList=10^seq(-2,4,by=1); # just a few, to try it out quickly. 

for(fileNum in 111:111){  # just one file for an illustration 

datName=paste('Art data ', fileNum,': ',typeOfCurve,sep="") 
savName=paste("fit",fileNum,sep=""); 

########################################################################
# Read in the data 
#######################################################################
home=getwd(); 
setwd(sprintf("Artdata%s",fileNum)); # move to the dataset directory
fileName=paste("Artdata",as.character(fileNum),".Rdata",sep="")
load(fileName);

########################################################################################
# Making a place to save our figures and results 
# NOTE: THIS WILL OVERWRITE ANY PREVIOUS WORK WITH SAME DATA AND "ANALYSISNUM" IDENTIFIER
########################################################################################
dirName=paste("Dat",fileNum,"run",analysisNum,sep=''); # move to the dataset directory
unlink(dirName,recursive=TRUE)
dir.create(dirName)
setwd(dirName)   # THIS is where everything gets saved. NOTE: depends on setwd() above.  

########################################################
# CALL FITTING FUNCTION 
########################################################
report=TRUE; 
res<-fit_spline_aic(sppData.s=sppData.s,sppData.g=sppData.g,
               nonCompLength.s=nonCompLength.s, nonCompLength.g=nonCompLength.g,
               lambdaList=lambdaList,breaks=breaks,savName=savName,
               description=description,midRings=midRings)

Y=res$Y; Pars=res$Pars; aic=res$aic; Ymid=res$Ymid; dof=res$dof; dfit2=res$dfit2; emptyCols=res$emptyCols;

save(list=c("Y","Ymid","Pars","lambdaList","dfit2","aic","dof","emptyCols","warnings","breaks","description",
            "midRings","cur_dist_wts"), file=paste(savName,".Rdata",sep="")); 

##########################################################
# Everything past this point is just for making figures.  
# NOTE: if you fit at just a few lambda values, plotting might not work. 
##########################################################

closeDist=20 #maximum distance to plot for the "close" plots
colfun=colorRampPalette(c("blue","red"))
colors=colfun(ncol(Y))

# choose the best fits based on AIC 
aicPeak=max(which(peaks(-aic)==TRUE)); # smoothest local min 
span = seq(aicPeak-2,aicPeak+2,by=1); 
span= span[span>0]; span=span[span<=length(aic)]; 
N=min(3,length(span)); 
bestAICPos=match(sort(aic[span])[1:N],aic[span]) 
bestAICPos=sort(span[bestAICPos]); 

YBest=Ymid[,bestAICPos]; aicBest=aic[bestAICPos];
YMean=apply(YBest,1,mean); 

## scale to plot relative weights starting at d=1
kmid = cur_dist_wts(midRings,alpha); kmid=kmid/kmid[1]; 
for(k in 1:ncol(Ymid)) Ymid[,k]=Ymid[,k]/Ymid[1,k]; 
    
graphics.off(); 
dev.nu <- function(height,width) {
    dev.new(height=height,width=width); 
    par(mgp=c(2.25,1,0),cex.axis=1.4,cex.lab=1.4)
}

#### Make a clean plot for all distances
dev.nu(height=6,width=9); 
 
matplot(midRings,Ymid,type='l',lty=1,col=colors,main=paste("Actual fits for",datName),
        xlab='Distance',ylab='Relative interaction strength')
points(midRings,kmid,type="p",lwd=2,cex=1.5); 
legend("topright",inset=.05,legend=as.character(rev(lambdaList)),pch=19,col=rev(colors),title="Lambda",cex=.8)
abline(h=0,lty=2)
dev.copy2pdf(file=sprintf("SplineDat%s%sActual_all.pdf",doSpp,analysisNum)); 

#####Clean plot, only the near distances
graphics.off(); dev.nu(height=6,width=9); 
matplot(midRings,Ymid,type='l',lty=1,col=colors,main=paste("Actual fits for",datName),
        xlab='Distance',ylab='Relative interaction strength',xlim=c(0,closeDist))
legend("topright",inset=.05,legend=as.character(lambdaList),pch=19,col=colors,title="Lambda",cex=.8)
points(midRings,kmid,type="b",lwd=2,cex=1.25,lty=2); 
abline(h=0,lty=2)
dev.copy2pdf(file=sprintf("SplineDat%s%sActual_close.pdf",doSpp,analysisNum)); 

#### Plot of aic vs lambda
graphics.off(); dev.nu(height=6,width=9); 
plot(log10(lambdaList),aic,type='o',ylab='Approximate AIC',xlab='Lambda (log10 scale)',main='AIC vs. lambda')
dev.copy2pdf(file=sprintf("aic%s%s.pdf",doSpp,analysisNum));

#### Plot best N number of fits at near distances 
colors=colfun(N); 

graphics.off(); dev.nu(height=6,width=9); 
matplot(midRings,YBest,type='l',lty=1,col=colors,main=paste('Best fitting kernel for',datName),
        xlab='Distance (cm)',ylab='Relative interaction strength',xlim=c(0,closeDist),cex.lab=1.5,cex.main=1.5,lwd=3)
abline(h=0,lty=2)
legend("topright",inset=.05,legend=as.character(round(aicBest,digits=1)),pch=19,col=colors,title="AIC of fit",cex=1.5)
points(midRings,kmid,type="b",lwd=2,cex=1.25,lty=2); 
dev.copy2pdf(file=sprintf("SplineDat%s%sBest%s_close.pdf",doSpp,analysisNum,N)); 

### restore working directory to what it was 
setwd(home); 
}

