#Figures and individual life history traits estimated from posterior samples of multistate model for Isotria medioloides Alton NH population
#Data provided by Bill Brumback
#Coding by Ailene Ettinger with help frmo Andy Roly and Elizabeth Crone
#this file has code for all figures in the manuscript, and for estimating life history traits 
#(lifepsan, proportion dormant, length of dormancy)
#setwd("~/isotria") #at usgs
#setwd("/Users/aileneettinger/git/isotria/analyses")
rm(list=ls()) 
options(stringsAsFactors=FALSE)
###Figure 1
isoinds<-read.csv("Isotria_Stage3_2016.csv", header=T)#
isoinds<-isoinds[-which(isoinds$UniqueID=="X-01-244"),]#individual has to be removed because of errors in its monitoring
#head(isoinds)
#Add column for emergent/nonemergent
isoinds$Emerg<-NA
isoinds[which(isoinds$TotNoStems>0),]$Emerg=1
isoinds[which(isoinds$TotNoStems==0),]$Emerg=0
##Select out just groups X and Y for this analysis
isoindsX=isoinds[isoinds$Group=="X",]
isoindsY=isoinds[isoinds$Group=="Y",]
isoindsXY=isoinds[isoinds$Group=="X"|isoinds$Group=="Y",]
isoindsXY$UniqueID=factor(isoindsXY$UniqueID)
#dim(isoindsXY)
isoindsXY$Group=factor(isoindsXY$Group)
##get isotria data into format such that 1=veg, 2=rep and 3=not seen
#to do this, Add column for reproductive (=arrested, flowering, or fruiting)/not rep
isoindsXY$Repro<-NA
isoindsXY[which(isoindsXY$NoFrStems>0|isoindsXY$NoFlStems>0|isoindsXY$NoArrStems>0),]$Repro=1
isoindsXY[which(isoindsXY$NoFrStems==0&isoindsXY$NoFlStems==0&isoindsXY$NoArrStems==0),]$Repro=0
isoindsX$Repro<-NA
isoindsX[which(isoindsX$NoFrStems>0|isoindsX$NoFlStems>0|isoindsX$NoArrStems>0),]$Repro=1
isoindsX[which(isoindsX$NoFrStems==0&isoindsX$NoFlStems==0&isoindsX$NoArrStems==0&isoindsX$Stage!="D"),]$Repro=0
isoindsY$Repro<-NA
isoindsY[which(isoindsY$NoFrStems>0|isoindsY$NoFlStems>0|isoindsY$NoArrStems>0),]$Repro=1
isoindsY[which(isoindsY$NoFrStems==0 & isoindsY$NoFlStems==0 & isoindsY$NoArrStems==0&isoindsY$Stage!="D"),]$Repro=0

quartz(height=6,width=10)
par(mfrow=c(1,1),mar=c(1,5,1,.5), oma=c(5,.5,.5,.5))
pop_x<-tapply(isoindsX$UniqueID,list(isoindsX$Year,isoindsX$Repro),length)#how does pop change over time, by group?
pop_y<-tapply(isoindsY$UniqueID,list(isoindsY$Year,isoindsY$Repro),length)#how does pop change over time, by group?
pop_x[which(is.na(pop_x))]=0#replace NAs with 0s
pop_y[which(is.na(pop_y))]=0#replace NAs with 0s

plot(pop_x[,1]~rownames(pop_x),type="l",ylab="# Individuals Observed", xlab="Year", xaxt="n", xlim=c(1985,2025),ylim=c(0,60), bty="l", lty=3,col="black", lwd=2, cex.axis=1.3, cex.lab=1.5)
lines(pop_x[,2]~rownames(pop_x), lty=1, lwd=2)
lines(pop_y[,1]~rownames(pop_y), lty=3,col="darkgray", lwd=2)
lines(pop_y[,2]~rownames(pop_y), lty=1, col="darkgray", lwd=2)
text(2015, pop_x[31,1]+1,labels="Control Group, X",adj=0,cex=1.1)
text(2015.2, pop_x[31,1]-1.5,labels="(Vegetative)",adj=0,cex=1.1)

text(2015, pop_y[31,1],labels="Cleared Group, Y",adj=0,cex=1.1)
text(2015.2, pop_y[31,1]-2.5,labels="(Vegetative)",adj=0,cex=1.1)

text(2015, pop_x[31,2],labels="Control Group, X",adj=0,cex=1.1)
text(2015.2, pop_x[31,2]-2.5,labels="(Reproductive)",adj=0,cex=1.1)

text(2015,pop_y[31,2],labels="Cleared Group, Y",adj=0,cex=1.1)
text(2015.2, pop_y[31,2]-2.5,labels="(Reproductive)",adj=0,cex=1.1)

abline(v=1997,lty=2,col="gray", lwd=3)
axis(side=1,at=rownames(pop_x), labels=TRUE, cex.axis=1.3)
mtext("Year",side=1, adj=.35, cex=1.5, line=2.5)
### select out vital rates to calculate dwell times, etc
library(popbio)
mod.samples<-read.csv("msmod_samples_complex.csv", header=T)
# vital rates for group X prior to clearing
phiV_Xpre<-mod.samples[,which(colnames(mod.samples)=="phiA0.1.")]#if not read in, colname=psiA0[1]
phiR_Xpre<-mod.samples[,which(colnames(mod.samples)=="phiB0.1.")]#if not read in, colname=psiB0[1]
pdormV_Xpre<-1-mod.samples[,which(colnames(mod.samples)=="pA0.1.")]#if not read in, colname=pA0[1]
pdormR_Xpre<-1-mod.samples[,which(colnames(mod.samples)=="pB0.1.")]#if not read in, colname=pB0[1]
veg.rep_Xpre<- mod.samples[,which(colnames(mod.samples)=="psiA0.1.")]#if not read in, colname=phiA0[1]
rep.veg_Xpre<-mod.samples[,which(colnames(mod.samples)=="psiB0.1.")]#if not read in, colname=phiB0[1]

# vital rates for group Y prior to clearing
phiV_Ypre<-mod.samples[,which(colnames(mod.samples)=="phiA0.2.")]#if not read in, colname=psiA0[2]
phiR_Ypre<-mod.samples[,which(colnames(mod.samples)=="phiB0.2.")]#if not read in, colname=psiB0[2]
pdormV_Ypre<-1-mod.samples[,which(colnames(mod.samples)=="pA0.2.")]#if not read in, colname=pA0[2]
pdormR_Ypre<-1-mod.samples[,which(colnames(mod.samples)=="pB0.2.")]#if not read in, colname=pB0[2]
veg.rep_Ypre<- mod.samples[,which(colnames(mod.samples)=="psiA0.2.")]#if not read in, colname=phiA0[2]
rep.veg_Ypre<-mod.samples[,which(colnames(mod.samples)=="psiB0.2.")]#if not read in, colname=phiB0[2]

# vital rates for group X after clearing
phiV_Xpost<-mod.samples[,which(colnames(mod.samples)=="phiA1.1.")]#if not read in, colname=psiA1[1]
phiR_Xpost<-mod.samples[,which(colnames(mod.samples)=="phiB1.1.")]#if not read in, colname=psiB1[1]
pdormV_Xpost<-1-mod.samples[,which(colnames(mod.samples)=="pA1.1.")]#if not read in, colname=pA1[1]
pdormR_Xpost<-1-mod.samples[,which(colnames(mod.samples)=="pB1.1.")]#if not read in, colname=pB1[1]
veg.rep_Xpost<- mod.samples[,which(colnames(mod.samples)=="psiA1.1.")]#if not read in, colname=phiA1[1]
rep.veg_Xpost<-mod.samples[,which(colnames(mod.samples)=="psiB1.1.")]#if not read in, colname=phiB1[1]

# vital rates for group Y after clearing
phiV_Ypost<-mod.samples[,which(colnames(mod.samples)=="phiA1.2.")]#if not read in, colname=psiA1[2]
phiR_Ypost<-mod.samples[,which(colnames(mod.samples)=="phiB1.2.")]#if not read in, colname=psiB1[2]
pdormV_Ypost<-1-mod.samples[,which(colnames(mod.samples)=="pA1.2.")]#if not read in, colname=pA1[2]
pdormR_Ypost<-1-mod.samples[,which(colnames(mod.samples)=="pB1.2.")]#if not read in, colname=pB1[2]
veg.rep_Ypost<- mod.samples[,which(colnames(mod.samples)=="psiA1.2.")]#if not read in, colname=phiA1[2]
rep.veg_Ypost<-mod.samples[,which(colnames(mod.samples)=="psiB1.2.")]#if not read in, colname=phiB1[2]

###Porportion of plants dormant in each condition
get.propdorm <- function(phiV,veg.rep,pdormV,phiR,rep.veg,pdormR) {
  prop.dorm= array()
  for (i in 1:length(phiV)){
    tmx = c(phiV[i]*(1-veg.rep[i])*pdormV[i], phiR[i]*rep.veg[i]*pdormV[i], phiV[i]*(1-veg.rep[i])*pdormV[i], phiR[i]*rep.veg[i]*pdormV[i],
            phiV[i]*veg.rep[i]*pdormR[i], phiR[i]*(1-rep.veg[i])*pdormR[i], phiV[i]*veg.rep[i]*pdormR[i], phiR[i]*(1-rep.veg[i])*pdormR[i],
            phiV[i]*(1-veg.rep[i])*(1-pdormV[i]), phiR[i]*rep.veg[i]*(1-pdormV[i]), phiV[i]*(1-veg.rep[i])*(1-pdormV[i]), phiR[i]*rep.veg[i]*(1-pdormV[i]),
            phiV[i]*veg.rep[i]*(1-pdormR[i]), phiR[i]*(1-rep.veg[i])*(1-pdormR[i]), phiV[i]*veg.rep[i]*(1-pdormR[i]), phiR[i]*(1-rep.veg[i])*(1-pdormR[i]))
    tmx = matrix(tmx, nrow = 4, byrow = T)
    eigen.analysis(tmx)$stable.stage
    prop.dorm[i] = sum(eigen.analysis(tmx)$stable.stage[1:2])
  }
  return(prop.dorm)#
  
}
#even though clearing does not change the probability of dormancy per se, it could change the expected proportion of dormant plants via changes in other vital rates.
propdorm_Xpre<-get.propdorm(phiV_Xpre,veg.rep_Xpre,pdormV_Xpre,phiR_Xpre,rep.veg_Xpre,pdormR_Xpre)
propdorm_Ypre<-get.propdorm(phiV_Ypre,veg.rep_Ypre,pdormV_Ypre,phiR_Ypre,rep.veg_Ypre,pdormR_Ypre)
propdorm_Xpost<-get.propdorm(phiV_Xpost,veg.rep_Xpost,pdormV_Xpost,phiR_Xpost,rep.veg_Xpost,pdormR_Xpost)
propdorm_Ypost<-get.propdorm(phiV_Ypost,veg.rep_Ypost,pdormV_Ypost,phiR_Ypost,rep.veg_Ypost,pdormR_Ypost)
#windows(height=6,width=10)
#quartz(height=6,width=10)
#par(mfrow=c(2,2))
#hist(propdorm_Xpre, xlim=c(0,1))
#hist(propdorm_Ypre,xlim=c(0,1))
hist(propdorm_Xpost,xlim=c(0,1))
hist(propdorm_Ypost,xlim=c(0,1))

mean(propdorm_Xpre);sd(propdorm_Xpre)#0.257 (0.052 plants dormant in uncleared prior to clearing
mean(propdorm_Ypre);sd(propdorm_Ypre)#0.219 (0.050)plants dormant in cleared prior to clearing
mean(propdorm_Xpost);sd(propdorm_Xpost)#0.10 (0.06) plants dormant in uncleared post clearing
mean(propdorm_Ypost);sd(propdorm_Ypost)#0.094 (0.069) plants dormant in cleared post clearing

####Now life expectancy:
##test:
#phiV=phiV_Ypost
#veg.rep=veg.rep_Ypost
#pdormV=pdormV_Ypost
#phiR=phiR_Ypost
#rep.veg=rep.veg_Ypost
#pdormR=pdormR_Ypost

#to figure out effect of survival on lifepsand estimates, plug in mean values for everything then try changing phi:
#phiV=0.9999
#veg.rep=0.47
#pdormV=0.30
#phiR=0.99
#rep.veg=0.02
#pdormR=0.024
#with theabove mean parameters, lifespan_med is 30. 
#if i change phiV to 0.98, lifepsan_med is 35
#to 0.99, liefepsan=40; change of phiV frmo .99 to .999 moves lifespan from  69 to 74
#both phiV and phiR changed to .99; lifespan goes up to 69;
#with PhiV at .99 and when phiR changed frmo .99 to .999-.9933, med lifespan=Inf
#with PhiV at .99 and when phiR .991, med lifespan=77
#with PhiV at .99 and when phiR .992, med lifespan=85
#with PhiV at .99 and when phiR .993, med lifespan=97
#with PhiV at .99 and when phiR .9931, med lifespan=98
#with PhiV at .99 and when phiR .9932, med lifespan=99
#with PhiV at .99 and when phiR .99325-8, med lifespan=100

get.lifespan<- function(phiV,veg.rep,pdormV,phiR,rep.veg,pdormR){
  lifespan_med= array()
  #lifespan_95th= array()
  #lifespan_rep_med= array()
  #lifespan_rep_95th= array()
  #nyrs_fl= array()
  for (i in 1:length(phiV)){
    tmx = c(phiV[i]*(1-veg.rep[i])*pdormV[i], phiR[i]*rep.veg[i]*pdormV[i], phiV[i]*(1-veg.rep[i])*pdormV[i], phiR[i]*rep.veg[i]*pdormV[i],
            phiV[i]*veg.rep[i]*pdormR[i], phiR[i]*(1-rep.veg[i])*pdormR[i], phiV[i]*veg.rep[i]*pdormR[i], phiR[i]*(1-rep.veg[i])*pdormR[i],
            phiV[i]*(1-veg.rep[i])*(1-pdormV[i]), phiR[i]*rep.veg[i]*(1-pdormV[i]), phiV[i]*(1-veg.rep[i])*(1-pdormV[i]), phiR[i]*rep.veg[i]*(1-pdormV[i]),
            phiV[i]*veg.rep[i]*(1-pdormR[i]), phiR[i]*(1-rep.veg[i])*(1-pdormR[i]), phiV[i]*veg.rep[i]*(1-pdormR[i]), phiR[i]*(1-rep.veg[i])*(1-pdormR[i]))
    tmx = matrix(tmx, nrow = 4, byrow = T)
    ##### one way to calculate life span - calculate the probability of still being alive i years into the future
    n0 = c(0,0,1000,0)#lifespan starting from vegetative
    nsum = array()
    flwrsum = array()
    for(j in 1:1800){
      n1 = tmx%*%n0
      nsum[j] = sum(n1)
      flwrsum[j] = n1[4]
      n0 = n1
    }#
    lifespan_med[i]= min(which(nsum <900)) # this is actually the median survival time
    #nyrs_fl[i]=sum(flwrsum)/1000 # number of years flowering, over an average lifetime = 1.9 without clearing, 11.6 with
  }
  return (lifespan_med)
}
lifespan_Xpre<-get.lifespan(phiV_Xpre,veg.rep_Xpre,pdormV_Xpre,phiR_Xpre,rep.veg_Xpre,pdormR_Xpre)
lifespan_Ypre<-get.lifespan(phiV_Ypre,veg.rep_Ypre,pdormV_Ypre,phiR_Ypre,rep.veg_Ypre,pdormR_Ypre)
lifespan_Xpost<-get.lifespan(phiV_Xpost,veg.rep_Xpost,pdormV_Xpost,phiR_Xpost,rep.veg_Xpost,pdormR_Xpost)
lifespan_Ypost<-get.lifespan(phiV_Ypost,veg.rep_Ypost,pdormV_Ypost,phiR_Ypost,rep.veg_Ypost,pdormR_Ypost)
lifespan_Xpost2<-lifespan_Xpost[-(which(lifespan_Xpost=="Inf"))]
lifespan_Ypost2<-lifespan_Ypost[-(which(lifespan_Ypost=="Inf"))]

#Alternatively, life expectancy can be calculated as -1/ln(s)
LEV_Xpre<--1/(log(phiV_Xpre))#median=5.5
LEV_Xpost<--1/(log(phiV_Xpost))#median=7.69
LEV_Ypre<--1/(log(phiV_Ypre))#median=6.31
LEV_Ypost<--1/(log(phiV_Ypost))#median=60

get.lifespan_flow<- function(phiV,veg.rep,pdormV,phiR,rep.veg,pdormR){
  lifespan_med= array()
  #lifespan_95th= array()
  #lifespan_rep_med= array()
  #lifespan_rep_95th= array()
  #nyrs_fl= array()
  for (i in 1:length(phiV)){
    tmx = c(phiV[i]*(1-veg.rep[i])*pdormV[i], phiR[i]*rep.veg[i]*pdormV[i], phiV[i]*(1-veg.rep[i])*pdormV[i], phiR[i]*rep.veg[i]*pdormV[i],
            phiV[i]*veg.rep[i]*pdormR[i], phiR[i]*(1-rep.veg[i])*pdormR[i], phiV[i]*veg.rep[i]*pdormR[i], phiR[i]*(1-rep.veg[i])*pdormR[i],
            phiV[i]*(1-veg.rep[i])*(1-pdormV[i]), phiR[i]*rep.veg[i]*(1-pdormV[i]), phiV[i]*(1-veg.rep[i])*(1-pdormV[i]), phiR[i]*rep.veg[i]*(1-pdormV[i]),
            phiV[i]*veg.rep[i]*(1-pdormR[i]), phiR[i]*(1-rep.veg[i])*(1-pdormR[i]), phiV[i]*veg.rep[i]*(1-pdormR[i]), phiR[i]*(1-rep.veg[i])*(1-pdormR[i]))
    tmx = matrix(tmx, nrow = 4, byrow = T)
    ##### one way to calculate life span - calculate the probability of still being alive i years into the future
    n0 = c(0,0,0,1000)#lifespan starting from vegetative
    nsum = array()
    flwrsum = array()
    for(j in 1:1800){
      n1 = tmx%*%n0
      nsum[j] = sum(n1)
      flwrsum[j] = n1[4]
      n0 = n1
    }#
    lifespan_med[i]= min(which(nsum <900)) # 
    #nyrs_fl[i]=sum(flwrsum)/1000 # number of years flowering, over an average lifetime = 1.9 without clearing, 11.6 with
  }
  return (lifespan_med)
}
lifespan_flow_Xpre<-get.lifespan_flow(phiV_Xpre,veg.rep_Xpre,pdormV_Xpre,phiR_Xpre,rep.veg_Xpre,pdormR_Xpre)
lifespan_flow_Ypre<-get.lifespan_flow(phiV_Ypre,veg.rep_Ypre,pdormV_Ypre,phiR_Ypre,rep.veg_Ypre,pdormR_Ypre)
lifespan_flow_Xpost<-get.lifespan_flow(phiV_Xpost,veg.rep_Xpost,pdormV_Xpost,phiR_Xpost,rep.veg_Xpost,pdormR_Xpost)
lifespan_flow_Ypost<-get.lifespan_flow(phiV_Ypost,veg.rep_Ypost,pdormV_Ypost,phiR_Ypost,rep.veg_Ypost,pdormR_Ypost)
lifespan_flow_Xpost2<-lifespan_flow_Xpost[-(which(lifespan_flow_Xpost=="Inf"))]
lifespan_flow_Ypost2<-lifespan_flow_Ypost[-(which(lifespan_flow_Ypost=="Inf"))]

###Length of each bout of dormancy
# even though clearing does not change the probability of dormancy per se, it could change the expected proportion of dormant plants via changes in other vital rates.
get.lengthdorm <- function(phiV,veg.rep,pdormV,phiR,rep.veg,pdormR){
  mydorm_all=matrix(data=NA,nrow=length(phiV),ncol=11,byrow=TRUE)
  mnlengthdor=array()
  for (i in 1:length(phiV)){
    tmx.dorm = c(phiV[i]*(1-veg.rep[i])*pdormV[i], phiR[i]*rep.veg[i]*pdormV[i], phiV[i]*(1-veg.rep[i])*pdormV[i], phiR[i]*rep.veg[i]*pdormV[i],
                 phiV[i]*veg.rep[i]*pdormR[i], phiR[i]*(1-rep.veg[i])*pdormR[i], phiV[i]*veg.rep[i]*pdormR[i], phiR[i]*(1-rep.veg[i])*pdormR[i],
                 0,0,0,0,0,0,0,0)
    tmx.dorm = matrix(tmx.dorm, nrow = 4, byrow = T)
    # length of dormancy starting from dormant veg
    n0 = c(1000,0,0,0)
    nsum = array()
    for(j in 1:100){
      n1 = tmx.dorm%*%n0
      nsum[j] = sum(n1)
      n0 = n1
    }
    mydorm = c(1, nsum[1:10]/1000)/(1+sum(nsum)/1000)
    mydorm_all[i,]= mydorm
    numinds<-mydorm*1000
    dormls<-array()
    for (k in 1:length(numinds)){
      inddormls<-c(rep(k,times=numinds[k]))
      dormls<-c(dormls,inddormls)
    }
    mnlengthdor[i]<-mean(dormls, na.rm=T)
  }
  return(mnlengthdor)#
}
lengthdorm_Xpre<-get.lengthdorm(phiV_Xpre,veg.rep_Xpre,pdormV_Xpre,phiR_Xpre,rep.veg_Xpre,pdormR_Xpre)
lengthdorm_Ypre<-get.lengthdorm(phiV_Ypre,veg.rep_Ypre,pdormV_Ypre,phiR_Ypre,rep.veg_Ypre,pdormR_Ypre)
lengthdorm_Xpost<-get.lengthdorm(phiV_Xpost,veg.rep_Xpost,pdormV_Xpost,phiR_Xpost,rep.veg_Xpost,pdormR_Xpost)
lengthdorm_Ypost<-get.lengthdorm(phiV_Ypost,veg.rep_Ypost,pdormV_Ypost,phiR_Ypost,rep.veg_Ypost,pdormR_Ypost)

###Now, calculate length of each bout of dormancy and proportion dormant plant, starting with reproductive plants
get.lengthdorm_flow <- function(phiV,veg.rep,pdormV,phiR,rep.veg,pdormR){
  mnlengthdor=array()
  mydorm_all=matrix(data=NA,nrow=length(phiV),ncol=11,byrow=TRUE)
  for (i in 1:length(phiV)){
    tmx.dorm = c(phiV[i]*(1-veg.rep[i])*pdormV[i], phiR[i]*rep.veg[i]*pdormV[i], phiV[i]*(1-veg.rep[i])*pdormV[i], phiR[i]*rep.veg[i]*pdormV[i],
                 phiV[i]*veg.rep[i]*pdormR[i], phiR[i]*(1-rep.veg[i])*pdormR[i], phiV[i]*veg.rep[i]*pdormR[i], phiR[i]*(1-rep.veg[i])*pdormR[i],
                 0,0,0,0,0,0,0,0)
    tmx.dorm = matrix(tmx.dorm, nrow = 4, byrow = T)
    # length of dormancy starting from dormant flowering
    n0_flow = c(0,1000,0,0)
    nsum_flow = array()
    for(j in 1:100){
      n1_flow = tmx.dorm%*%n0_flow
      nsum_flow[j] = sum(n1_flow)
      n0_flow = n1_flow
    }
    mydorm_flow = c(1, nsum_flow[1:10]/1000)/(1+sum(nsum_flow)/1000)
    #prop_dorm1yr_flow[i]= mydorm_flow[1]
    mydorm_all[i,]= mydorm_flow
    numinds<-mydorm_flow*1000
    dormls<-array()
    for (k in 1:length(numinds)){
      inddormls<-c(rep(k,times=numinds[k]))
      dormls<-c(dormls,inddormls)
    }
    mnlengthdor[i]<-mean(dormls, na.rm=T)
  }
  return(mnlengthdor)#
}
lengthdorm_flow_Xpre<-get.lengthdorm_flow(phiV_Xpre,veg.rep_Xpre,pdormV_Xpre,phiR_Xpre,rep.veg_Xpre,pdormR_Xpre)
lengthdorm_flow_Ypre<-get.lengthdorm_flow(phiV_Ypre,veg.rep_Ypre,pdormV_Ypre,phiR_Ypre,rep.veg_Ypre,pdormR_Ypre)
lengthdorm_flow_Xpost<-get.lengthdorm_flow(phiV_Xpost,veg.rep_Xpost,pdormV_Xpost,phiR_Xpost,rep.veg_Xpost,pdormR_Xpost)
lengthdorm_flow_Ypost<-get.lengthdorm_flow(phiV_Ypost,veg.rep_Ypost,pdormV_Ypost,phiR_Ypost,rep.veg_Ypost,pdormR_Ypost)

##Figures
#2x2table for each vital rate with first column control, second column logged
#if model not loaded, then use model sample files to get estinat
ms3a<-read.csv("isotria2stagemodsum_complex.csv", header=T)
rownames(ms3a)<-ms3a[,1]
surv_veg<-as.data.frame(rbind(ms3a$mean[grep("phiA0",substr(rownames(ms3a),1,5))],ms3a$mean[grep("phiA1",substr(rownames(ms3a),1,5))]))
surv_rep<-as.data.frame(rbind(ms3a$mean[grep("phiB0",substr(rownames(ms3a),1,5))],ms3a$mean[grep("phiB1",substr(rownames(ms3a),1,5))]))
emer_veg<-as.data.frame(rbind(ms3a$mean[grep("pA0",substr(rownames(ms3a),1,3))],ms3a$mean[grep("pA1",substr(rownames(ms3a),1,3))]))
emer_rep<-as.data.frame(rbind(ms3a$mean[grep("pB0",substr(rownames(ms3a),1,3))],ms3a$mean[grep("pB1",substr(rownames(ms3a),1,3))]))
trans_vr<-as.data.frame(rbind(ms3a$mean[grep("psiA0",substr(rownames(ms3a),1,5))],ms3a$mean[grep("psiA1",substr(rownames(ms3a),1,5))]))
trans_rv<-as.data.frame(rbind(ms3a$mean[grep("psiB0",substr(rownames(ms3a),1,5))],ms3a$mean[grep("psiB1",substr(rownames(ms3a),1,5))]))

surv_veg_med<-as.data.frame(rbind(ms3a$X50.[grep("phiA0",substr(rownames(ms3a),1,5))],ms3a$X50.[grep("phiA1",substr(rownames(ms3a),1,5))]))
surv_rep_med<-as.data.frame(rbind(ms3a$X50.[grep("phiB0",substr(rownames(ms3a),1,5))],ms3a$X50.[grep("phiB1",substr(rownames(ms3a),1,5))]))
emer_veg_med<-as.data.frame(rbind(ms3a$X50.[grep("pA0",substr(rownames(ms3a),1,3))],ms3a$X50.[grep("pA1",substr(rownames(ms3a),1,3))]))
emer_rep_med<-as.data.frame(rbind(ms3a$X50.[grep("pB0",substr(rownames(ms3a),1,3))],ms3a$X50.[grep("pB1",substr(rownames(ms3a),1,3))]))
trans_vr_med<-as.data.frame(rbind(ms3a$X50.[grep("psiA0",substr(rownames(ms3a),1,5))],ms3a$X50.[grep("psiA1",substr(rownames(ms3a),1,5))]))
trans_rv_med<-as.data.frame(rbind(ms3a$X50.[grep("psiB0",substr(rownames(ms3a),1,5))],ms3a$X50.[grep("psiB1",substr(rownames(ms3a),1,5))]))

colnames(surv_veg)<-c("control","logged")
colnames(surv_rep)<-c("control","logged")
colnames(emer_veg)<-c("control","logged")
colnames(emer_rep)<-c("control","logged")
colnames(trans_vr)<-c("control","logged")
colnames(trans_rv)<-c("control","logged")

colnames(surv_veg_med)<-c("control","logged")
colnames(surv_rep_med)<-c("control","logged")
colnames(emer_veg_med)<-c("control","logged")
colnames(emer_rep_med)<-c("control","logged")
colnames(trans_vr_med)<-c("control","logged")
colnames(trans_rv_med)<-c("control","logged")

##use code below if model not loaded:
surv_veg_q2.5<-c(ms3a$X2.5.[grep("phiA0",substr(rownames(ms3a),1,5))],ms3a$X2.5.[grep("phiA1",substr(rownames(ms3a),1,5))])
surv_rep_q2.5<-c(ms3a$X2.5.[grep("phiB0",substr(rownames(ms3a),1,5))],ms3a$X2.5.[grep("phiB1",substr(rownames(ms3a),1,5))])
trans_vr_q2.5<-c(ms3a$X2.5.[grep("psiA0",substr(rownames(ms3a),1,5))],ms3a$X2.5.[grep("psiA1",substr(rownames(ms3a),1,5))])
trans_rv_q2.5<-c(ms3a$X2.5.[grep("psiB0",substr(rownames(ms3a),1,5))],ms3a$X2.5.[grep("psiB1",substr(rownames(ms3a),1,5))])
emer_veg_q2.5<-c(ms3a$X2.5.[grep("pA0",substr(rownames(ms3a),1,3))],ms3a$X2.5.[grep("pA1",substr(rownames(ms3a),1,3))])
emer_rep_q2.5<-c(ms3a$X2.5.[grep("pB0",substr(rownames(ms3a),1,3))],ms3a$X2.5.[grep("pB1",substr(rownames(ms3a),1,3))])

surv_veg_q97.5<-c(ms3a$X97.5.[grep("phiA0",substr(rownames(ms3a),1,5))],ms3a$X97.5.[grep("phiA1",substr(rownames(ms3a),1,5))])
surv_rep_q97.5<-c(ms3a$X97.5.[grep("phiB0",substr(rownames(ms3a),1,5))],ms3a$X97.5.[grep("phiB1",substr(rownames(ms3a),1,5))])
trans_vr_q97.5<-c(ms3a$X97.5.[grep("psiA0",substr(rownames(ms3a),1,5))],ms3a$X97.5.[grep("psiA1",substr(rownames(ms3a),1,5))])
trans_rv_q97.5<-c(ms3a$X97.5.[grep("psiB0",substr(rownames(ms3a),1,5))],ms3a$X97.5.[grep("psiB1",substr(rownames(ms3a),1,5))])
emer_veg_q97.5<-c(ms3a$X97.5.[grep("pA0",substr(rownames(ms3a),1,3))],ms3a$X97.5.[grep("pA1",substr(rownames(ms3a),1,3))])
emer_rep_q97.5<-c(ms3a$X97.5.[grep("pB0",substr(rownames(ms3a),1,3))],ms3a$X97.5.[grep("pB1",substr(rownames(ms3a),1,3))])
#Figure 3, of vital rates
x<-c(1,2,1,2)
#x<-c(1,2,1.05,2.05)#jittered
xerror<-c(1,1,2,2)
#xerror<-c(1,1.05,2,2.05)#jittered
x2<-c(3,4,3,4)
#x2<-c(3,4,3.05,4.05)#jittered
x2error<-c(3,3,4,4)
#x2error<-c(3,3.05,4,4.05)#jittered
windows(height=6,width=7)
quartz(height=6,width=7)
par(mfrow=c(3,1),mar=c(.5,4.1,1,.5), oma=c(3,.6,.5,.5))
#survival
plot(x,c(surv_veg$control,surv_veg$logged), pch=21, bg=c("black","black","white","white"), ylim=c(0,1), ylab="Survival", xaxt="n", cex=1.5, xlab="", xlim=c(0.75,4.25), cex.lab=1.5, las=1,, cex.axis=1.3)
lines(x[1:2],c(surv_veg$control), lty=1)
lines(x[3:4],c(surv_veg$logged), lty=3)
abline(v=2.5,lty=1, lwd=2)
abline(v=1.5,lty=2,col="gray", lwd=2)
abline(v=3.5,lty=2,col="gray", lwd=2)
arrows(xerror,surv_veg_q2.5,xerror,surv_veg_q97.5, code=0,angle=90, length=0.1)
points(x,c(surv_veg$control,surv_veg$logged),pch=21, bg=c("black","black","white","white"), cex=1.5)
arrows(x2error,surv_rep_q2.5,x2error,surv_rep_q97.5, code=0,angle=90, length=0.1)
lines(x2[1:2],c(surv_rep$control), lty=1)
lines(x2[3:4],c(surv_rep$logged), lty=3)

points(x2,c(surv_rep$control,surv_rep$logged),pch=21, bg=c("black","black","white","white"), cex=1.5)
axis(side=1,at=c(1.5,3.5),labels=c("Vegetative","Reproductive" ),line=-15, tick=F, cex.axis=1.2)
legend("bottomright",legend=c("Control", "Cleared"),pch=21,pt.cex=1.5,pt.bg=c("black","white"), bty="n", cex=1.2)
#cbind(rownames(ms3a[25:32,]),ms3a$mean[25:32],ms3a$X2.5[25:32],ms3a$X97.5.[25:32])
#Dormancy(=1-)Emergence
plot(x,c(1-emer_veg$control,1-emer_veg$logged), pch=21, bg="black", ylim=c(0,1), ylab="Dormancy", xaxt="n", cex=1.5, xlab="", xlim=c(0.75,4.25), cex.lab=1.5, las=1, cex.axis=1.3)
lines(x[1:2],c(1-emer_veg$control), lty=1)
lines(x[3:4],c(1-emer_veg$logged), lty=3)
abline(v=2.5,lty=1, lwd=2)
abline(v=1.5,lty=2,col="gray", lwd=2)
abline(v=3.5,lty=2,col="gray", lwd=2)
arrows(xerror,1-emer_veg_q2.5,xerror,1-emer_veg_q97.5, code=0,angle=90, length=0.1)
points(x,c(1-emer_veg$control,1-emer_veg$logged),pch=21, bg=c("black","black","white","white"), cex=1.5)
arrows(x2error,1-emer_rep_q2.5,x2error,1-emer_rep_q97.5, code=0,angle=90, length=0.1)
lines(x2[1:2],c(1-emer_rep$control), lty=1)
lines(x2[3:4],c(1-emer_rep$logged), lty=3)
points(x2,c(1-emer_rep$control,1-emer_rep$logged),pch=21, bg=c("black","black","white","white"), cex=1.5)
cbind(rownames(ms3a[41:48,]),1-ms3a$mean[41:48],1-ms3a$X2.5[41:48],1-ms3a$X97.5.[41:48])#check error bars
#transition
plot(x,c(trans_vr$control,trans_vr$logged), pch=21, bg=c("black","black","white","white"), ylim=c(0,1), ylab="Transition", xaxt="n", cex=1.6, xlab="", xlim=c(0.75,4.25), cex.lab=1.5, las=1, cex.axis=1.3)
lines(x[1:2],c(trans_vr$control), lty=1)
lines(x[3:4],c(trans_vr$logged), lty=3)
abline(v=2.5,lty=1, lwd=2)
abline(v=1.5,lty=2,col="gray", lwd=2)
abline(v=3.5,lty=2,col="gray", lwd=2)
arrows(xerror,trans_vr_q2.5,xerror,trans_vr_q97.5, code=0,angle=90, length=0.1)
points(x,c(trans_vr$control,trans_vr$logged),pch=21, bg=c("black","black","white","white"), cex=1.5)
arrows(x2error,trans_rv_q2.5,x2error,trans_rv_q97.5, code=0,angle=90, length=0.1)
lines(x2[1:2],c(trans_rv$control), lty=1)
lines(x2[3:4],c(trans_rv$logged), lty=3)
points(x2,c(trans_rv$control,trans_rv$logged),pch=21, bg=c("black","black","white","white"), cex=1.5)
axis(side=1,at=c(x[1:2],x2[1:2]),labels=c("pre","post","pre","post"), cex.axis=1.3)
axis(side=1,at=c(x[1:2],x2[1:2]),labels=c("(1982-1997)","(1998-2015)","(1982-1997)","(1998-2015)"), line=1.2,tick=F, cex.axis=1.3)

#####Figure 4,  length of dormancy, lifepsan, etc

#x<-c(1,2,1,2)
windows(height=6,width=7)
quartz(height=6,width=7)
par(mfrow=c(3,1),mar=c(.5,4.1,1,.5), oma=c(3,.6,.5,.5))
#lifespan
plot(x,c(mean(lifespan_Xpre, na.rm=T),mean(lifespan_Xpost2, na.rm=T),mean(lifespan_Ypre, na.rm=T),mean(lifespan_Ypost2, na.rm=T)), pch=21, bg=c("black","black","white","white"), ylim=c(0,100), ylab="Lifespan (yrs)", xaxt="n", cex=1.5, xlab="", xlim=c(0.75,4.25), cex.lab=1.5, cex.axis=1.3, las=1)
abline(v=2.5,lty=1, lwd=2)
abline(v=1.5,lty=2,col="gray", lwd=2)
abline(v=3.5,lty=2,col="gray", lwd=2)
lines(x[1:2],c(mean(lifespan_Xpre, na.rm=T),mean(lifespan_Xpost, na.rm=T)), lty=1)
lines(x[3:4],c(mean(lifespan_Ypre, na.rm=T),mean(lifespan_Ypost2, na.rm=T)), lty=3)
arrows(xerror,c(mean(lifespan_Xpre, na.rm=T)-sd(lifespan_Xpre, na.rm=T),mean(lifespan_Ypre, na.rm=T)-sd(lifespan_Ypre, na.rm=T),mean(lifespan_Xpost2, na.rm=T)-sd(lifespan_Xpost2, na.rm=T),mean(lifespan_Ypost2, na.rm=T)-sd(lifespan_Ypost2, na.rm=T)),xerror,c(mean(lifespan_Xpre, na.rm=T)+sd(lifespan_Xpre, na.rm=T),mean(lifespan_Ypre, na.rm=T)+sd(lifespan_Ypre, na.rm=T),mean(lifespan_Xpost2, na.rm=T)+sd(lifespan_Xpost2, na.rm=T),mean(lifespan_Ypost2, na.rm=T)+sd(lifespan_Ypost2, na.rm=T)), code=0,angle=90, length=0.1)
points(x,c(mean(lifespan_Xpre, na.rm=T),mean(lifespan_Xpost2, na.rm=T),mean(lifespan_Ypre, na.rm=T),mean(lifespan_Ypost2, na.rm=T)), pch=21, bg=c("black","black","white","white"), cex=1.5)
lines(x2[1:2],c(mean(lifespan_flow_Xpre, na.rm=T),mean(lifespan_flow_Xpost2, na.rm=T)), lty=1)
lines(x2[3:4],c(mean(lifespan_flow_Ypre, na.rm=T),mean(lifespan_flow_Ypost2, na.rm=T)), lty=3)
arrows(x2error,c(mean(lifespan_flow_Xpre, na.rm=T)-sd(lifespan_flow_Xpre, na.rm=T),mean(lifespan_flow_Ypre, na.rm=T)-sd(lifespan_flow_Ypre, na.rm=T),mean(lifespan_flow_Xpost2, na.rm=T)-sd(lifespan_flow_Xpost2, na.rm=T),mean(lifespan_flow_Ypost2, na.rm=T)-sd(lifespan_flow_Ypost2, na.rm=T)),x2error,c(mean(lifespan_flow_Xpre, na.rm=T)+sd(lifespan_flow_Xpre, na.rm=T),mean(lifespan_flow_Ypre, na.rm=T)+sd(lifespan_flow_Ypre, na.rm=T),mean(lifespan_flow_Xpost2, na.rm=T)+sd(lifespan_flow_Xpost2, na.rm=T),mean(lifespan_flow_Ypost2, na.rm=T)+sd(lifespan_flow_Ypost2, na.rm=T)), code=0,angle=90, length=0.1)
points(x2,c(mean(lifespan_flow_Xpre, na.rm=T),mean(lifespan_flow_Xpost2, na.rm=T),mean(lifespan_flow_Ypre, na.rm=T),mean(lifespan_flow_Ypost2, na.rm=T)), pch=21, bg=c("black","black","white","white"), cex=1.5)
axis(side=1,at=c(1.5,3.5),labels=c("Vegetative","Reproductive" ),line=-15, tick=F, cex.axis=1.2)

#Length of dormancy, starting from veg (black) or rep (white)
plot(x,c(mean(lengthdorm_Xpre, na.rm=T),mean(lengthdorm_Xpost, na.rm=T),mean(lengthdorm_Ypre, na.rm=T),mean(lengthdorm_Ypost, na.rm=T)), pch=21, bg=c("black","black","white","white"), ylim=c(0,2), ylab="Dormancy length (yrs)", xaxt="n", cex=1.5, xlab="", xlim=c(0.75,4.25), cex.lab=1.2,cex.lab=1.5, cex.axis=1.3, las=1)
abline(v=2.5,lty=1, lwd=2)
abline(v=1.5,lty=2,col="gray", lwd=2)
abline(v=3.5,lty=2,col="gray", lwd=2)
lines(x[1:2],c(mean(lengthdorm_Xpre, na.rm=T),mean(lengthdorm_Xpost, na.rm=T)), lty=1)
lines(x[3:4],c(mean(lengthdorm_Ypre, na.rm=T),mean(lengthdorm_Ypost, na.rm=T)), lty=3)
arrows(xerror,c(mean(lengthdorm_Xpre, na.rm=T)-sd(lengthdorm_Xpre, na.rm=T),mean(lengthdorm_Ypre, na.rm=T)-sd(lengthdorm_Ypre, na.rm=T),mean(lengthdorm_Xpost, na.rm=T)-sd(lengthdorm_Xpost, na.rm=T),mean(lengthdorm_Ypost, na.rm=T)-sd(lengthdorm_Ypost, na.rm=T)),xerror,c(mean(lengthdorm_Xpre, na.rm=T)+sd(lengthdorm_Xpre, na.rm=T),mean(lengthdorm_Ypre, na.rm=T)+sd(lengthdorm_Ypre, na.rm=T),mean(lengthdorm_Xpost, na.rm=T)+sd(lengthdorm_Xpost, na.rm=T),mean(lengthdorm_Ypost, na.rm=T)+sd(lengthdorm_Ypost, na.rm=T)), code=0,angle=90, length=0.1)
points(x,c(mean(lengthdorm_Xpre, na.rm=T),mean(lengthdorm_Xpost, na.rm=T),mean(lengthdorm_Ypre, na.rm=T),mean(lengthdorm_Ypost, na.rm=T)), pch=21, bg=c("black","black","white","white"), cex=1.5)
arrows(x2error,c(mean(lengthdorm_flow_Xpre, na.rm=T)-sd(lengthdorm_flow_Xpre, na.rm=T),mean(lengthdorm_flow_Ypre, na.rm=T)-sd(lengthdorm_flow_Ypre, na.rm=T),mean(lengthdorm_flow_Xpost, na.rm=T)-sd(lengthdorm_flow_Xpost, na.rm=T),mean(lengthdorm_flow_Ypost, na.rm=T)-sd(lengthdorm_flow_Ypost, na.rm=T)),x2error,c(mean(lengthdorm_flow_Xpre, na.rm=T)+sd(lengthdorm_flow_Xpre, na.rm=T),mean(lengthdorm_flow_Ypre, na.rm=T)+sd(lengthdorm_flow_Ypre, na.rm=T),mean(lengthdorm_flow_Xpost, na.rm=T)+sd(lengthdorm_flow_Xpost, na.rm=T),mean(lengthdorm_flow_Ypost, na.rm=T)+sd(lengthdorm_flow_Ypost, na.rm=T)), code=0,angle=90, length=0.1)
lines(x2[1:2],c(mean(lengthdorm_flow_Xpre, na.rm=T),mean(lengthdorm_flow_Xpost, na.rm=T)), lty=1)
lines(x2[3:4],c(mean(lengthdorm_flow_Ypre, na.rm=T),mean(lengthdorm_flow_Ypost, na.rm=T)), lty=3)
points(x2,c(mean(lengthdorm_flow_Xpre, na.rm=T),mean(lengthdorm_flow_Xpost, na.rm=T),mean(lengthdorm_flow_Ypre, na.rm=T),mean(lengthdorm_flow_Ypost, na.rm=T)), pch=21, bg=c("black","black","white","white"), cex=1.5)
axis(side=1,at=c(x[1:2],x2[1:2]),labels=c("pre","post","pre","post"), cex.axis=1.3)
axis(side=1,at=c(x[1:2],x2[1:2]),labels=c("(1982-1997)","(1998-2015)","(1982-1997)","(1998-2015)"), line=1.2,tick=F, cex.axis=1.3)

#proportion of plants dormant
plot(x,c(mean(propdorm_Xpre, na.rm=T),mean(propdorm_Xpost, na.rm=T),mean(propdorm_Ypre, na.rm=T),mean(propdorm_Ypost, na.rm=T)), pch=21, bg=c("black","black","white","white"), ylim=c(0,1), ylab="Proportion Dormant", xaxt="n", cex=1.5, xlab="", xlim=c(0.75,4.25), cex.lab=1.5, cex.axis=1.3, las=1)
abline(v=2.5,lty=1, lwd=2)
abline(v=1.5,lty=2,col="gray", lwd=2)
lines(x[1:2],c(mean(propdorm_Xpre, na.rm=T),mean(propdorm_Xpost, na.rm=T)), lty=1)
lines(x[3:4],c(mean(propdorm_Ypre, na.rm=T),mean(propdorm_Ypost, na.rm=T)), lty=3)
arrows(xerror,c(mean(propdorm_Xpre, na.rm=T)-sd(propdorm_Xpre, na.rm=T),mean(propdorm_Ypre, na.rm=T)-sd(propdorm_Ypre, na.rm=T),mean(propdorm_Xpost, na.rm=T)-sd(propdorm_Xpost, na.rm=T),mean(propdorm_Ypost, na.rm=T)-sd(propdorm_Ypost, na.rm=T)),xerror,c(mean(propdorm_Xpre, na.rm=T)+sd(propdorm_Xpre, na.rm=T),mean(propdorm_Ypre, na.rm=T)+sd(propdorm_Ypre, na.rm=T),mean(propdorm_Xpost, na.rm=T)+sd(propdorm_Xpost, na.rm=T),mean(propdorm_Ypost, na.rm=T)+sd(propdorm_Ypost, na.rm=T)), code=0,angle=90, length=0.1)
points(x,c(mean(propdorm_Xpre, na.rm=T),mean(propdorm_Xpost, na.rm=T),mean(propdorm_Ypre, na.rm=T),mean(propdorm_Ypost, na.rm=T)), pch=21, bg=c("black","black","white","white"), cex=1.5)
axis(side=1,at=c(x[1:2]),labels=c("pre","post"), cex.axis=1.3)
axis(side=1,at=c(x[1:2]),labels=c("(1982-1997)","(1998-2015)"), line=1.2,tick=F, cex.axis=1.3)


##Figuring out why error bars are so wide for lifespan Ypost
inf.params<-cbind(phiV_Ypost[which(lifespan_Ypost=="Inf")],veg.rep_Ypost[which(lifespan_Ypost=="Inf")],pdormV_Ypost[which(lifespan_Ypost=="Inf")],phiR_Ypost[which(lifespan_Ypost=="Inf")],rep.veg_Ypost[which(lifespan_Ypost=="Inf")],pdormR_Ypost[which(lifespan_Ypost=="Inf")])
noninf.params<-cbind(phiV_Ypost[which(lifespan_Ypost!="Inf")],veg.rep_Ypost[which(lifespan_Ypost!="Inf")],pdormV_Ypost[which(lifespan_Ypost!="Inf")],phiR_Ypost[which(lifespan_Ypost!="Inf")],rep.veg_Ypost[which(lifespan_Ypost!="Inf")],pdormR_Ypost[which(lifespan_Ypost!="Inf")])
t.test(noninf.params[,6],inf.params[,6])
high.params<-cbind(phiV_Ypost[which(lifespan_Ypost>200)],veg.rep_Ypost[which(lifespan_Ypost>200)],pdormV_Ypost[which(lifespan_Ypost>200)],phiR_Ypost[which(lifespan_Ypost>200)],rep.veg_Ypost[which(lifespan_Ypost>200)],pdormR_Ypost[which(lifespan_Ypost>200)])
low.params<-cbind(phiV_Ypost[which(lifespan_Ypost<200)],veg.rep_Ypost[which(lifespan_Ypost<200)],pdormV_Ypost[which(lifespan_Ypost<200)],phiR_Ypost[which(lifespan_Ypost<200)],rep.veg_Ypost[which(lifespan_Ypost<200)],pdormR_Ypost[which(lifespan_Ypost<200)])
lowlow.params<-cbind(phiV_Ypost[which(lifespan_Ypost<10)],veg.rep_Ypost[which(lifespan_Ypost<10)],pdormV_Ypost[which(lifespan_Ypost<10)],phiR_Ypost[which(lifespan_Ypost<10)],rep.veg_Ypost[which(lifespan_Ypost<10)],pdormR_Ypost[which(lifespan_Ypost<10)])
t.test(lowlow.params[,5],low.params[,5])
#vital rates from high lifespan estimates (>200 years) have the following differences from low lifespan estimates:
#1) higher phis for both reproductive and veg plants
#2) lower transition from reproductive to vegetative rep.veg_Ypost
#vital rates frmo low low lifespand estimates (<10 years) have the following
#1) higher phis for both reproductive and veg plants 
#2) higher transition from reproductive to vegetative rep.veg_Ypost
