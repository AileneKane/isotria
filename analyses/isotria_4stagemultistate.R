#Multistate model for Isotria medioloides Alton NH population
#Data provided by Bill Brumback
#Coding by Ailene Ettinger with help frmo Andy Royle and Elizabeth Crone
#this file has code for a 4 state multistate model to analyze two groups (x and y) from a population of Isotria medeoloides
#the model here allows vital rate estimates in the multistate model to vary independently among years in their random effect structure
#started September 26, 2016
#setwd("~/isotria") #at usgs
setwd("~/git/isotria/analyses")
rm(list=ls()) 
options(stringsAsFactors=FALSE)

library(rjags)
library(jagsUI)
library(lattice)
library(coda)
library(boot)
#update.packages()
isoinds<-read.csv("Isotria_Stage3_2016.csv", header=T)#
isoinds<-isoinds[-which(isoinds$UniqueID=="X-01-244"),]#individual has to be removed because of errors in its monitoring
#Add column for emergent/nonemergent
isoinds$Emerg<-NA
isoinds[which(isoinds$TotNoStems>0),]$Emerg=1
isoinds[which(isoinds$TotNoStems==0),]$Emerg=0
##Select out just groups X and Y for this analysis
isoindsX=isoinds[isoinds$Group=="X",]
isoindsY=isoinds[isoinds$Group=="Y",]
isoindsXY=isoinds[isoinds$Group=="X"|isoinds$Group=="Y",]
isoindsXY$UniqueID=factor(isoindsXY$UniqueID)
isoindsXY$Group=factor(isoindsXY$Group)
##get isotria data into format such that 1=veg, 2=rep and 3=not seen
#to do this, Add column for reproductive (=arrested, flowering, or fruiting)/not rep
isoindsXY$Repro<-NA
isoindsXY[which(isoindsXY$NoFrStems>0|isoindsXY$NoFlStems>0|isoindsXY$NoArrStems>0),]$Repro=1
isoindsXY[which(isoindsXY$NoFrStems==0&isoindsXY$NoFlStems==0&isoindsXY$NoArrStems==0),]$Repro=0
isoindsXY<-isoindsXY[-which(isoindsXY$Year<isoindsXY$YearFound),]
# need to use previous state to get 2 dormant stages
past = isoindsXY[isoindsXY$Year < 2015,]; past = past[order(past$UniqueID),]; past = past[order(past$Year),]
dim(past)# 4315 (dim(isoinds) is 4570 rows)
past$YearPrev=past$Year+1
head(past)
colnames(past)[5]<-"YearCurrent"
#Merge present and past to get fates
isoinds_prev<-merge(isoindsXY,past,by.x=c("UniqueID","Year","Group",  "Block", "Plant","YearFound"),by.y=c("UniqueID","YearPrev","Group",  "Block", "Plant","YearFound"),all.x = TRUE)#this includes all plants, without matches in past, as well
#delete unncessary columns (we only want stage)
isoinds_prev2=subset(isoinds_prev,select=c("UniqueID","Year", "Group",  "Block", "Plant",  "YearFound","Emerg.x","Repro.x","Emerg.y","Repro.y","Stage.x","Stage.y"))
#I want to add a column with a 1 in it if the plant is dormant AND the last observed state was reproductive
isoinds_prev2$dormstage=NA
isoinds_prev2[which(isoinds_prev2$Stage.y=="Veg" & isoinds_prev2$Stage.x=="D"),]$dormstage="UV"
isoinds_prev2[which(isoinds_prev2$Stage.y=="Flower" & isoinds_prev2$Stage.x=="D"),]$dormstage="UF"
isoinds_prev2[which(isoinds_prev2$Stage.y=="Arr" & isoinds_prev2$Stage.x=="D"),]$dormstage="UF"
isoinds_prev2[which(isoinds_prev2$Stage.y=="Fruit" & isoinds_prev2$Stage.x=="D"),]$dormstage="UF"
for(i in 2:dim(isoinds_prev2)[1]){
  ind<-isoinds_prev2[i,which(colnames(isoinds_prev2)=="UniqueID")]
  currentstage<-isoinds_prev2[i,which(colnames(isoinds_prev2)=="Stage.x")]
  previousdormstage<-isoinds_prev2[i-1,which(colnames(isoinds_prev2)=="dormstage")]
  currentdormstage<-isoinds_prev2[i,which(colnames(isoinds_prev2)=="dormstage")]
  if(currentstage!="D"){next}
  if(currentstage=="D"&!is.na(currentdormstage)){next}
  if(isoinds_prev2$UniqueID[i]==ind & currentstage=="D" & previousdormstage=="UV"){isoinds_prev2[i,which(colnames(isoinds_prev2)=="dormstage")]<-"UV"}
  else if (isoinds_prev2$UniqueID[i]==ind & currentstage=="D" & previousdormstage=="UF"){isoinds_prev2[i,which(colnames(isoinds_prev2)=="dormstage")]<-"UF"}
  }
#Make new column with a 1 for when dormant stage=UF
isoinds_prev2$dormstageUF<-0
isoinds_prev2[which(isoinds_prev2$dormstage=="UF"),]$dormstageUF<-1
isoalleme_CH=tapply(isoinds_prev2$Emerg.x,list(isoinds_prev2$UniqueID,isoinds_prev2$Year),sum)#emergence
isoallrep_CH=tapply(isoinds_prev2$Repro.x,list(isoinds_prev2$UniqueID,isoinds_prev2$Year),sum)#reproductive status
isoall_CH=isoalleme_CH+isoallrep_CH
isoall_CH[which(isoall_CH==0)]=3#replace 0s wih 3
#isoall_CH.ms[73,23]<-NA#remove case where first observation is a zero/3- makes no sense
#add 4s for when dormant stage =UF
isoallUF_CH=tapply(isoinds_prev2$dormstageUF,list(isoinds_prev2$UniqueID,isoinds_prev2$Year),sum)#reproductive status
isoall_CH.ms=isoallUF_CH+isoall_CH
isoall_CH.ms[which(isoall_CH.ms==4)]<-3#BAsed on what Nate said, this should not be explicitly coded this way because we can't see this stage
n.occasions<- dim(isoall_CH.ms)[2]
#head(isoall_CH.ms)
get.first <- function(x) min(which(x!=0))
f<-apply(isoall_CH.ms,1,get.first)#first occasion of marking
###################################
####Add variable with difference among groups x and y before and after logging occurring
group<-c(rep(1,length(unique(isoinds[isoinds$Group=="X",]$UniqueID))),rep(2,length(unique(isoinds[isoinds$Group=="Y",]$UniqueID))))#group x=control=1, group y=logged=2
logged_x<-c(rep(0,30))#logging never occurred for group x- this is one less than the number of years 
logged_y<-c(rep(0,12),rep(1,18))#for group y, logging occurred between 1997 and 1998, so first 13 years= no logging treatment
logged_yrs<-c(rep(0,13),seq(1:17))#add decay time- years since logging this is Beta1, this is x1, add a year after the initial year of logging to allow for a bump up
logged_yrs2<-rbind(logged_yrs,logged_yrs)
time_log<-rbind(logged_y,logged_y)#to test differene between groups before and after logging

#Fit a multistate model in JAGS, with random effects of year####
#The random effect allows time-variance for all vital rates (i.e. random effect of time), plus fixed effect of group (X vs Y)
sink("ms-ranef4stages.jags")
cat("

model {
    # -------------------------------------------------
    # Parameters:
    # sV: survival probability for nonreproductive plants
    # sF: survival probability for reproductive plants
    # sUV: survival probability for dormant plants that were vegetative above ground
    # sUF: survival probability for dormant plants that were reproductive above ground
    # fV: reproduction probability for vegetative plants
    # fF: reproduction probability for reproductive plants
    # fUV: reproduction probability for dormant plants that were vegetative above ground
    # fUF: reproduction probability for dormant plants that were reproductive above ground
    # dV: dormancy probability for nonreproductive plants
    # dF: dormancy probability for reproductive plants
    # dUV: dormancy probability for dormant plants that were vegetative above ground
    # dUF: dormancy probability for dormant plants that were reproductive above ground
    # -------------------------------------------------
    # States (S):
    # 1 vegetative
    # 2 reproductive
    # 3 dormant, previously vegetative
    # 4 dormant, previously reproductive
    # 5 dead
    # Observations (O):  
    # 1 seen vegetative
    # 2 seen reproductive
    # 3 not seen
    # -------------------------------------------------
    
    # Model, priors and constraints
    for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
    logit(sV[i,t]) <- eta.sV[group[i],t]
    logit(sF[i,t]) <- eta.sF[group[i],t]
    logit(sUV[i,t]) <- eta.sUV[group[i],t]
    logit(sUF[i,t]) <- eta.sUF[group[i],t]
    logit(fV[i,t]) <- eta.fV[group[i],t]
    logit(fF[i,t]) <- eta.fF[group[i],t]
    logit(fUV[i,t]) <- eta.fUV[group[i],t]
    logit(fUF[i,t]) <- eta.fUF[group[i],t]
    logit(dV[i,t]) <- eta.dV[group[i],t]
    logit(dF[i,t]) <- eta.dF[group[i],t]
    logit(dUV[i,t]) <- eta.dUV[group[i],t]
    logit(dUF[i,t]) <- eta.dUF[group[i],t]
    }#t
    }#i
    for (g in 1:2){#group  
    for (t in 1:(n.occasions-1)){
    eta.sV[g,t]  <-mu.sV[g]+beta.sV[g]*x[g,t]+beta1.sV[g]*x1[g,t]+epsilon.sV[g,t]
    eta.sF[g,t]  <-mu.sF[g]+beta.sF[g]*x[g,t]+beta1.sF[g]*x1[g,t]+epsilon.sF[g,t]
    eta.sUV[g,t]  <-mu.sUV[g]+beta.sUV[g]*x[g,t]+beta1.sUV[g]*x1[g,t]+epsilon.sUV[g,t]
    eta.sUF[g,t]  <-mu.sUF[g]+beta.sUF[g]*x[g,t]+beta1.sUF[g]*x1[g,t]+epsilon.sUF[g,t]
    eta.fV[g,t]  <-mu.fV[g]+beta.fV[g]*x[g,t]+beta1.fV[g]*x1[g,t]+epsilon.fV[g,t]
    eta.fF[g,t]  <-mu.fF[g]+beta.fF[g]*x[g,t]+beta1.fF[g]*x1[g,t]+epsilon.fF[g,t]
    eta.fUV[g,t]  <-mu.fUV[g]+beta.fUV[g]*x[g,t]+beta1.fUV[g]*x1[g,t]+epsilon.fUV[g,t]
    eta.fUF[g,t]  <-mu.fUF[g]+beta.fUF[g]*x[g,t]+beta1.fUF[g]*x1[g,t]+epsilon.fUF[g,t]
    eta.dV[g,t]  <-mu.dV[g]+beta.dV[g]*x[g,t]+beta1.dV[g]*x1[g,t]+epsilon.dV[g,t]
    eta.dF[g,t]  <-mu.dF[g]+beta.dF[g]*x[g,t]+beta1.dF[g]*x1[g,t]+epsilon.dF[g,t]
    eta.dUV[g,t]  <-mu.dUV[g]+beta.dUV[g]*x[g,t]+beta1.dUV[g]*x1[g,t]+epsilon.dUV[g,t]
    eta.dUF[g,t]  <-mu.dUF[g]+beta.dUF[g]*x[g,t]+beta1.dUF[g]*x1[g,t]+epsilon.dUF[g,t]
    epsilon.sV[g,t]~dnorm(0,tau.sV[g])#could move all the mean stuff to this part, might help model mix better? hierarchical centering
    epsilon.sF[g,t]~dnorm(0,tau.sF[g])
    epsilon.sUV[g,t]~dnorm(0,tau.sUV[g])
    epsilon.sUF[g,t]~dnorm(0,tau.sUF[g])
    epsilon.fV[g,t]~dnorm(0,tau.fV[g])
    epsilon.fF[g,t]~dnorm(0,tau.fF[g])
    epsilon.fUV[g,t]~dnorm(0,tau.fUV[g])
    epsilon.fUF[g,t]~dnorm(0,tau.fUF[g])
    epsilon.dV[g,t]~dnorm(0,tau.dV[g])
    epsilon.dF[g,t]~dnorm(0,tau.dF[g])
    epsilon.dUV[g,t]~dnorm(0,tau.dUV[g])
    epsilon.dUF[g,t]~dnorm(0,tau.dUF[g])
      }#t
    mean.sV[g]~dunif(0,1)# Priors for mean group-specific survival for veg plants
    mu.sV[g]<-log(mean.sV[g]/(1-mean.sV[g]))  
    mean.sF[g]~dunif(0,1)# Priors for mean group-specific survivalsurvival for rep plants
    mu.sF[g]<-log(mean.sF[g]/(1-mean.sF[g]))    
    mean.sUV[g]~dunif(0,1)# Priors for mean group-specific survival for dorm(prev veg) plants
    mu.sUV[g]<-log(mean.sUV[g]/(1-mean.sUV[g])) 
    mean.sUF[g]~dunif(0,1)# Priors for mean group-specific survival dorm (prev rep) plants
    mu.sUF[g]<-log(mean.sUF[g]/(1-mean.sUF[g])) 
    mean.fV[g]~dunif(0,1)# Priors for mean group-specific prob of rep for veg plants
    mu.fV[g]<-log(mean.fV[g]/(1-mean.fV[g]))       
    mean.fF[g]~dunif(0,1)# Priors for mean group-specific prob of rep for rep plants
    mu.fF[g]<-log(mean.fF[g]/(1-mean.fF[g]))    
    mean.fUV[g]~dunif(0,1)# Priors for mean group-specific prob of rep for dorm (prev veg) plants
    mu.fUV[g]<-log(mean.fUV[g]/(1-mean.fUV[g])) 
    mean.fUF[g]~dunif(0,1)# Priors for mean group-specific prob of rep for dorm (prev rep) plants
    mu.fUF[g]<-log(mean.fUF[g]/(1-mean.fUF[g])) 
    mean.dV[g]~dunif(0,1)# Priors for mean group-specific dormancy for veg plants
    mu.dV[g]<-log(mean.dV[g]/(1-mean.dV[g]))       
    mean.dF[g]~dunif(0,1)# Priors for mean group-specific dormancy for rep plants
    mu.dF[g]<-log(mean.dF[g]/(1-mean.dF[g]))    
    mean.dUV[g]~dunif(0,1)# Priors for mean group-specific dormancy for dorm (prev veg) plants
    mu.dUV[g]<-log(mean.dUV[g]/(1-mean.dUV[g]))
    mean.dUF[g]~dunif(0,1)# Priors for mean group-specific dormancy for dorm (prev rep) plants
    mu.dUF[g]<-log(mean.dUF[g]/(1-mean.dUV[g])) 
    sigma.sV[g]~dunif(0,10)#temporal variance for veg plants survival
    tau.sV[g]<-pow(sigma.sV[g],-2)
    sigma.sV2[g]<-pow(sigma.sV[g],2)
    sigma.sF[g]~dunif(0,10)#temporal variance for rep plants survival
    tau.sF[g]<-pow(sigma.sF[g],-2)
    sigma.sF2[g]<-pow(sigma.sF[g],2)
    sigma.sUV[g]~dunif(0,10)#temporal variance for dorm plants survival
    tau.sUV[g]<-pow(sigma.sUV[g],-2)
    sigma.sUV2[g]<-pow(sigma.sUV[g],2)
    sigma.sUF[g]~dunif(0,10)#temporal variance for dorm plants survival
    tau.sUF[g]<-pow(sigma.sUF[g],-2)
    sigma.sUF2[g]<-pow(sigma.sUF[g],2)
    sigma.fV[g]~dunif(0,10)#temporal variance for prob of rep
    tau.fV[g]<-pow(sigma.fV[g],-2)
    sigma.fV2[g]<-pow(sigma.fV[g],2)
    sigma.fF[g]~dunif(0,10)#
    tau.fF[g]<-pow(sigma.fF[g],-2)
    sigma.fF2[g]<-pow(sigma.fF[g],2)
    sigma.fUV[g]~dunif(0,10)
    tau.fUV[g]<-pow(sigma.fUV[g],-2)
    sigma.fUV2[g]<-pow(sigma.fUV[g],2)
    sigma.fUF[g]~dunif(0,10)
    tau.fUF[g]<-pow(sigma.fUF[g],-2)
    sigma.fUF2[g]<-pow(sigma.fUF[g],2)
    sigma.dV[g]~dunif(0,10)#temporal variance of dormancy for veg plants 
    tau.dV[g]<-pow(sigma.dV[g],-2)
    sigma.dV2[g]<-pow(sigma.dV[g],2)
    sigma.dF[g]~dunif(0,10)#temporal variance of dormancy for rep plants
    tau.dF[g]<-pow(sigma.dF[g],-2)
    sigma.dF2[g]<-pow(sigma.dF[g],2)
    sigma.dUV[g]~dunif(0,10)#temporal variance of dormancy for dormant plants
    tau.dUV[g]<-pow(sigma.dUV[g],-2)
    sigma.dUV2[g]<-pow(sigma.dUV[g],2)
     sigma.dUF[g]~dunif(0,10)#temporal variance of dormancy for dormant plants
    tau.dUF[g]<-pow(sigma.dUF[g],-2)
    sigma.dUF2[g]<-pow(sigma.dUF[g],2)
    beta.sV[g]~dnorm(0,0.001)I(-10,10)
    beta.sF[g]~dnorm(0,0.001)I(-10,10)
    beta.sUV[g]~dnorm(0,0.001)I(-10,10)
    beta.sUF[g]~dnorm(0,0.001)I(-10,10)
    beta.fV[g]~dnorm(0,0.001)I(-10,10)
    beta.fF[g]~dnorm(0,0.001)I(-10,10)
    beta.fUV[g]~dnorm(0,0.001)I(-10,10)
    beta.fUF[g]~dnorm(0,0.001)I(-10,10)
    beta.dV[g]~dnorm(0,0.001)I(-10,10)
    beta.dF[g]~dnorm(0,0.001)I(-10,10)
    beta.dUV[g]~dnorm(0,0.001)I(-10,10)
    beta.dUF[g]~dnorm(0,0.001)I(-10,10)
    beta1.sV[g]~dnorm(0,0.001)I(-10,10)
    beta1.sF[g]~dnorm(0,0.001)I(-10,10)
    beta1.sUV[g]~dnorm(0,0.001)I(-10,10)
    beta1.sUF[g]~dnorm(0,0.001)I(-10,10)
    beta1.fV[g]~dnorm(0,0.001)I(-10,10)
    beta1.fF[g]~dnorm(0,0.001)I(-10,10)
    beta1.fUV[g]~dnorm(0,0.001)I(-10,10)
    beta1.fUF[g]~dnorm(0,0.001)I(-10,10)
    beta1.dV[g]~dnorm(0,0.001)I(-10,10)
    beta1.dF[g]~dnorm(0,0.001)I(-10,10)
    beta1.dUV[g]~dnorm(0,0.001)I(-10,10)
    beta1.dUF[g]~dnorm(0,0.001)I(-10,10)
    #calculate probabilities of 2 different groups and 2 different time periods (=4 group-time periods) to look at later... 
    #survival
    logit(sV0[g])<- mu.sV[g]
    logit(sV1[g])<- mu.sV[g] + beta.sV[g]
    logit(sF0[g])<- mu.sF[g]
    logit(sF1[g])<- mu.sF[g] + beta.sF[g]
    logit(sUV0[g])<- mu.sUV[g]
    logit(sUV1[g])<- mu.sUV[g] + beta.sUV[g]
    logit(sUF0[g])<- mu.sUF[g]
    logit(sUF1[g])<- mu.sUF[g] + beta.sUF[g]
    #reproduction
    logit(fV0[g])<- mu.fV[g]
    logit(fV1[g])<- mu.fV[g] + beta.fV[g]
    logit(fF0[g])<- mu.fF[g]
    logit(fF1[g])<- mu.fF[g] + beta.fF[g]
    logit(fUV0[g])<- mu.fUV[g]
    logit(fUV1[g])<- mu.fUV[g] + beta.fUV[g]
    logit(fUF0[g])<- mu.fUF[g]
    logit(fUF1[g])<- mu.fUF[g] + beta.fUF[g]
    #dormancy
    logit(dV0[g])<- mu.dV[g]
    logit(dV1[g])<- mu.dV[g] + beta.dV[g]
    logit(dF0[g])<- mu.dF[g]
    logit(dF1[g])<- mu.dF[g] + beta.dF[g]
    logit(dUV0[g])<- mu.dUV[g]
    logit(dUV1[g])<- mu.dUV[g] + beta.dUV[g]
    logit(dUF0[g])<- mu.dUF[g]
    logit(dUF1[g])<- mu.dUF[g] + beta.dUF[g]
    }#g
    
    # Define state-transition and observation matrices
    for (i in 1:nind){  
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    ps[1,i,t,1] <- sV[i,t] * (1-fV[i,t])*(1-dV[i,t])
    ps[1,i,t,2] <- sV[i,t] * fV[i,t]*(1-dV[i,t])
    ps[1,i,t,3] <- sV[i,t]*dV[i,t]
    ps[1,i,t,4] <- 0
    ps[1,i,t,5] <- 1-sV[i,t]
    ps[2,i,t,1] <- sF[i,t] * (1-fF[i,t])*(1-dF[i,t])
    ps[2,i,t,2] <- sF[i,t] * fF[i,t]*(1-dF[i,t])
    ps[2,i,t,3] <- 0
    ps[2,i,t,4] <- sF[i,t]*dF[i,t]
    ps[2,i,t,5] <- 1-sF[i,t]
    ps[3,i,t,1] <- sUV[i,t] * (1-fUV[i,t])*(1-dUV[i,t])
    ps[3,i,t,2] <- sUV[i,t] * fUV[i,t]*(1-dUV[i,t])
    ps[3,i,t,3] <- sUV[i,t]* dUV[i,t]
    ps[3,i,t,4] <- 0
    ps[3,i,t,5] <- 1-sUV[i,t]
    ps[4,i,t,1] <- sUF[i,t] * (1-fUF[i,t])*(1-dUF[i,t])
    ps[4,i,t,2] <- sUF[i,t] * fUF[i,t]*(1-dUF[i,t])
    ps[4,i,t,3] <- 0
    ps[4,i,t,4] <- sUF[i,t]* dUF[i,t]
    ps[4,i,t,5] <- 1-sUF[i,t]
    ps[5,i,t,1] <- 0
    ps[5,i,t,2] <- 0
    ps[5,i,t,3] <- 0
    ps[5,i,t,4] <- 0
    ps[5,i,t,5] <- 1
    # Define probabilities of O(t) given S(t) ##first column is observed state, last is true state
    po[1,i,t,1] <- 1
    po[1,i,t,2] <- 0
    po[1,i,t,3] <-0
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- 1
    po[2,i,t,3] <- 0
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- 1
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 1
    po[5,i,t,1] <- 0
    po[5,i,t,2] <- 0
    po[5,i,t,3] <- 1

    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

#Function to create known latent state z
########
known.state.ms <- function(ch){##removes 3s and replaces them with NAs, and replaces first observation with NA
  state <- ch
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,]<3))#
    state[i,n1] <- NA
  }
  state[state==3] <- NA
  #state[state==4] <- NA
  return(state)
}
##I need to recode this so that the starting value is always a 2 for unobserved instances when the last seen state was 2 and a 1 for unobserved instances when the last seen state was 1
ms.init.z <- function(ch, f){#ms.init.z gives starting values of 1 or 2 to all unknown states
  for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
  n1 <- min(which(ch[i,]<3))#
  states <- max(ch, na.rm = TRUE)
  known.states <- 1:2
  v <- which(ch>=3)#which occurences are unknown states (=3 and 4)
  ch[-v] <- NA#everything else besides the 3 gets an NA
  ch[v] <- sample(known.states, length(v), replace = TRUE)##gives all unseen occurrences (3s and 4s) a starting value that is a random sampling of either 1 or 2
  ch[i,n1] <- NA#make first observance an NA (not a 3)
  return(ch)
}

ch=isoall_CH.ms
#y = isoall_CH.ms
#y[73,24]=NA
zst1=ms.init.z(ch,f)
#fixes that i tried but did not help the model run!
zst1[119,21:31]<-NA#a work around, until i fix the starting value code!
zst1[4,24:31]<-NA
zst1[6,24:31]<-NA
zst1[5,24:31]<-NA
zst1[60,31]<-NA
#zst1[255,10:31]<-NA#because of "cannot normalize density" error


# Bundle data
jags.data <- list(y = isoall_CH.ms, f = f, n.occasions = dim(isoall_CH.ms)[2], nind = dim(isoall_CH.ms)[1], z = known.state.ms(isoall_CH.ms), group=group, x=time_log, x1=logged_yrs2)

# Initial values
inits<-function(){list(mean.sV=c(.5,.5), mean.sF=c(.5,.5),mean.sUV=c(.5,.5),mean.sUF=c(.5,.5),mean.fV=c(.5,.5),mean.fF=c(.5,.5),mean.fUV=c(.5,.5),mean.fUF=c(.5,.5),mean.dV=c(.5,.5),mean.dF=c(.5,.5),mean.dUV=c(.5,.5),mean.dUF=c(.5,.5),sigma.sV=c(1,1),sigma.sF=c(1,1),sigma.sUV=c(1,1),sigma.sUF=c(1,1),sigma.fV=c(1,1),sigma.fF=c(1,1),sigma.fUV=c(1,1),sigma.fUF=c(1,1),sigma.dV=c(1,1),sigma.dF=c(1,1),sigma.dUV=c(1,1),sigma.dUF=c(1,1),z = zst1)}
# Parameters monitored
parameters <- c("mean.sV","mean.sF","mean.sUV","mean.sUF", "mean.fV","mean.fF","mean.fUV","mean.fUF","mean.dV","mean.dF","mean.dUV","mean.dUF","beta.sV","beta.sF", "beta.sUV","beta.sUF","beta.fV","beta.fF","beta.fUV","beta.fUF","beta.dV", "beta.dF", "beta.dUV","beta.dUF","sV0","sV1","sF0","sF1","sUV0","sUF0","sUV1","sUF1","fV0","fV1","fF0","fF1","fUV0","fUV1","fUF0","fUF1","dV0","dV1","dF0","dF1","dUV0","dUF0","dUV1","dUF1","beta1.sV","beta1.sF","beta1.sUV","beta1.sUF","beta1.fV","beta1.fF","beta1.fUV","beta1.fUF","beta1.dV", "beta1.dF","beta1.dUV","beta1.dUF","mu.sV","mu.sF","mu.sUV","mu.sUF","mu.fV","mu.fF","mu.fUV","mu.fUF","mu.dV","mu.dF","mu.dUV","mu.dUF","sigma.sV2","sigma.sF2","sigma.sUV2","sigma.sUF2","sigma.fV2","sigma.fF2","sigma.fUV2","sigma.fUF2","sigma.dV2","sigma.dF2","sigma.dUV2","sigma.dUF2")

# MCMC settings
ni <- 2500
nt <- 5
nb <- 500
nc <- 3

# Call JAGS from R #
#complex model
ms4.rf <- jags(jags.data, inits, parameters, "ms-ranef4stages.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T)
print(ms4.rf, digits=3)

###save all samples
mod.samples_4stage<- as.data.frame(do.call("rbind", ms3.rf$samples))
write.csv(mod.samples_4stage,"msmod_samples_4stage.csv",row.names=T)

##Table of model summary, with betas
mod.sum<-ms4.rf$summary
write.csv(mod.sum,"isotria4stagemodsum.csv", row.names=T)
mod.table<-round(subset(mod.sum, select=c("mean","sd","Rhat", "n.eff")), digits=3)
mod.table2<-rbind(mod.table[which(substr(rownames(mod.table),1,2)=="mu"),],mod.table[which(substr(rownames(mod.table),1,4)=="beta"),])
write.csv(mod.table2,"isotria4stage_modtable.csv")
