#Multistate model for Isotria medioloides Alton, NH population
#Data provided by Bill Brumback
#Coding by Ailene Ettinger with help frmo Andy Royle and Elizabeth Crone
#This file is to try to debug my 4 state multistate model 
#Try a really simple model (no random effects, etc, no groups) with four states
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
#add 4s for when dormant stage =UF
isoallUF_CH=tapply(isoinds_prev2$dormstageUF,list(isoinds_prev2$UniqueID,isoinds_prev2$Year),sum)#reproductive status
isoall_CH.ms=isoallUF_CH+isoall_CH
isoall_CH.ms[which(isoall_CH.ms==4)]<-3#BAsed on what Nate said, this should not be explicitly coded this way because we can't see this stage
n.occasions<- dim(isoall_CH.ms)[2]
get.first <- function(x) min(which(x!=0))

f<-apply(isoall_CH.ms,1,get.first)#first occasion of marking
#Fit simple model: no groups, no random effects
sink("ms-simple4stages.jags")
cat("
    
    model {
    # -------------------------------------------------
    # Parameters:
    # sV: survival probability for nonreproductive plants
    # sF: survival probability for reproductive plants
    # sUV: survival probability for dormant plants that were vegetative above ground
    # sUF: survival probability for dormant plants that were reproductive above ground
    #  dUF: dormancy probability for dormant plants that were reproductive above ground
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
    sV[i,t] <- mean.sV
    sF[i,t] <- mean.sF
    sUV[i,t] <- mean.sUV
    sUF[i,t] <- mean.sUF
    fV[i,t] <- mean.fV
    fF[i,t] <- mean.fF
    fUV[i,t] <- mean.fUV
    fUF[i,t] <- mean.fUF
    dV[i,t] <- mean.dV
    dF[i,t] <-mean.dF
    dUV[i,t] <- mean.dUV
    dUF[i,t] <- mean.dUF
    }#t
    }#i
    mean.sV~dunif(0,1)# Priors for mean survival for veg plants
    mean.sF~dunif(0,1)# Priors for mean survival for rep plants
    mean.sUV~dunif(0,1)# Priors for mean survival for dorm(prev veg) plants
    mean.sUF~dunif(0,1)# Priors for mean survival dorm (prev rep) plants
    mean.fV~dunif(0,1)# Priors for mean prob of rep for veg plants
    mean.fF~dunif(0,1)# Priors for mean  prob of rep for rep plants
    mean.fUV~dunif(0,1)# Priors for mean prob of rep for dorm (prev veg) plants
    mean.fUF~dunif(0,1)# Priors for mean prob of rep for dorm (prev rep) plants
    mean.dV~dunif(0,1)# Priors for mean dormancy for veg plants
    mean.dF~dunif(0,1)# Priors for mean dormancy for rep plants
    mean.dUV~dunif(0,1)# Priors for mean  dormancy for dorm (prev veg) plants
    mean.dUF~dunif(0,1)# Priors for mean  dormancy for dorm (prev rep) plants

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
  return(state)
}
##I need to recode this so that the starting value is always a 2 for unobserved instances when the last seen state was 2 and a 1 for unobserved instances when the last seen state was 1
ms.init.z <- function(ch, f){#ms.init.z gives starting values of 1 or 2 to all unknown states
  for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
  n1 <- min(which(ch[i,]<3))#
  states <- max(ch, na.rm = TRUE)
  #the following is the beginning of my work to get the starting values correct:
  #anything that was last observed as a 1 should not get a 2 for the starting value when it is unobserved
  #anything last observed as a 2 should not get a 1 for the starting value when it is unobserved
  for (i in 1:dim(st[1]){
   rle(is.na(ch[i,]))$values
    s<-rle(ch[i,])$values[which(!is.na(rle(ch[i,])$values))]
    l<-rle(ch[i,])$lengths[which(!is.na(rle(ch[i,])$value))]
    
    for(j in 1:length(s)){
      ss<-s[j]
      ls<-rep(ss,times=l[j+1])
    }
    
  }
        # known.states <- 1:2
      v <- which(ch>=3)#which occurences are unknown states (=3)
      st<-ch
      st[-v] <- NA#everything else besides the 3 gets an NA
      st[v] <- sample(known.states, length(v), replace = TRUE)##gives all unseen occurrences (3s and 4s) a starting value that is a random sampling of either 1 or 2
      st[i,n1] <- NA#make first observance an NA (not a 3)
  return(st)
}

#ch[1]

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
jags.data <- list(y = isoall_CH.ms, f = f, n.occasions = dim(isoall_CH.ms)[2], nind = dim(isoall_CH.ms)[1], z = known.state.ms(isoall_CH.ms))

# Initial values
inits<-function(){list(mean.sV=c(.5), mean.sF=c(.5),mean.sUV=c(.5),mean.sUF=c(.5),mean.fV=c(.5),mean.fF=c(.5),mean.fUV=c(.5),mean.fUF=c(.5),mean.dV=c(.5),mean.dF=c(.5),mean.dUV=c(.5),mean.dUF=c(.5),z = zst1)}
# Parameters monitored
parameters <- c("mean.sV","mean.sF","mean.sUV","mean.sUF", "mean.fV","mean.fF","mean.fUV","mean.fUF","mean.dV","mean.dF","mean.dUV","mean.dUF")

# MCMC settings
ni <- 2500
nt <- 5
nb <- 500
nc <- 3

# Call JAGS from R #
#complex model
ms4.simple <- jags(jags.data, inits, parameters, "ms-simple4stages.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T)
print(ms4.rf, digits=3)
