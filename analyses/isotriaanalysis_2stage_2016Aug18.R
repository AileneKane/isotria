#Population model for Isotria medioloides Alton NH population
#Data provided by Bill Brumback
#Coding by Ailene Ettinger with help frmo Andy Roly and Elizabeth Crone
#this file has code for a multistate models to analyze two groups (x and y) from a population of Isotria medeoloides
#the model here allows vital rate estimates in the multistate model to vary independently among years in their random effect structure
#
#setwd("~/isotria") #at usgs
#setwd("~/Dropbox/Documents/Work/projects_inprog/isotria/analyses/MultistateModels")
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
isoindsX$Repro<-NA
isoindsX[which(isoindsX$NoFrStems>0|isoindsX$NoFlStems>0|isoindsX$NoArrStems>0),]$Repro=1
isoindsX[which(isoindsX$NoFrStems==0&isoindsX$NoFlStems==0&isoindsX$NoArrStems==0&isoindsX$Stage!="D"),]$Repro=0
isoindsY$Repro<-NA
isoindsY[which(isoindsY$NoFrStems>0|isoindsY$NoFlStems>0|isoindsY$NoArrStems>0),]$Repro=1
isoindsY[which(isoindsY$NoFrStems==0 & isoindsY$NoFlStems==0 & isoindsY$NoArrStems==0&isoindsY$Stage!="D"),]$Repro=0

isoalleme_CH=tapply(isoindsXY$Emerg,list(isoindsXY$UniqueID,isoindsXY$Year),sum)#emergence
isoallrep_CH=tapply(isoindsXY$Repro,list(isoindsXY$UniqueID,isoindsXY$Year),sum)#reproductive status
isoall_CH.ms=isoalleme_CH+isoallrep_CH
isoall_CH.ms[which(isoall_CH.ms==0)]=3#replace 0s wih 3
isoall_CH.ms[73,23]<-NA#remove case where first observation is a zero/3- makes no sense
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

sink("ms-ranef3a.jags")
cat("

model {
    # -------------------------------------------------
    # Parameters:
    # phiA: survival probability for nonreproductive plants
    # phiB: survival probability for reproductive plants
    # psiAB: transition probability from nonreproductive to fruiting
    # psiBA: transition probability from reproductive to nonfruiting
    # pA: emergence probability for nonreproductive plants
    # pB: emergence probability for reproductive plants
    # -------------------------------------------------
    # States (S):
    # 1 alive and nonreproductive
    # 2 alive and reproductive
    # 3 dead
    # Observations (O):  
    # 1 seen nonreproductive
    # 2 seen reproductive
    # 3 not seen
    # -------------------------------------------------
    
    # Model, priors and constraints
    for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
    logit(phiA[i,t]) <- eta.phiA[group[i],t]
    logit(phiB[i,t]) <- eta.phiB[group[i],t]
    logit(psiAB[i,t]) <- eta.psiAB[group[i],t]
    logit(psiBA[i,t]) <- eta.psiBA[group[i],t]
    logit(pA[i,t]) <- eta.pA[group[i],t]
    logit(pB[i,t]) <- eta.pB[group[i],t]    
    }#t
    }#i
    for (g in 1:2){#group  
    for (t in 1:(n.occasions-1)){
    eta.phiA[g,t]  <-mu.phiA[g]+beta.phiA[g]*x[g,t]+beta1.phiA[g]*x1[g,t]+epsilon.phiA[g,t]
    eta.phiB[g,t]  <-mu.phiB[g]+beta.phiB[g]*x[g,t]+beta1.phiB[g]*x1[g,t]+epsilon.phiB[g,t]
    eta.psiAB[g,t]  <-mu.psiAB[g]+beta.psiAB[g]*x[g,t]+beta1.psiAB[g]*x1[g,t]+epsilon.psiAB[g,t]
    eta.psiBA[g,t]  <-mu.psiBA[g]+beta.psiBA[g]*x[g,t]+beta1.psiBA[g]*x1[g,t]+epsilon.psiBA[g,t]
    eta.pA[g,t]  <-mu.pA[g]+beta.pA[g]*x[g,t]+beta1.pA[g]*x1[g,t]+epsilon.pA[g,t]
    eta.pB[g,t]  <-mu.pB[g]+beta.pB[g]*x[g,t]+beta1.pB[g]*x1[g,t]+epsilon.pB[g,t]
    epsilon.phiA[g,t]~dnorm(0,tau.phiA[g])#could move all the mean stuff to this part, might help model mix better? hierarchical centering
    epsilon.phiB[g,t]~dnorm(0,tau.phiB[g])
    epsilon.psiAB[g,t]~dnorm(0,tau.psiAB[g])
    epsilon.psiBA[g,t]~dnorm(0,tau.psiBA[g])
    epsilon.pA[g,t]~dnorm(0,tau.pA[g])
    epsilon.pB[g,t]~dnorm(0,tau.pB[g])
    }#t
    mean.phiA[g]~dunif(0,1)# Priors for mean group-specific survival for veg plants
    mu.phiA[g]<-log(mean.phiA[g]/(1-mean.phiA[g]))  
    mean.phiB[g]~dunif(0,1)# Priors for mean group-specific survivalsurvival for rep plants
    mu.phiB[g]<-log(mean.phiB[g]/(1-mean.phiB[g]))    
    mean.psiAB[g]~dunif(0,1)# Priors for mean group-specific transition from veg to rep plants
    mu.psiAB[g]<-log(mean.psiAB[g]/(1-mean.psiAB[g]))       
    mean.psiBA[g]~dunif(0,1)# Priors for mean group-specific transition from rep to veg plants
    mu.psiBA[g]<-log(mean.psiBA[g]/(1-mean.psiBA[g]))    
    mean.pA[g]~dunif(0,1)#Priors for mean group-specific detection (emergence) probability or veg plants
    mu.pA[g]<-log(mean.pA[g]/(1-mean.pA[g]))   
    mean.pB[g]~dunif(0,1)#Priors for mean group-specific detection (emergence) probability or veg plants
    mu.pB[g]<-log(mean.pB[g]/(1-mean.pB[g]))    
    sigma.phiA[g]~dunif(0,10)#temporal variance for veg plants survival
    tau.phiA[g]<-pow(sigma.phiA[g],-2)
    sigma.phiA2[g]<-pow(sigma.phiA[g],2)
    sigma.phiB[g]~dunif(0,10)#temporal variance for rep plants survival
    tau.phiB[g]<-pow(sigma.phiB[g],-2)
    sigma.phiB2[g]<-pow(sigma.phiB[g],2)
    sigma.psiAB[g]~dunif(0,10)#temporal variance for veg plants transition to rep
    tau.psiAB[g]<-pow(sigma.psiAB[g],-2)
    sigma.psiAB2[g]<-pow(sigma.psiAB[g],2)
    sigma.psiBA[g]~dunif(0,10)#temporal variance for rep plants transition to veg
    tau.psiBA[g]<-pow(sigma.psiBA[g],-2)
    sigma.psiBA2[g]<-pow(sigma.psiBA[g],2)
    sigma.pA[g]~dunif(0,10)#temporal variance for veg plants recapture (emergence)
    tau.pA[g]<-pow(sigma.pA[g],-2)
    sigma.pA2[g]<-pow(sigma.pA[g],2)
    sigma.pB[g]~dunif(0,10)#temporal variance for rep plants recapture (emergence)
    tau.pB[g]<-pow(sigma.pB[g],-2)
    sigma.pB2[g]<-pow(sigma.pB[g],2)
    beta.phiA[g]~dnorm(0,0.001)I(-10,10)
    beta.phiB[g]~dnorm(0,0.001)I(-10,10)
    beta.psiAB[g]~dnorm(0,0.001)I(-10,10)
    beta.psiBA[g]~dnorm(0,0.001)I(-10,10)
    beta.pA[g]~dnorm(0,0.001)I(-10,10)
    beta.pB[g]~dnorm(0,0.001)I(-10,10)
    beta1.phiA[g]~dnorm(0,0.001)I(-10,10)
    beta1.phiB[g]~dnorm(0,0.001)I(-10,10)
    beta1.psiAB[g]~dnorm(0,0.001)I(-10,10)
    beta1.psiBA[g]~dnorm(0,0.001)I(-10,10)
    beta1.pA[g]~dnorm(0,0.001)I(-10,10)
    beta1.pB[g]~dnorm(0,0.001)I(-10,10)
    #calculate probabilities of four differnt groups to look at later... 
    logit(phiA0[g])<- mu.phiA[g]
    logit(phiA1[g])<- mu.phiA[g] + beta.phiA[g]
    logit(phiB0[g])<- mu.phiB[g]
    logit(phiB1[g])<- mu.phiB[g] + beta.phiB[g]
    logit(psiA0[g])<- mu.psiAB[g]
    logit(psiA1[g])<- mu.psiAB[g] + beta.psiAB[g]
    logit(psiB0[g])<- mu.psiBA[g]
    logit(psiB1[g])<- mu.psiBA[g] + beta.psiBA[g]
    logit(pA0[g])<- mu.pA[g]
    logit(pA1[g])<- mu.pA[g] + beta.pA[g]
    logit(pB0[g])<- mu.pB[g]
    logit(pB1[g])<- mu.pB[g] + beta.pB[g]
    }#g
    
    # Define state-transition and observation matrices
    for (i in 1:nind){  
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    ps[1,i,t,1] <- phiA[i,t] * (1-psiAB[i,t])
    ps[1,i,t,2] <- phiA[i,t] * psiAB[i,t]
    ps[1,i,t,3] <- 1-phiA[i,t]
    ps[2,i,t,1] <- phiB[i,t] * psiBA[i,t]
    ps[2,i,t,2] <- phiB[i,t] * (1-psiBA[i,t])
    ps[2,i,t,3] <- 1-phiB[i,t]
    ps[3,i,t,1] <- 0
    ps[3,i,t,2] <- 0
    ps[3,i,t,3] <- 1
    
    # Define probabilities of O(t) given S(t) ##first coumn is observed state, last is true state
    po[1,i,t,1] <- pA[i,t]
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 1-pA[i,t]
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- pB[i,t]
    po[2,i,t,3] <- 1-pB[i,t]
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- 1
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
ms.init.z <- function(ch, f){#ms.init.z gives starting values of 1 or 2 to all unknown states
  for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
  n1 <- min(which(ch[i,]<3))#
  states <- max(ch, na.rm = TRUE)
  known.states <- 1:(states-1)
  v <- which(ch==states)#which occurences are unknown states (=3)
  ch[-v] <- NA#everyhing else besides the 3s gets an NA
  ch[v] <- sample(known.states, length(v), replace = TRUE)##gives all unseen occurrences (3s) a starting value that is a random sampling of either 1 or 2
  #ch[i,n1] <- NA#make first observance an NA (ot a 3)
  return(ch)
}

ch=isoall_CH.ms
y = isoall_CH.ms
y[73,24]=NA
zst1=ms.init.z(ch,f)
#fixes necessary to get model to run with zst1:
zst1[4,24:31]<-NA
zst1[6,24:31]<-NA
zst1[5,24:31]<-NA
zst1[60,31]<-NA
#zst1[73,25:26]<-NA
zst1[119,21:31]<-NA

# Bundle data
jags.data <- list(y = isoall_CH.ms, f = f, n.occasions = dim(isoall_CH.ms)[2], nind = dim(isoall_CH.ms)[1], z = known.state.ms(isoall_CH.ms), group=group, x=time_log, x1=logged_yrs2)

# Initial values
inits3a<-function(){list(mean.phiA=c(.5,.5), mean.phiB=c(.5,.5),mean.psiAB=c(.5,.5),mean.psiBA=c(.5,.5),mean.pA=c(.5,.5),mean.pB=c(.5,.5),sigma.phiB=c(1,1),sigma.psiAB=c(1,1),sigma.psiBA=c(1,1),sigma.pA=c(1,1),sigma.pB=c(1,1),z = zst1)}#Error in jags.model(file = model.file, data = data, inits = inits, n.chains = n.chains,  : #Error in node y[1,3]

# Parameters monitored
parameters <- c("mean.phiA","mean.phiB", "mean.psiAB","mean.psiBA","mean.pA","mean.pB","beta.phiA","beta.phiB", "beta.pA", "beta.pB", "beta.psiAB","beta.psiBA","phiA0","phiA1","phiB0","phiB1","psiA0","psiA1","psiB0","psiB1","pA0","pA1","pB0","pB1","beta1.phiA","beta1.phiB", "beta1.pA", "beta1.pB", "beta1.psiAB","beta1.psiBA","mu.phiA","mu.phiB","mu.pA","mu.pB","mu.psiAB","mu.psiBA","sigma.phiA2","sigma.phiB2","sigma.pA2","sigma.pB2","sigma.phiB2","sigma.psiAB2","sigma.psiBA2")

# MCMC settings
ni <- 25000
nt <- 10
nb <- 5000
nc <- 3

# Call JAGS from R #
#complex model
ms.rf3a <- jags(jags.data, inits3a, parameters, "ms-ranef3a.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T)
print(ms.rf3a, digits=3)

###save all samples
mod.samples_comp <- as.data.frame(do.call("rbind", ms.rf3a$samples))
write.csv(mod.samples_comp,"msmod_samples_complex.csv",row.names=T)
mod.samples_comp<-read.csv("msmod_samples_complex.csv", header=T)

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
windows(height=6,width=10)
#quartz(height=6,width=10)
par(mfrow=c(2,2))
hist(propdorm_Xpre, xlim=c(0,1))
hist(propdorm_Ypre,xlim=c(0,1))
hist(propdorm_Xpost,xlim=c(0,1))
hist(propdorm_Ypost,xlim=c(0,1))

mean(propdorm_Xpre);sd(propdorm_Xpre)#0.257 (0.052 plants dormant in uncleared prior to clearing
mean(propdorm_Ypre);sd(propdorm_Ypre)#0.219 (0.050)plants dormant in cleared prior to clearing
mean(propdorm_Xpost);sd(propdorm_Xpost)#0.10 (0.06) plants dormant in uncleared post clearing
mean(propdorm_Ypost);sd(propdorm_Ypost)#0.094 (0.069) plants dormant in cleared post clearing

####Now lifespan:
##test:
#phiV=phiV_Ypost
#veg.rep=veg.rep_Ypost
#pdormV=pdormV_Ypost
#phiR=phiR_Ypost
#rep.veg=rep.veg_Ypost
#pdormR=pdormR_Ypost
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
    n0 = c(0,0,1000,0)
    nsum = array()
    flwrsum = array()
    for(j in 1:200){
    n1 = tmx%*%n0
    nsum[j] = sum(n1)
    flwrsum[j] = n1[4]
    n0 = n1
    }#
  lifespan_med[i]= min(which(nsum <500)) # 
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

windows(height=6,width=10)
#quartz(height=6,width=10)
par(mfrow=c(2,2))
barplot(colMeans(lengthdorm_flow_Xpre), names.arg = 1:11, xlab = "length of dormancy (years)", ylab = "% of bouts, rep. plants", main="Control, pre")
barplot(colMeans(lengthdorm_flow_Ypre), names.arg = 1:11, xlab = "length of dormancy (years)", ylab = "% of bouts, rep. plants", main="Cleared, pre")
barplot(colMeans(lengthdorm_flow_Xpost), names.arg = 1:11, xlab = "length of dormancy (years)", ylab = "% of bouts, rep. plants", main="Control, post")
barplot(colMeans(lengthdorm_flow_Ypost), names.arg = 1:11, xlab = "length of dormancy (years)", ylab = "% of bouts, rep. plants", main="Cleared, post")

###Figures
#2x2table for each vital rate with first column control, second column logged
#first row before logging, second row after logging
surv_veg<-as.data.frame(rbind(ms.rf3a$mean$phiA0,ms.rf3a$mean$phiA1))
surv_rep<-as.data.frame(rbind(ms.rf3a$mean$phiB0,ms.rf3a$mean$phiB1))
trans_vr<-as.data.frame(rbind(ms.rf3a$mean$psiA0,ms.rf3a$mean$psiA1))
trans_rv<-as.data.frame(rbind(ms.rf3a$mean$psiB0,ms.rf3a$mean$psiB1))
emer_veg<-as.data.frame(rbind(ms.rf3a$mean$pA0,ms.rf3a$mean$pA1))
emer_rep<-as.data.frame(rbind(ms.rf3a$mean$pB0,ms.rf3a$mean$pB1))
#use below code if model not loaded
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

surv_veg_q2.5<-c(ms.rf3a$q2.5$phiA0,ms.rf3a$q2.5$phiA1)
surv_rep_q2.5<-c(ms.rf3a$q2.5$phiB0,ms.rf3a$q2.5$phiB1)
trans_vr_q2.5<-c(ms.rf3a$q2.5$psiA0,ms.rf3a$q2.5$psiA1)
trans_rv_q2.5<-c(ms.rf3a$q2.5$psiB0,ms.rf3a$q2.5$psiB1)
emer_veg_q2.5<-c(ms.rf3a$q2.5$pA0,ms.rf3a$q2.5$pA1)
emer_rep_q2.5<-c(ms.rf3a$q2.5$pB0,ms.rf3a$q2.5$pB1)
surv_veg_q97.5<-c(ms.rf3a$q97.5$phiA0,ms.rf3a$q97.5$phiA1)
surv_rep_q97.5<-c(ms.rf3a$q97.5$phiB0,ms.rf3a$q97.5$phiB1)
trans_vr_q97.5<-c(ms.rf3a$q97.5$psiA0,ms.rf3a$q97.5$psiA1)
trans_rv_q97.5<-c(ms.rf3a$q97.5$psiB0,ms.rf3a$q97.5$psiB1)
emer_veg_q97.5<-c(ms.rf3a$q97.5$pA0,ms.rf3a$q97.5$pA1)
emer_rep_q97.5<-c(ms.rf3a$q97.5$pB0,ms.rf3a$q97.5$pB1)

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
x<-c(1,2,3,4)
x2<-c(1,3,2,4)
windows(height=6,width=7)
quartz(height=6,width=7)
par(mfrow=c(3,1),mar=c(.5,4,1,.5), oma=c(3,.5,.5,.5))
#survival
plot(x-.05,c(surv_veg$control,surv_veg$logged), pch=21, bg="black", ylim=c(0,1), ylab="Survival", xaxt="n", cex=1.5, xlab="", xlim=c(0.75,4.25), cex.lab=1.2)
abline(v=2.5,lty=1, lwd=2)
abline(v=1.5,lty=2,col="gray", lwd=2)
abline(v=3.5,lty=2,col="gray", lwd=2)

arrows(x2-.05,surv_veg_q2.5,x2-.05,surv_veg_q97.5, code=0,angle=90, length=0.1)
arrows(x2+.05,surv_rep_q2.5,x2+.05,surv_rep_q97.5, code=0,angle=90, length=0.1)
points(x-.05,c(surv_veg$control,surv_veg$logged),pch=21, bg="black", cex=1.5)
points(x+.05,c(surv_rep$control,surv_rep$logged),pch=21, bg="white", cex=1.5)

axis(side=1,at=c(1.5,3.5),labels=c("Control","Logged" ),line=-15, tick=F, cex.axis=1.2)
legend("bottomright",legend=c("Vegetative", "Reproductive"),pch=21,pt.cex=1.5,pt.bg=c("black","white"), bty="n")
#Dormancy(=1-)Emergence
plot(x-.05,c(1-emer_veg$control,1-emer_veg$logged), pch=21, bg="black", ylim=c(0,1), ylab="Dormancy", xaxt="n", cex=1.5, xlab="", xlim=c(0.75,4.25), cex.lab=1.2)
abline(v=2.5,lty=1, lwd=2)
abline(v=1.5,lty=2,col="gray", lwd=2)
abline(v=3.5,lty=2,col="gray", lwd=2)
arrows(x2-.05,1-emer_veg_q2.5,x2-.05,1-emer_veg_q97.5, code=0,angle=90, length=0.1)
arrows(x2+.05,1-emer_rep_q2.5,x2+.05,1-emer_rep_q97.5, code=0,angle=90, length=0.1)
points(x+.05,c(1-emer_rep$control,1-emer_rep$logged),pch=21, bg="white", cex=1.5)
points(x-.05,c(1-emer_veg$control,1-emer_veg$logged),pch=21, bg="black", cex=1.5)
#transition
plot(x-.05,c(trans_vr$control,trans_vr$logged), pch=21, bg="black", ylim=c(0,1), ylab="Transition", xaxt="n", cex=1.5, xlab="", xlim=c(0.75,4.25), cex.lab=1.2)
abline(v=2.5,lty=1, lwd=2)
arrows(x2-.05,trans_vr_q2.5,x2-.05,trans_vr_q97.5, code=0,angle=90, length=0.1)
arrows(x2+.05,trans_rv_q2.5,x2+.05,trans_rv_q97.5, code=0,angle=90, length=0.1)
points(x-.05,c(trans_vr$control,trans_vr$logged),pch=21, bg="black", cex=1.5)
points(x+.05,c(trans_rv$control,trans_rv$logged),pch=21, bg="white", cex=1.5)
abline(v=1.5,lty=2,col="gray", lwd=2)
abline(v=3.5,lty=2,col="gray", lwd=2)
legend("topright",legend=c("Vegetative->Reproductive", "Reproductive->Vegetative"),pch=21,pt.cex=1.5,pt.bg=c("black","white"), bty="n")
axis(side=1,at=x,labels=c("pre","post","pre","post"))
axis(side=1,at=x,labels=c("(1982-1997)","(1998-2015)","(1982-1997)","(1998-2015)"), line=1.2,tick=F)

#####Figure 4,  length of dormancy, lifepsan, etc
x<-c(1,2,3,4)
windows(height=6,width=7)
quartz(height=6,width=7)
par(mfrow=c(3,1),mar=c(.5,4,1,.5), oma=c(3,.5,.5,.5))

#lifespan
plot(x,c(mean(lifespan_Xpre, na.rm=T),mean(lifespan_Xpost2, na.rm=T),mean(lifespan_Ypre, na.rm=T),mean(lifespan_Ypost2, na.rm=T)), pch=21, bg="gray", ylim=c(0,100), ylab="Lifespan (yrs)", xaxt="n", cex=1.5, xlab="", xlim=c(0.75,4.25), cex.lab=1.2)
abline(v=2.5,lty=1, lwd=2)
abline(v=1.5,lty=2,col="gray", lwd=2)
abline(v=3.5,lty=2,col="gray", lwd=2)
arrows(x,c(mean(lifespan_Xpre, na.rm=T)-sd(lifespan_Xpre, na.rm=T),mean(lifespan_Xpost2, na.rm=T)-sd(lifespan_Xpost2, na.rm=T),mean(lifespan_Ypre, na.rm=T)-sd(lifespan_Ypre, na.rm=T),mean(lifespan_Ypost2, na.rm=T)-sd(lifespan_Ypost2, na.rm=T)),x,c(mean(lifespan_Xpre, na.rm=T)+sd(lifespan_Xpre, na.rm=T),mean(lifespan_Xpost2, na.rm=T)+sd(lifespan_Xpost2, na.rm=T),mean(lifespan_Ypre, na.rm=T)+sd(lifespan_Ypre, na.rm=T),mean(lifespan_Ypost2, na.rm=T)+sd(lifespan_Ypost2, na.rm=T)), code=0,angle=90, length=0.1)
points(x,c(mean(lifespan_Xpre, na.rm=T),mean(lifespan_Xpost2, na.rm=T),mean(lifespan_Ypre, na.rm=T),mean(lifespan_Ypost2, na.rm=T)), pch=21, bg="gray", cex=1.5)
axis(side=1,at=c(1.5,3.5),labels=c("Control","Logged" ),line=-14.5, tick=F, cex.axis=1.2)

#proportion of plants dormant
plot(x,c(mean(propdorm_Xpre, na.rm=T),mean(propdorm_Xpost, na.rm=T),mean(propdorm_Ypre, na.rm=T),mean(propdorm_Ypost, na.rm=T)), pch=21, bg="gray", ylim=c(0,1), ylab="Proportion Dormant", xaxt="n", cex=1.5, xlab="", xlim=c(0.75,4.25), cex.lab=1.2)
abline(v=2.5,lty=1, lwd=2)
abline(v=1.5,lty=2,col="gray", lwd=2)
abline(v=3.5,lty=2,col="gray", lwd=2)
arrows(x,c(mean(propdorm_Xpre, na.rm=T)-sd(propdorm_Xpre, na.rm=T),mean(propdorm_Xpost, na.rm=T)-sd(propdorm_Xpost, na.rm=T),mean(propdorm_Ypre, na.rm=T)-sd(propdorm_Ypre, na.rm=T),mean(propdorm_Ypost, na.rm=T)-sd(propdorm_Ypost, na.rm=T)),x,c(mean(propdorm_Xpre, na.rm=T)+sd(propdorm_Xpre, na.rm=T),mean(propdorm_Xpost, na.rm=T)+sd(propdorm_Xpost, na.rm=T),mean(propdorm_Ypre, na.rm=T)+sd(propdorm_Ypre, na.rm=T),mean(propdorm_Ypost, na.rm=T)+sd(propdorm_Ypost, na.rm=T)), code=0,angle=90, length=0.1)
points(x,c(mean(propdorm_Xpre, na.rm=T),mean(propdorm_Xpost, na.rm=T),mean(propdorm_Ypre, na.rm=T),mean(propdorm_Ypost, na.rm=T)), pch=21, bg="gray", cex=1.5)
#Length of dormancy, starting from veg (black) or rep (white)
plot(x-.05,c(mean(lengthdorm_Xpre, na.rm=T),mean(lengthdorm_Xpost, na.rm=T),mean(lengthdorm_Ypre, na.rm=T),mean(lengthdorm_Ypost, na.rm=T)), pch=21, bg="black", ylim=c(0,2), ylab="Dormancy length (yrs)", xaxt="n", cex=1.5, xlab="", xlim=c(0.75,4.25), cex.lab=1.2)
abline(v=2.5,lty=1, lwd=2)
abline(v=1.5,lty=2,col="gray", lwd=2)
abline(v=3.5,lty=2,col="gray", lwd=2)
arrows(x-.05,c(mean(lengthdorm_Xpre, na.rm=T)-sd(lengthdorm_Xpre, na.rm=T),mean(lengthdorm_Xpost, na.rm=T)-sd(lengthdorm_Xpost, na.rm=T),mean(lengthdorm_Ypre, na.rm=T)-sd(lengthdorm_Ypre, na.rm=T),mean(lengthdorm_Ypost, na.rm=T)-sd(lengthdorm_Ypost, na.rm=T)),x-.05,c(mean(lengthdorm_Xpre, na.rm=T)+sd(lengthdorm_Xpre, na.rm=T),mean(lengthdorm_Xpost, na.rm=T)+sd(lengthdorm_Xpost, na.rm=T),mean(lengthdorm_Ypre, na.rm=T)+sd(lengthdorm_Ypre, na.rm=T),mean(lengthdorm_Ypost, na.rm=T)+sd(lengthdorm_Ypost, na.rm=T)), code=0,angle=90, length=0.1)
points(x-.05,c(mean(lengthdorm_Xpre, na.rm=T),mean(lengthdorm_Xpost, na.rm=T),mean(lengthdorm_Ypre, na.rm=T),mean(lengthdorm_Ypost, na.rm=T)), pch=21, bg="black", cex=1.5)
arrows(x+.05,c(mean(lengthdorm_flow_Xpre, na.rm=T)-sd(lengthdorm_flow_Xpre, na.rm=T),mean(lengthdorm_flow_Xpost, na.rm=T)-sd(lengthdorm_flow_Xpost, na.rm=T),mean(lengthdorm_flow_Ypre, na.rm=T)-sd(lengthdorm_flow_Ypre, na.rm=T),mean(lengthdorm_flow_Ypost, na.rm=T)-sd(lengthdorm_flow_Ypost, na.rm=T)),x+.05,c(mean(lengthdorm_flow_Xpre, na.rm=T)+sd(lengthdorm_flow_Xpre, na.rm=T),mean(lengthdorm_flow_Xpost, na.rm=T)+sd(lengthdorm_flow_Xpost, na.rm=T),mean(lengthdorm_flow_Ypre, na.rm=T)+sd(lengthdorm_flow_Ypre, na.rm=T),mean(lengthdorm_flow_Ypost, na.rm=T)+sd(lengthdorm_flow_Ypost, na.rm=T)), code=0,angle=90, length=0.1)
points(x+.05,c(mean(lengthdorm_flow_Xpre, na.rm=T),mean(lengthdorm_flow_Xpost, na.rm=T),mean(lengthdorm_flow_Ypre, na.rm=T),mean(lengthdorm_flow_Ypost, na.rm=T)), pch=21, bg="white", cex=1.5)
axis(side=1,at=x,labels=c("pre","post","pre","post"))
axis(side=1,at=x,labels=c("(1982-1997)","(1998-2015)","(1982-1997)","(1998-2015)"), line=1.2,tick=F)
legend("bottomright",legend=c("Vegetative", "Reproductive"),pch=21,pt.cex=1.5,pt.bg=c("black","white"), bty="n")

##Table of model summary, with betas
mod.sum<-ms.rf3a$summary
write.csv(mod.sum,"isotria2stagemodsum_complex.csv", row.names=T)
mod.table<-round(subset(mod.sum, select=c("mean","sd","Rhat", "n.eff")), digits=3)
mod.table2<-rbind(mod.table[which(substr(rownames(mod.table),1,2)=="mu"),],mod.table[which(substr(rownames(mod.table),1,4)=="beta"),])
write.csv(mod.table2,"isotria_modtable.csv")

##Table of model summary, with betas
mod.sum<-ms.rf3a$summary
write.csv(mod.sum,"isotria2stagemodsum_complex.csv", row.names=T)
mod.table<-round(subset(mod.sum, select=c("mean","sd","Rhat", "n.eff")), digits=3)
mod.table2<-rbind(mod.table[which(substr(rownames(mod.table),1,2)=="mu"),],mod.table[which(substr(rownames(mod.table),1,4)=="beta"),])
write.csv(mod.table2,"isotria_modtable.csv")

#Simulate data with known parameters to test the ability of my model to recover them
#start with simple dataset with no variation between groups x and y
#first, define mean survival, transisions, and emergence probability, 
#as well as the number of occassions, states, observations, and marked individuals
phiA0<-0.5
#phiA1<-0.7
phiB0<-0.7
#phiB1<-0.9
pA0<-0.5
#pA1<-0.7
pB0<-0.3
#pB1<-0.5
phiA0<-0.2
#phiA1<-0.4
n.occasions<-30
n.states<-2
n.obs<-15
marked<-matrix(0,ncol=n.states,nrow=n.occasions)
marked[,1]<-rep(200,n.occasions)
#Now define matrices with survival, transition, and emergence probailites
#these are 4-dimensions amtraicies with:
#Dimension 1: state of departure
#Dimension 2: state of arrival
#Dimension 3: individual
#Dimension 4:time
totrel<-sum(marked)*(n.occasions-1)
PSI.STATE<-array(NA,dim=c(n.states,n.states,totrel,n.occasions-1))
for(i in 1:totrel){
  PSI.STATE[,,i,t]<-matric(c(
    s*
  ))
}

