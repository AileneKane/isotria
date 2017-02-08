#Multistate model for Isotria medioloides Alton NH population
#Data provided by Bill Brumback
#Coding by Ailene Ettinger with help frmo Andy Roly and Elizabeth Crone
#this file has code for a multistate models to analyze two groups (x and y) from a population of Isotria medeoloides
#the model here allows vital rate estimates in the multistate model to vary independently among years in their random effect structure
#
#setwd("~/isotria") #at usgs
setwd("~/git/isotria/analyses")
rm(list=ls()) 
options(stringsAsFactors=FALSE)

library(rjags)
library(jagsUI)
library(lattice)
library(coda)
library(boot)
#yupdate.packages()
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

sink("ms-ranef2stages.jags")
cat("

model {
    # -------------------------------------------------
    # Parameters:
    # sV: survival probability for nonreproductive plants
    # sF: survival probability for reproductive plants
    # fV: reproduction probability for vegetative plants
    # fF: reproduction probability for reproductive plants
    
    # sV: survival probability for nonreproductive plants
    # sF: survival probability for reproductive plants
    # fV: transition probability from nonreproductive to fruiting
    # fFA: transition probability from reproductive to nonfruiting
    # pV: emergence probability for nonreproductive plants
    # pF: emergence probability for reproductive plants
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
    logit(sV[i,t]) <- eta.sV[group[i],t]
    logit(sF[i,t]) <- eta.sF[group[i],t]
    logit(fV[i,t]) <- eta.fV[group[i],t]
    logit(fFA[i,t]) <- eta.fFA[group[i],t]
    logit(pV[i,t]) <- eta.pV[group[i],t]
    logit(pF[i,t]) <- eta.pF[group[i],t]    
    }#t
    }#i
    for (g in 1:2){#group  
    for (t in 1:(n.occasions-1)){
    eta.sV[g,t]  <-mu.sV[g]+beta.sV[g]*x[g,t]+beta1.sV[g]*x1[g,t]+epsilon.sV[g,t]
    eta.sF[g,t]  <-mu.sF[g]+beta.sF[g]*x[g,t]+beta1.sF[g]*x1[g,t]+epsilon.sF[g,t]
    eta.fV[g,t]  <-mu.fV[g]+beta.fV[g]*x[g,t]+beta1.fV[g]*x1[g,t]+epsilon.fV[g,t]
    eta.fFA[g,t]  <-mu.fFA[g]+beta.fFA[g]*x[g,t]+beta1.fFA[g]*x1[g,t]+epsilon.fFA[g,t]
    eta.pV[g,t]  <-mu.pV[g]+beta.pV[g]*x[g,t]+beta1.pV[g]*x1[g,t]+epsilon.pV[g,t]
    eta.pF[g,t]  <-mu.pF[g]+beta.pF[g]*x[g,t]+beta1.pF[g]*x1[g,t]+epsilon.pF[g,t]
    epsilon.sV[g,t]~dnorm(0,tau.sV[g])#could move all the mean stuff to this part, might help model mix better? hierarchical centering
    epsilon.sF[g,t]~dnorm(0,tau.sF[g])
    epsilon.fV[g,t]~dnorm(0,tau.fV[g])
    epsilon.fFA[g,t]~dnorm(0,tau.fFA[g])
    epsilon.pV[g,t]~dnorm(0,tau.pV[g])
    epsilon.pF[g,t]~dnorm(0,tau.pF[g])
    }#t
    mean.sV[g]~dunif(0,1)# Priors for mean group-specific survival for veg plants
    mu.sV[g]<-log(mean.sV[g]/(1-mean.sV[g]))  
    mean.sF[g]~dunif(0,1)# Priors for mean group-specific survivalsurvival for rep plants
    mu.sF[g]<-log(mean.sF[g]/(1-mean.sF[g]))    
    mean.fV[g]~dunif(0,1)# Priors for mean group-specific transition from veg to rep plants
    mu.fV[g]<-log(mean.fV[g]/(1-mean.fV[g]))       
    mean.fFA[g]~dunif(0,1)# Priors for mean group-specific transition from rep to veg plants
    mu.fFA[g]<-log(mean.fFA[g]/(1-mean.fFA[g]))    
    mean.pV[g]~dunif(0,1)#Priors for mean group-specific detection (emergence) probability or veg plants
    mu.pV[g]<-log(mean.pV[g]/(1-mean.pV[g]))   
    mean.pF[g]~dunif(0,1)#Priors for mean group-specific detection (emergence) probability or veg plants
    mu.pF[g]<-log(mean.pF[g]/(1-mean.pF[g]))    
    sigma.sV[g]~dunif(0,10)#temporal variance for veg plants survival
    tau.sV[g]<-pow(sigma.sV[g],-2)
    sigma.sV2[g]<-pow(sigma.sV[g],2)
    sigma.sF[g]~dunif(0,10)#temporal variance for rep plants survival
    tau.sF[g]<-pow(sigma.sF[g],-2)
    sigma.sF2[g]<-pow(sigma.sF[g],2)
    sigma.fV[g]~dunif(0,10)#temporal variance for veg plants transition to rep
    tau.fV[g]<-pow(sigma.fV[g],-2)
    sigma.fV2[g]<-pow(sigma.fV[g],2)
    sigma.fFA[g]~dunif(0,10)#temporal variance for rep plants transition to veg
    tau.fFA[g]<-pow(sigma.fFA[g],-2)
    sigma.fFA2[g]<-pow(sigma.fFA[g],2)
    sigma.pV[g]~dunif(0,10)#temporal variance for veg plants recapture (emergence)
    tau.pV[g]<-pow(sigma.pV[g],-2)
    sigma.pV2[g]<-pow(sigma.pV[g],2)
    sigma.pF[g]~dunif(0,10)#temporal variance for rep plants recapture (emergence)
    tau.pF[g]<-pow(sigma.pF[g],-2)
    sigma.pF2[g]<-pow(sigma.pF[g],2)
    beta.sV[g]~dnorm(0,0.001)I(-10,10)
    beta.sF[g]~dnorm(0,0.001)I(-10,10)
    beta.fV[g]~dnorm(0,0.001)I(-10,10)
    beta.fFA[g]~dnorm(0,0.001)I(-10,10)
    beta.pV[g]~dnorm(0,0.001)I(-10,10)
    beta.pF[g]~dnorm(0,0.001)I(-10,10)
    beta1.sV[g]~dnorm(0,0.001)I(-10,10)
    beta1.sF[g]~dnorm(0,0.001)I(-10,10)
    beta1.fV[g]~dnorm(0,0.001)I(-10,10)
    beta1.fFA[g]~dnorm(0,0.001)I(-10,10)
    beta1.pV[g]~dnorm(0,0.001)I(-10,10)
    beta1.pF[g]~dnorm(0,0.001)I(-10,10)
    #calculate probabilities of four differnt groups to look at later... 
    logit(sV0[g])<- mu.sV[g]
    logit(sV1[g])<- mu.sV[g] + beta.sV[g]
    logit(sF0[g])<- mu.sF[g]
    logit(sF1[g])<- mu.sF[g] + beta.sF[g]
    logit(fV0[g])<- mu.fV[g]
    logit(fV1[g])<- mu.fV[g] + beta.fV[g]
    logit(fF0[g])<- mu.fFA[g]
    logit(fF1[g])<- mu.fFA[g] + beta.fFA[g]
    logit(pV0[g])<- mu.pV[g]
    logit(pV1[g])<- mu.pV[g] + beta.pV[g]
    logit(pF0[g])<- mu.pF[g]
    logit(pF1[g])<- mu.pF[g] + beta.pF[g]
    }#g
    
    # Define state-transition and observation matrices
    for (i in 1:nind){  
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    ps[1,i,t,1] <- sV[i,t] * (1-fV[i,t])
    ps[1,i,t,2] <- sV[i,t] * fV[i,t]
    ps[1,i,t,3] <- 1-sV[i,t]
    ps[2,i,t,1] <- sF[i,t] * (1-fF[i,t])
    ps[2,i,t,2] <- sF[i,t] * (fF[i,t])
    ps[2,i,t,3] <- 1-sF[i,t]
    ps[3,i,t,1] <- 0
    ps[3,i,t,2] <- 0
    ps[3,i,t,3] <- 1
    
    # Define probabilities of O(t) given S(t) ##first coumn is observed state, last is true state
    po[1,i,t,1] <- pV[i,t]
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 1-pV[i,t]
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- pF[i,t]
    po[2,i,t,3] <- 1-pF[i,t]
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
inits3a<-function(){list(mean.sV=c(.5,.5), mean.sF=c(.5,.5),mean.fV=c(.5,.5),mean.fF=c(.5,.5),mean.pV=c(.5,.5),mean.pF=c(.5,.5),sigma.sF=c(1,1),sigma.fV=c(1,1),sigma.fF=c(1,1),sigma.pV=c(1,1),sigma.pF=c(1,1),z = zst1)}#Error in jags.model(file = model.file, data = data, inits = inits, n.chains = n.chains,  : #Error in node y[1,3]

# Parameters monitored
parameters <- c("mean.sV","mean.sF", "mean.fV","mean.fF","mean.pV","mean.pF","beta.sV","beta.sF", "beta.pV", "beta.pF", "beta.fV","beta.fF","sV0","sV1","sF0","sF1","fV0","fV1","fF0","fF1","pV0","pV1","pF0","pF1","beta1.sV","beta1.sF", "beta1.pV", "beta1.pF", "beta1.fV","beta1.fFA","mu.sV","mu.sF","mu.pV","mu.pF","mu.fV","mu.fFA","sigma.sV2","sigma.sF2","sigma.pV2","sigma.pF2","sigma.sF2","sigma.fV2","sigma.fFA2")

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

##Figures
#2x2table for each vital rate with first column control, second column logged
#first row before logging, second row after logging
surv_veg<-as.data.frame(rbind(ms.rf3a$mean$sV0,ms.rf3a$mean$sV1))
surv_rep<-as.data.frame(rbind(ms.rf3a$mean$sF0,ms.rf3a$mean$sF1))
trans_vr<-as.data.frame(rbind(ms.rf3a$mean$fV0,ms.rf3a$mean$fV1))
trans_rv<-as.data.frame(rbind(ms.rf3a$mean$fF0,ms.rf3a$mean$fF1))
emer_veg<-as.data.frame(rbind(ms.rf3a$mean$pV0,ms.rf3a$mean$pV1))
emer_rep<-as.data.frame(rbind(ms.rf3a$mean$pF0,ms.rf3a$mean$pF1))

surv_veg_q2.5<-c(ms.rf3a$q2.5$sV0,ms.rf3a$q2.5$sV1)
surv_rep_q2.5<-c(ms.rf3a$q2.5$sF0,ms.rf3a$q2.5$sF1)
trans_vr_q2.5<-c(ms.rf3a$q2.5$fV0,ms.rf3a$q2.5$fV1)
trans_rv_q2.5<-c(ms.rf3a$q2.5$fF0,ms.rf3a$q2.5$fF1)
emer_veg_q2.5<-c(ms.rf3a$q2.5$pV0,ms.rf3a$q2.5$pV1)
emer_rep_q2.5<-c(ms.rf3a$q2.5$pF0,ms.rf3a$q2.5$pF1)
surv_veg_q97.5<-c(ms.rf3a$q97.5$sV0,ms.rf3a$q97.5$sV1)
surv_rep_q97.5<-c(ms.rf3a$q97.5$sF0,ms.rf3a$q97.5$sF1)
trans_vr_q97.5<-c(ms.rf3a$q97.5$fV0,ms.rf3a$q97.5$fV1)
trans_rv_q97.5<-c(ms.rf3a$q97.5$fF0,ms.rf3a$q97.5$fF1)
emer_veg_q97.5<-c(ms.rf3a$q97.5$pV0,ms.rf3a$q97.5$pV1)
emer_rep_q97.5<-c(ms.rf3a$q97.5$pF0,ms.rf3a$q97.5$pF1)

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

