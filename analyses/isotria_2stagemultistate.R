#Multistate model for Isotria medioloides Alton NH population
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

##Figures
#2x2table for each vital rate with first column control, second column logged
#first row before logging, second row after logging
surv_veg<-as.data.frame(rbind(ms.rf3a$mean$phiA0,ms.rf3a$mean$phiA1))
surv_rep<-as.data.frame(rbind(ms.rf3a$mean$phiB0,ms.rf3a$mean$phiB1))
trans_vr<-as.data.frame(rbind(ms.rf3a$mean$psiA0,ms.rf3a$mean$psiA1))
trans_rv<-as.data.frame(rbind(ms.rf3a$mean$psiB0,ms.rf3a$mean$psiB1))
emer_veg<-as.data.frame(rbind(ms.rf3a$mean$pA0,ms.rf3a$mean$pA1))
emer_rep<-as.data.frame(rbind(ms.rf3a$mean$pB0,ms.rf3a$mean$pB1))

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

