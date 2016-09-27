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
#Try with 3 stages- as in example 9.7 page 307 in Kery & Schaub, and using probabilities of flowering and dormancy (instend of transitions and p)
sink("ms-ranef3stages.jags")
cat("

model {
    # -------------------------------------------------
    # Parameters:
    # sV: survival probability for nonreproductive plants
    # sF: survival probability for reproductive plants
    # sU: survival probability for dormant plants
    # fV: reproduction probability for vegetative plants
    # fF: reproduction probability for reproductive plants
    # fU: reproduction probability for dormant plants
    # dV: dormancy probability for nonreproductive plants
    # dF: dormancy probability for reproductive plants
    # dU: dormancy probability for dormant plants
    # -------------------------------------------------
    # States (S):
    # 1 vegetation
    # 2 reproductive
    # 3 dormant
    # 4 dead
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
    logit(sU[i,t]) <- eta.sU[group[i],t]
    logit(fV[i,t]) <- eta.fV[group[i],t]
    logit(fF[i,t]) <- eta.fF[group[i],t]
    logit(fU[i,t]) <- eta.fU[group[i],t]
    logit(dV[i,t]) <- eta.dV[group[i],t]
    logit(dF[i,t]) <- eta.dF[group[i],t]
    logit(dU[i,t]) <- eta.dU[group[i],t]
    }#t
    }#i
    for (g in 1:2){#group  
    for (t in 1:(n.occasions-1)){
    eta.sV[g,t]  <-mu.sV[g]+beta.sV[g]*x[g,t]+beta1.sV[g]*x1[g,t]+epsilon.sV[g,t]
    eta.sF[g,t]  <-mu.sF[g]+beta.sF[g]*x[g,t]+beta1.sF[g]*x1[g,t]+epsilon.sF[g,t]
    eta.sU[g,t]  <-mu.sU[g]+beta.sU[g]*x[g,t]+beta1.sU[g]*x1[g,t]+epsilon.sU[g,t]
    eta.fV[g,t]  <-mu.fV[g]+beta.fV[g]*x[g,t]+beta1.fV[g]*x1[g,t]+epsilon.fV[g,t]
    eta.fF[g,t]  <-mu.fF[g]+beta.fF[g]*x[g,t]+beta1.fF[g]*x1[g,t]+epsilon.fF[g,t]
    eta.fU[g,t]  <-mu.fU[g]+beta.fU[g]*x[g,t]+beta1.fU[g]*x1[g,t]+epsilon.fU[g,t]
    eta.dV[g,t]  <-mu.dV[g]+beta.dV[g]*x[g,t]+beta1.dV[g]*x1[g,t]+epsilon.dV[g,t]
    eta.dF[g,t]  <-mu.dF[g]+beta.dF[g]*x[g,t]+beta1.dF[g]*x1[g,t]+epsilon.dF[g,t]
    eta.dU[g,t]  <-mu.dU[g]+beta.dU[g]*x[g,t]+beta1.dU[g]*x1[g,t]+epsilon.dU[g,t]
    epsilon.sV[g,t]~dnorm(0,tau.sV[g])#could move all the mean stuff to this part, might help model mix better? hierarchical centering
    epsilon.sF[g,t]~dnorm(0,tau.sF[g])
    epsilon.sU[g,t]~dnorm(0,tau.sU[g])
    epsilon.fV[g,t]~dnorm(0,tau.fV[g])
    epsilon.fF[g,t]~dnorm(0,tau.fF[g])
    epsilon.fU[g,t]~dnorm(0,tau.fU[g])
    epsilon.dV[g,t]~dnorm(0,tau.dV[g])
    epsilon.dF[g,t]~dnorm(0,tau.dF[g])
    epsilon.dU[g,t]~dnorm(0,tau.dU[g])
      }#t
    mean.sV[g]~dunif(0,1)# Priors for mean group-specific survival for veg plants
    mu.sV[g]<-log(mean.sV[g]/(1-mean.sV[g]))  
    mean.sF[g]~dunif(0,1)# Priors for mean group-specific survivalsurvival for rep plants
    mu.sF[g]<-log(mean.sF[g]/(1-mean.sF[g]))    
    mean.sU[g]~dunif(0,1)# Priors for mean group-specific survivalsurvival for rep plants
    mu.sU[g]<-log(mean.sU[g]/(1-mean.sU[g])) 
    mean.fV[g]~dunif(0,1)# Priors for mean group-specific transition from veg to rep plants
    mu.fV[g]<-log(mean.fV[g]/(1-mean.fV[g]))       
    mean.fF[g]~dunif(0,1)# Priors for mean group-specific transition from rep to veg plants
    mu.fF[g]<-log(mean.fF[g]/(1-mean.fF[g]))    
    mean.fU[g]~dunif(0,1)# Priors for mean group-specific transition from rep to veg plants
    mu.fU[g]<-log(mean.fU[g]/(1-mean.fU[g])) 
    mean.dV[g]~dunif(0,1)# Priors for mean group-specific transition from veg to rep plants
    mu.dV[g]<-log(mean.dV[g]/(1-mean.dV[g]))       
    mean.dF[g]~dunif(0,1)# Priors for mean group-specific transition from rep to veg plants
    mu.dF[g]<-log(mean.dF[g]/(1-mean.dF[g]))    
    mean.dU[g]~dunif(0,1)# Priors for mean group-specific transition from rep to veg plants
    mu.dU[g]<-log(mean.dU[g]/(1-mean.dU[g])) 
    sigma.sV[g]~dunif(0,10)#temporal variance for veg plants survival
    tau.sV[g]<-pow(sigma.sV[g],-2)
    sigma.sV2[g]<-pow(sigma.sV[g],2)
    sigma.sF[g]~dunif(0,10)#temporal variance for rep plants survival
    tau.sF[g]<-pow(sigma.sF[g],-2)
    sigma.sF2[g]<-pow(sigma.sF[g],2)
    sigma.sU[g]~dunif(0,10)#temporal variance for rep plants survival
    tau.sU[g]<-pow(sigma.sU[g],-2)
    sigma.sU2[g]<-pow(sigma.sU[g],2)
    sigma.fV[g]~dunif(0,10)#temporal variance for veg plants transition to rep
    tau.fV[g]<-pow(sigma.fV[g],-2)
    sigma.fV2[g]<-pow(sigma.fV[g],2)
    sigma.fF[g]~dunif(0,10)#temporal variance for rep plants transition to veg
    tau.fF[g]<-pow(sigma.fF[g],-2)
    sigma.fF2[g]<-pow(sigma.fF[g],2)
    sigma.fU[g]~dunif(0,10)#temporal variance for rep plants transition to veg
    tau.fU[g]<-pow(sigma.fU[g],-2)
    sigma.fU2[g]<-pow(sigma.fU[g],2)
    sigma.dV[g]~dunif(0,10)#temporal variance of dormancy for veg plants 
    tau.dV[g]<-pow(sigma.dV[g],-2)
    sigma.dV2[g]<-pow(sigma.dV[g],2)
    sigma.dF[g]~dunif(0,10)#temporal variance of dormancy for rep plants
    tau.dF[g]<-pow(sigma.dF[g],-2)
    sigma.dF2[g]<-pow(sigma.dF[g],2)
    sigma.dU[g]~dunif(0,10)#temporal variance of dormancy for dormant plants
    tau.dU[g]<-pow(sigma.dU[g],-2)
    sigma.dU2[g]<-pow(sigma.dU[g],2)
    beta.sV[g]~dnorm(0,0.001)I(-10,10)
    beta.sF[g]~dnorm(0,0.001)I(-10,10)
    beta.sU[g]~dnorm(0,0.001)I(-10,10)
    beta.fV[g]~dnorm(0,0.001)I(-10,10)
    beta.fF[g]~dnorm(0,0.001)I(-10,10)
    beta.fU[g]~dnorm(0,0.001)I(-10,10)
    beta.dV[g]~dnorm(0,0.001)I(-10,10)
    beta.dF[g]~dnorm(0,0.001)I(-10,10)
    beta.dU[g]~dnorm(0,0.001)I(-10,10)
    beta1.sV[g]~dnorm(0,0.001)I(-10,10)
    beta1.sF[g]~dnorm(0,0.001)I(-10,10)
    beta1.sU[g]~dnorm(0,0.001)I(-10,10)
    beta1.fV[g]~dnorm(0,0.001)I(-10,10)
    beta1.fF[g]~dnorm(0,0.001)I(-10,10)
    beta1.fU[g]~dnorm(0,0.001)I(-10,10)
    beta1.dV[g]~dnorm(0,0.001)I(-10,10)
    beta1.dF[g]~dnorm(0,0.001)I(-10,10)
    beta1.dU[g]~dnorm(0,0.001)I(-10,10)
    #calculate probabilities of 2 different groups and 2 different time periods (=4 group-time periods) to look at later... 
    #survival
    logit(sV0[g])<- mu.sV[g]
    logit(sV1[g])<- mu.sV[g] + beta.sV[g]
    logit(sF0[g])<- mu.sF[g]
    logit(sF1[g])<- mu.sF[g] + beta.sF[g]
    logit(sU0[g])<- mu.sU[g]
    logit(sU1[g])<- mu.sU[g] + beta.sU[g]
    #reproduction
    logit(fV0[g])<- mu.fV[g]
    logit(fV1[g])<- mu.fV[g] + beta.fV[g]
    logit(fF0[g])<- mu.fF[g]
    logit(fF1[g])<- mu.fF[g] + beta.fF[g]
    logit(fU0[g])<- mu.fU[g]
    logit(fU1[g])<- mu.fU[g] + beta.fU[g]
    #dormancy
    logit(dV0[g])<- mu.dV[g]
    logit(dV1[g])<- mu.dV[g] + beta.dV[g]
    logit(dF0[g])<- mu.dF[g]
    logit(dF1[g])<- mu.dF[g] + beta.dF[g]
    logit(dU0[g])<- mu.dU[g]
    logit(dU1[g])<- mu.dU[g] + beta.dU[g]
    }#g
    
    # Define state-transition and observation matrices
    for (i in 1:nind){  
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    ps[1,i,t,1] <- sV[i,t] * (1-fV[i,t])*(1-dV[i,t])
    ps[1,i,t,2] <- sV[i,t] * fV[i,t]*(1-dV[i,t])
    ps[1,i,t,3] <- sV[i,t]*dV[i,t]
    ps[1,i,t,4] <- 1-sV[i,t]

    ps[2,i,t,1] <- sF[i,t] * (1-fF[i,t])*(1-dF[i,t])
    ps[2,i,t,2] <- sF[i,t] * fF[i,t]*(1-dF[i,t])
    ps[2,i,t,3] <- sF[i,t]*dF[i,t]
    ps[2,i,t,4] <- 1-sF[i,t]

    ps[3,i,t,1] <- sU[i,t] * (1-fU[i,t])*(1-dU[i,t])
    ps[3,i,t,2] <- sU[i,t] * fU[i,t]*(1-dU[i,t])
    ps[3,i,t,3] <- sU[i,t]* dU[i,t]
    ps[3,i,t,4] <- 1-sU[i,t]

    ps[4,i,t,1] <- 0
    ps[4,i,t,2] <- 0
    ps[4,i,t,3] <- 0
    ps[4,i,t,4] <- 1

    # Define probabilities of O(t) given S(t) ##first coumn is observed state, last is true state
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
inits<-function(){list(mean.sV=c(.5,.5), mean.sF=c(.5,.5),mean.sU=c(.5,.5),mean.fV=c(.5,.5),mean.fF=c(.5,.5),mean.fU=c(.5,.5),mean.dV=c(.5,.5),mean.dF=c(.5,.5),mean.dU=c(.5,.5),sigma.sV=c(1,1),sigma.sF=c(1,1),sigma.sU=c(1,1),sigma.fV=c(1,1),sigma.fF=c(1,1),sigma.fU=c(1,1),sigma.dV=c(1,1),sigma.dF=c(1,1),sigma.dU=c(1,1),z = zst1)}
# Parameters monitored
parameters <- c("mean.sV","mean.sF","mean.sU", "mean.fV","mean.fF","mean.fU","mean.dV","mean.dF","mean.dU","beta.sV","beta.sF", "beta.sU","beta.fV","beta.fF","beta.fU","beta.dV", "beta.dF", "beta.dU","sV0","sV1","sF0","sF1","sU0","sU1","fV0","fV1","fF0","fF1","fU0","fU1","dV0","dV1","dF0","dF1","dU0","dU1","beta1.sV","beta1.sF","beta1.sU", "beta1.fV","beta1.fF","beta1.fU","beta1.dV", "beta1.dF","beta1.dU","mu.sV","mu.sF","mu.sU","mu.fV","mu.fF","mu.fU","mu.dV","mu.dF","mu.dU","sigma.sV2","sigma.sF2","sigma.sU2","sigma.fV2","sigma.fF2","sigma.fU2","sigma.dV2","sigma.dF2","sigma.dU2")

# MCMC settings
ni <- 10000
nt <- 3
nb <- 500
nc <- 3

# Call JAGS from R #
#complex model
ms3.rf <- jags(jags.data, inits, parameters, "ms-ranef3stages.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T)
print(ms3.rf, digits=3)

###save all samples
mod.samples_3stage<- as.data.frame(do.call("rbind", ms3.rf$samples))
write.csv(mod.samples_3stage,"msmod_samples_3stage.csv",row.names=T)
mod.samples_3stage<-read.csv("msmod_samples_complex.csv", header=T)

##Figures
#2x2table for each vital rate with first column control, second column logged
#first row before logging, second row after logging
surv_veg<-as.data.frame(rbind(ms3.rf$mean$sV0,ms3.rf$mean$sV1))
surv_rep<-as.data.frame(rbind(ms3.rf$mean$sF0,ms3.rf$mean$sF1))
surv_dorm<-as.data.frame(rbind(ms3.rf$mean$sU0,ms3.rf$mean$sU1))

prep_veg<-as.data.frame(rbind(ms3.rf$mean$fV0,ms3.rf$mean$fV1))
prep_rep<-as.data.frame(rbind(ms3.rf$mean$fF0,ms3.rf$mean$fF1))
prep_dorm<-as.data.frame(rbind(ms3.rf$mean$fU0,ms3.rf$mean$fU1))

pdorm_veg<-as.data.frame(rbind(ms3.rf$mean$dV0,ms3.rf$mean$dV1))
pdorm_rep<-as.data.frame(rbind(ms3.rf$mean$dF0,ms3.rf$mean$dF1))
pdorm_dorm<-as.data.frame(rbind(ms3.rf$mean$dU0,ms3.rf$mean$dU1))

surv_veg_q2.5<-c(ms3.rf$q2.5$sV0,ms3.rf$q2.5$sV1)
surv_rep_q2.5<-c(ms3.rf$q2.5$sF0,ms3.rf$q2.5$sF1)
surv_dorm_q2.5<-c(ms3.rf$q2.5$sU0,ms3.rf$q2.5$sU1)

prep_veg_q2.5<-c(ms3.rf$q2.5$fV0,ms3.rf$q2.5$fV1)
prep_rep_q2.5<-c(ms3.rf$q2.5$fF0,ms3.rf$q2.5$fF1)
prep_dorm_q2.5<-c(ms3.rf$q2.5$fU0,ms3.rf$q2.5$fU1)

pdorm_veg_q2.5<-c(ms3.rf$q2.5$dV0,ms3.rf$q2.5$dV1)
pdorm_rep_q2.5<-c(ms3.rf$q2.5$dF0,ms3.rf$q2.5$dF1)
pdorm_rep_q2.5<-c(ms3.rf$q2.5$dU0,ms3.rf$q2.5$dU1)

surv_veg_q97.5<-c(ms3.rf$q97.5$sV0,ms3.rf$q97.5$sV1)
surv_rep_q97.5<-c(ms3.rf$q97.5$sF0,ms3.rf$q97.5$sF1)
surv_dorm_q97.5<-c(ms3.rf$q97.5$sU0,ms3.rf$q97.5$sU1)

prep_veg_q97.5<-c(ms3.rf$q97.5$fV0,ms3.rf$q97.5$fV1)
prep_rep_q97.5<-c(ms3.rf$q97.5$fF0,ms3.rf$q97.5$fF1)
prep_dorm_q97.5<-c(ms3.rf$q97.5$fU0,ms3.rf$q97.5$fU1)

pdorm_veg_q97.5<-c(ms3.rf$q97.5$dV0,ms3.rf$q97.5$dV1)
pdorm_rep_q97.5<-c(ms3.rf$q97.5$dF0,ms3.rf$q97.5$dF1)
pdorm_rep_q97.5<-c(ms3.rf$q97.5$dU0,ms3.rf$q97.5$dU1)
##Table of model summary, with betas
mod.sum<-ms3.rf$summary
write.csv(mod.sum,"isotria3stagemodsum.csv", row.names=T)
mod.table<-round(subset(mod.sum, select=c("mean","sd","Rhat", "n.eff")), digits=3)
mod.table2<-rbind(mod.table[which(substr(rownames(mod.table),1,2)=="mu"),],mod.table[which(substr(rownames(mod.table),1,4)=="beta"),])
write.csv(mod.table2,"isotria3stage_modtable.csv")
