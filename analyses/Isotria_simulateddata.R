#Isotria multistate model
#Simulate data with known parameters to test the ability of my model to recover them
rm(list=ls()) 
options(stringsAsFactors=FALSE)

library(rjags)
library(jagsUI)
library(lattice)
library(coda)
library(boot)
library(dplyr)
#start with simple dataset with no variation between groups x and y
#first, define mean survival, flowering probability, and emergence probability, 
#as well as the number of occassions, states, observations, and marked individuals
sV<-0.8#survival=phi
sF<-0.85
sUV<-0.7
sUF<-0.75
fV<-0.2#flowering probability
fF<- 0.5
fUV<-0.2
fUF<- 0.5
dV<-0.3#probability of dormancy
dF<-0.2
dUV<-0.6
dUF<-0.4
n.occasions<-31
n.states<-5
n.obs<-3
marked<-matrix(NA,ncol=n.states,nrow=n.occasions)
marked[,1]<-c(100,rep(10,n.occasions-1))
marked[,2]<-c(100,rep(10,n.occasions-1))
marked[,3]<-rep(0,n.occasions)
marked[,4]<-rep(0,n.occasions)
marked[,5]<-rep(0,n.occasions)

#State process matrix.This has 4 -dimensions: 
#d1: state of departure
#d2: state or arrival
#d3: individual
#d4: time
totrel<-sum(marked)*(n.occasions-1)
STATE<-array(NA,dim=c(n.states,n.states,totrel,n.occasions-1))
for(i in 1:totrel){
  for(t in 1:(n.occasions-1)){
    STATE[,,i,t]<-matrix(c(
    sV*(1-fV)*(1-dV),sV*fV*(1-dV),sV*dV,0,1-sV,
    sF*(1-fF)*(1-dF),sF*fF*(1-dF),0,sF*dF,1-sF,
    sUV*(1-fUV)*(1-dUV),sUV*fUV*(1-dUV),sUV*dUV,0,1-sUV,
   sUF*(1-fUF)*(1-dUF),sUF*fUF*(1-dUF),0,sUF*dUF,1-sUF,
    0,0,0,0,1),nrow=n.states, byrow=TRUE)
  }#t
}#i
#Observation process matrix
# Define probabilities of O(t) given S(t) ##first column is observed state, last is true state
OBS<-array(NA,dim=c(n.states,n.obs,totrel,n.occasions-1))#replaced NA with 0 becuase got error "
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    OBS[,,i,t]<-matrix(c(
      1,0,0,
      0,1,0,
      0,0,1,
      0,0,1,
      0,0,1), nrow=n.states,byrow = TRUE)
  }#t
}#i
#function to simulate multistate capture-recapture data
simul.ms<-function(STATE,OBS,marked,unobservable=NA){
  #unobservable: number of states that are unobservable. in this
  n.occasions<-dim(STATE)[4]+1
  CH<-CH.TRUE<-matrix(NA,ncol=n.occasions,nrow=sum(marked))
  #define a vector with the occasion of marking
  mark.occ<-matrix(0,ncol=n.occasions, nrow=sum(marked))
  g<-colSums(marked)
  for (s in 1:dim(STATE)[1]){
    if (g[s]==0) {next}#to avoid error message if nothing to replace
    mark.occ[(cumsum(g[1:s])-g[s]+1) [s] :cumsum(g[1:s])[s],s] <-rep(1:n.occasions, marked[1:n.occasions,s])
  }#s
  for(i in 1:sum(marked)){
    for(s in 1:dim(STATE)[1]){
      if (mark.occ[i,s]==0) {next}
      first<-mark.occ[i,s]
      CH[i,first]<-s
      CH.TRUE[i,first]<-s
    }#s
    for(t in (first+1):n.occasions){
      #multinomial trials for state transitions
      if (first==n.occasions) {next}
      state<-which(rmultinom(1,1,STATE[CH.TRUE[i,t-1],,i,t-1])==1)
      CH.TRUE[i,t]<-state
      #multinomial trials for observation process
      event<-which(rmultinom(1,1,OBS[CH.TRUE[i,t],,i,t-1])==1)
      CH[i,t]<-event
    }#t
  }#i
  #replace the NA and the highest state number (dead) in the file by 0
  #CH[is.na(CH)]<-0
  CH[CH==dim(STATE)[1]]<-0
  CH[CH==unobservable]<-0
  id<-numeric(0)
  for(i in 1:dim(CH)[1]){
    z<-min(which(CH[i,]!=0))
    ifelse(z==dim(CH)[2], id<-c(id,i), id<-c(id))
  }
  return(list(CH=CH[-id,],CH.TRUE=CH.TRUE[-id,]))
  #CH: captire-histories to be used
  #CH.TRUE: capture histories with perfect observation
}
sim<-simul.ms(STATE,OBS,marked)
CH<-sim$CH
#Replace first observations that are 0 with NAs
CH[which(CH[,1]==0),1]<-NA
#compute vector with occasion of first capture
get.first<-function(x) min(which(x!=0))
CH2<-rbind(CH,CH)
f<-apply(CH2,1,get.first)
#recode CH matrix since 0 is not allowed in WinBUGS
rCH<-CH2
rCH[rCH==0]<-3
rCH[rCH==4]<-3
n=dim(rCH)[1]#number plants in sulated dataset
###now use my 4-stage model to estimate parameters of these data
group<-c(rep(1,n/2),rep(2,n/2))#group x=control=1, group y=logged=2#assign half plants to gorup 1 and half to group 2
logged_x<-c(rep(0,30))#logging never occurred for group x- this is one less than the number of years 
logged_y<-c(rep(0,12),rep(1,18))#for group y, logging occurred between 1997 and 1998, so first 13 years= no logging treatment
logged_yrs<-c(rep(0,13),seq(1:17))#add decay time- years since logging this is Beta1, this is x1, add a year after the initial year of logging to allow for a bump up
logged_yrs2<-rbind(logged_yrs,logged_yrs)
time_log<-rbind(logged_y,logged_y)#to test difference between groups before and after logging

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
known.state.ms <- function(ch){##removes 3s and replaces them with NAs, and replaces first observation with NA
  state <- ch
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,]<3))#
    state[i,n1] <- NA
  }
  state[state==3] <- NA
  return(state)
}

#give starting values of 3 or 4 to all unknown states
ms.init.z <- function(ch, f){#ms.init.z gives starting values of 3 or 4 to all unknown states
  zstart<-c()
  #matrix(data=NA,nrow=dim(ch)[1], ncol=dim(ch)[2]))
  for(i in 1:dim(ch)[1]){
    z<-ch[i,]
    v <- which(z>=3)# occurences that are unknown states get 3s by default (already coded)
    w<-which(z<3)
    x<-which(z==2)
    z[-v] <- NA#observed states get NA
    if(length(x)==0){z[w]<-NA}
    else if(length(which(z==3))> 0 & length(x)>0 & !('2' %in% ch[i,which(z==3)-1])){z[w]<-NA}
    else if(length(which(ch[i,]>=3))==0){
      z[w]<-NA
    } else if (length(which(z==3))> 0 & length(x)>0 & '2' %in% ch[i,which(z==3)-1])
    {
      #if(length(x)>=1){
        temp<-ch[i,]
        temp[which(is.na(temp))]<-0
        for (j in 1:sum(rle(temp)$values == 2)){
          if(is.na(rle(temp)$values[(which(rle(temp)$values==2)+1)[j]]==3)){next}
          if(rle(temp)$values[(which(rle(temp)$values==2)+1)[j]]==3){
          fourstart<-sum(rle(temp)$lengths[1:which(rle(temp)$values==2)[j]])+1
          threestofours<-which(rle(temp)$values==3)
          threestofours<-threestofours[threestofours>which(rle(temp)$values==2)[j]]
          fourend<-sum(rle(temp)$lengths[1:threestofours])
          z[fourstart:fourend]<-4}
          else {next}
          }
        }
    zstart<-rbind(zstart,z)
  }
      return(zstart)
}

#for (i in 1:length(sim_ds){
  
zst=ms.init.z(rCH,f)
# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = known.state.ms(rCH),group=group, x=time_log, x1=logged_yrs2)

# Initial values
inits<-function(){list(mean.sV=c(.5,.5), mean.sF=c(.5,.5),mean.sUV=c(.5,.5),mean.sUF=c(.5,.5),mean.fV=c(.5,.5),mean.fF=c(.5,.5),mean.fUV=c(.5,.5),mean.fUF=c(.5,.5),mean.dV=c(.5,.5),mean.dF=c(.5,.5),mean.dUV=c(.5,.5),mean.dUF=c(.5,.5),sigma.sV=c(1,1),sigma.sF=c(1,1),sigma.sUV=c(1,1),sigma.sUF=c(1,1),sigma.fV=c(1,1),sigma.fF=c(1,1),sigma.fUV=c(1,1),sigma.fUF=c(1,1),sigma.dV=c(1,1),sigma.dF=c(1,1),sigma.dUV=c(1,1),sigma.dUF=c(1,1),z = zst)}
# Parameters monitored
parameters <- c("mean.sV","mean.sF","mean.sUV","mean.sUF", "mean.fV","mean.fF","mean.fUV","mean.fUF","mean.dV","mean.dF","mean.dUV","mean.dUF","beta.sV","beta.sF", "beta.sUV","beta.sUF","beta.fV","beta.fF","beta.fUV","beta.fUF","beta.dV", "beta.dF", "beta.dUV","beta.dUF","sV0","sV1","sF0","sF1","sUV0","sUF0","sUV1","sUF1","fV0","fV1","fF0","fF1","fUV0","fUV1","fUF0","fUF1","dV0","dV1","dF0","dF1","dUV0","dUF0","dUV1","dUF1","beta1.sV","beta1.sF","beta1.sUV","beta1.sUF","beta1.fV","beta1.fF","beta1.fUV","beta1.fUF","beta1.dV", "beta1.dF","beta1.dUV","beta1.dUF","mu.sV","mu.sF","mu.sUV","mu.sUF","mu.fV","mu.fF","mu.fUV","mu.fUF","mu.dV","mu.dF","mu.dUV","mu.dUF","sigma.sV2","sigma.sF2","sigma.sUV2","sigma.sUF2","sigma.fV2","sigma.fF2","sigma.fUV2","sigma.fUF2","sigma.dV2","sigma.dF2","sigma.dUV2","sigma.dUF2")

# MCMC settings
ni <- 5000
nt <- 5
nb <- 500
nc <- 3

# Call JAGS from R #
#complex model
ms4.rf <- jags(jags.data, inits, parameters, "ms-ranef4stages.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T)
print(ms4.rf, digits=3)
###save all samples
mod.samples_4stage<- as.data.frame(do.call("rbind", ms4.rf$samples))
write.csv(mod.samples_4stage,"msmod_samples_4stagesim.csv",row.names=T)

##Table of model summary, with betas
mod.sum<-ms4.rf$summary
write.csv(mod.sum,"isotria4stagemodsumsim.csv", row.names=T)
mod.table<-round(subset(mod.sum, select=c("mean","sd","Rhat", "n.eff")), digits=3)
mod.table2<-rbind(mod.table[which(substr(rownames(mod.table),1,2)=="mu"),],mod.table[which(substr(rownames(mod.table),1,4)=="beta"),])
write.csv(mod.table2,"isotria4stage_modtablesim.csv")
#simulate 100 datasets, then fit the model to each dataset let run for 5000 iterations for each.
#Andy writes a function that saves output. 
#Results after comparing 4-stage model with 2-stage model (below):
#2-stage more accurately recovers survival; 4-stage model more accurately recovers d and f
#Survival estimates for 4-stage are high

#In meeting 9 Feb 2017, Andy suggests fitting model to simulated data with 2 survival parameters: 1 for V/UV  & one for F/UF
#Here is code for that model
sink("ms-ranef4stages_2surv.jags")
cat("
    
    model {
    # -------------------------------------------------
    # Parameters:
    # sV: survival probability for nonreproductive plants and for for dormant plants that were vegetative above ground
    # sF: survival probability for reproductive plants and for dormant plants that were reproductive above ground
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
    ps[3,i,t,1] <- sV[i,t] * (1-fUV[i,t])*(1-dUV[i,t])
    ps[3,i,t,2] <- sV[i,t] * fUV[i,t]*(1-dUV[i,t])
    ps[3,i,t,3] <- sV[i,t]* dUV[i,t]
    ps[3,i,t,4] <- 0
    ps[3,i,t,5] <- 1-sV[i,t]
    ps[4,i,t,1] <- sF[i,t] * (1-fUF[i,t])*(1-dUF[i,t])
    ps[4,i,t,2] <- sF[i,t] * fUF[i,t]*(1-dUF[i,t])
    ps[4,i,t,3] <- 0
    ps[4,i,t,4] <- sF[i,t]* dUF[i,t]
    ps[4,i,t,5] <- 1-sF[i,t]
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
known.state.ms <- function(ch){##removes 3s and replaces them with NAs, and replaces first observation with NA
  state <- ch
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,]<3))#
    state[i,n1] <- NA
  }
  state[state==3] <- NA
  return(state)
}

#give starting values of 3 or 4 to all unknown states
ms.init.z <- function(ch, f){#ms.init.z gives starting values of 3 or 4 to all unknown states
  zstart<-c()
  #matrix(data=NA,nrow=dim(ch)[1], ncol=dim(ch)[2]))
  for(i in 1:dim(ch)[1]){
    z<-ch[i,]
    v <- which(z>=3)# occurences that are unknown states get 3s by default (already coded)
    w<-which(z<3)
    x<-which(z==2)
    z[-v] <- NA#observed states get NA
    if(length(x)==0){z[w]<-NA}
    else if(length(which(z==3))> 0 & length(x)>0 & !('2' %in% ch[i,which(z==3)-1])){z[w]<-NA}
    else if(length(which(ch[i,]>=3))==0){
      z[w]<-NA
    } else if (length(which(z==3))> 0 & length(x)>0 & '2' %in% ch[i,which(z==3)-1])
    {
      #if(length(x)>=1){
      temp<-ch[i,]
      temp[which(is.na(temp))]<-0
      for (j in 1:sum(rle(temp)$values == 2)){
        if(is.na(rle(temp)$values[(which(rle(temp)$values==2)+1)[j]]==3)){next}
        if(rle(temp)$values[(which(rle(temp)$values==2)+1)[j]]==3){
          fourstart<-sum(rle(temp)$lengths[1:which(rle(temp)$values==2)[j]])+1
          threestofours<-which(rle(temp)$values==3)
          threestofours<-threestofours[threestofours>which(rle(temp)$values==2)[j]]
          fourend<-sum(rle(temp)$lengths[1:threestofours])
          z[fourstart:fourend]<-4}
        else {next}
      }
    }
    zstart<-rbind(zstart,z)
  }
  return(zstart)
}



#for (i in 1:length(sim_ds){

zst=ms.init.z(rCH,f)
# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = known.state.ms(rCH),group=group, x=time_log, x1=logged_yrs2)

# Initial values
inits<-function(){list(mean.sV=c(.5,.5), mean.sF=c(.5,.5),mean.fV=c(.5,.5),mean.fF=c(.5,.5),mean.fUV=c(.5,.5),mean.fUF=c(.5,.5),mean.dV=c(.5,.5),mean.dF=c(.5,.5),mean.dUV=c(.5,.5),mean.dUF=c(.5,.5),sigma.sV=c(1,1),sigma.sF=c(1,1),sigma.fV=c(1,1),sigma.fF=c(1,1),sigma.fUV=c(1,1),sigma.fUF=c(1,1),sigma.dV=c(1,1),sigma.dF=c(1,1),sigma.dUV=c(1,1),sigma.dUF=c(1,1),z = zst)}
# Parameters monitored
parameters <- c("mean.sV","mean.sF", "mean.fV","mean.fF","mean.fUV","mean.fUF","mean.dV","mean.dF","mean.dUV","mean.dUF","beta.sV","beta.sF", "beta.fV","beta.fF","beta.fUV","beta.fUF","beta.dV", "beta.dF", "beta.dUV","beta.dUF","sV0","sV1","sF0","sF1","fV0","fV1","fF0","fF1","fUV0","fUV1","fUF0","fUF1","dV0","dV1","dF0","dF1","dUV0","dUF0","dUV1","dUF1","beta1.sV","beta1.sF","beta1.fV","beta1.fF","beta1.fUV","beta1.fUF","beta1.dV", "beta1.dF","beta1.dUV","beta1.dUF","mu.sV","mu.sF","mu.fV","mu.fF","mu.fUV","mu.fUF","mu.dV","mu.dF","mu.dUV","mu.dUF","sigma.sV2","sigma.sF2","sigma.fV2","sigma.fF2","sigma.fUV2","sigma.fUF2","sigma.dV2","sigma.dF2","sigma.dUV2","sigma.dUF2")

# MCMC settings
ni <- 5000
nt <- 5
nb <- 500
nc <- 3

# Call JAGS from R #
#complex model
ms4_2surv.rf <- jags(jags.data, inits, parameters, "ms-ranef4stages_2surv.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T)
print(ms4_2surv.rf, digits=3)
###save all samples
mod.samples_4stage2surv<- as.data.frame(do.call("rbind", ms4_2surv.rf$samples))
write.csv(mod.samples_4stage2surv,"msmod_samples_4stage_2survsim.csv",row.names=T)

#Now try fitting 2-stage model to simulated data and see how it does:
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
    # fF: transition probability from reproductive to nonfruiting
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
    logit(fF[i,t]) <- eta.fF[group[i],t]
    logit(pV[i,t]) <- eta.pV[group[i],t]
    logit(pF[i,t]) <- eta.pF[group[i],t]    
    }#t
    }#i
    for (g in 1:2){#group  
    for (t in 1:(n.occasions-1)){
    eta.sV[g,t]  <-mu.sV[g]+beta.sV[g]*x[g,t]+beta1.sV[g]*x1[g,t]+epsilon.sV[g,t]
    eta.sF[g,t]  <-mu.sF[g]+beta.sF[g]*x[g,t]+beta1.sF[g]*x1[g,t]+epsilon.sF[g,t]
    eta.fV[g,t]  <-mu.fV[g]+beta.fV[g]*x[g,t]+beta1.fV[g]*x1[g,t]+epsilon.fV[g,t]
    eta.fF[g,t]  <-mu.fF[g]+beta.fF[g]*x[g,t]+beta1.fF[g]*x1[g,t]+epsilon.fF[g,t]
    eta.pV[g,t]  <-mu.pV[g]+beta.pV[g]*x[g,t]+beta1.pV[g]*x1[g,t]+epsilon.pV[g,t]
    eta.pF[g,t]  <-mu.pF[g]+beta.pF[g]*x[g,t]+beta1.pF[g]*x1[g,t]+epsilon.pF[g,t]
    epsilon.sV[g,t]~dnorm(0,tau.sV[g])#could move all the mean stuff to this part, might help model mix better? hierarchical centering
    epsilon.sF[g,t]~dnorm(0,tau.sF[g])
    epsilon.fV[g,t]~dnorm(0,tau.fV[g])
    epsilon.fF[g,t]~dnorm(0,tau.fF[g])
    epsilon.pV[g,t]~dnorm(0,tau.pV[g])
    epsilon.pF[g,t]~dnorm(0,tau.pF[g])
    }#t
    mean.sV[g]~dunif(0,1)# Priors for mean group-specific survival for veg plants
    mu.sV[g]<-log(mean.sV[g]/(1-mean.sV[g]))  
    mean.sF[g]~dunif(0,1)# Priors for mean group-specific survivalsurvival for rep plants
    mu.sF[g]<-log(mean.sF[g]/(1-mean.sF[g]))    
    mean.fV[g]~dunif(0,1)# Priors for mean group-specific transition from veg to rep plants
    mu.fV[g]<-log(mean.fV[g]/(1-mean.fV[g]))       
    mean.fF[g]~dunif(0,1)# Priors for mean group-specific transition from rep to veg plants
    mu.fF[g]<-log(mean.fF[g]/(1-mean.fF[g]))    
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
    sigma.fF[g]~dunif(0,10)#temporal variance for rep plants transition to veg
    tau.fF[g]<-pow(sigma.fF[g],-2)
    sigma.fF2[g]<-pow(sigma.fF[g],2)
    sigma.pV[g]~dunif(0,10)#temporal variance for veg plants recapture (emergence)
    tau.pV[g]<-pow(sigma.pV[g],-2)
    sigma.pV2[g]<-pow(sigma.pV[g],2)
    sigma.pF[g]~dunif(0,10)#temporal variance for rep plants recapture (emergence)
    tau.pF[g]<-pow(sigma.pF[g],-2)
    sigma.pF2[g]<-pow(sigma.pF[g],2)
    beta.sV[g]~dnorm(0,0.001)I(-10,10)
    beta.sF[g]~dnorm(0,0.001)I(-10,10)
    beta.fV[g]~dnorm(0,0.001)I(-10,10)
    beta.fF[g]~dnorm(0,0.001)I(-10,10)
    beta.pV[g]~dnorm(0,0.001)I(-10,10)
    beta.pF[g]~dnorm(0,0.001)I(-10,10)
    beta1.sV[g]~dnorm(0,0.001)I(-10,10)
    beta1.sF[g]~dnorm(0,0.001)I(-10,10)
    beta1.fV[g]~dnorm(0,0.001)I(-10,10)
    beta1.fF[g]~dnorm(0,0.001)I(-10,10)
    beta1.pV[g]~dnorm(0,0.001)I(-10,10)
    beta1.pF[g]~dnorm(0,0.001)I(-10,10)
    #calculate probabilities of four differnt groups to look at later... 
    logit(sV0[g])<- mu.sV[g]
    logit(sV1[g])<- mu.sV[g] + beta.sV[g]
    logit(sF0[g])<- mu.sF[g]
    logit(sF1[g])<- mu.sF[g] + beta.sF[g]
    logit(fV0[g])<- mu.fV[g]
    logit(fV1[g])<- mu.fV[g] + beta.fV[g]
    logit(fF0[g])<- mu.fF[g]
    logit(fF1[g])<- mu.fF[g] + beta.fF[g]
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
#Function to get starting values
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

zst=ms.init.z(rCH,f)
#
# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = known.state.ms(rCH),group=group, x=time_log, x1=logged_yrs2)

# Initial values
inits<-function(){list(mean.sV=c(.5,.5), mean.sF=c(.5,.5),mean.fV=c(.5,.5),mean.fF=c(.5,.5),mean.pV=c(.5,.5),mean.pF=c(.5,.5),sigma.sV=c(1,1),sigma.sF=c(1,1),sigma.fV=c(1,1),sigma.fF=c(1,1),sigma.pV=c(1,1),sigma.pF=c(1,1),z = zst)}
# Parameters monitored
parameters <- c("mean.sV","mean.sF", "mean.fV","mean.fF","mean.pV","mean.pF","beta.sV","beta.sF", "beta.fV","beta.fF","beta.pV","beta.pF","sV0","sV1","sF0","sF1","fV0","fV1","fF0","fF1","pV0","pV1","pF0","pF1","beta1.sV","beta1.sF","beta1.fV","beta1.fF","beta1.pV","beta1.pF","mu.sV","mu.sF","mu.fV","mu.fF","mu.pV","mu.pF","sigma.sV2","sigma.sF2","sigma.fV2","sigma.fF2","sigma.pV2","sigma.pF2")

# MCMC settings
ni <- 2500
nt <- 5
nb <- 500
nc <- 3

# Call JAGS from R #
#complex model
ms2.rf <- jags(jags.data, inits, parameters, "ms-ranef2stages.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T)
print(ms2.rf, digits=3)
###save all samples
mod.samples_2stage<- as.data.frame(do.call("rbind", ms2.rf$samples))
write.csv(mod.samples_2stage,"msmod_samples_2stagesim.csv",row.names=FALSE)

##Table of model summary, with betas
mod.sum<-ms2.rf$summary
write.csv(mod.sum,"isotria2stagemodsumsim.csv", row.names=T)
mod.table<-round(subset(mod.sum, select=c("mean","sd","Rhat", "n.eff")), digits=3)
mod.table2<-rbind(mod.table[which(substr(rownames(mod.table),1,2)=="mu"),],mod.table[which(substr(rownames(mod.table),1,4)=="beta"),])
write.csv(mod.table2,"isotria2stage_modtablesim.csv")
##Figures to compare estimates 2-stage model, 4-stage model, and true vital rates:

#mod.samples_2stage<-read.csv("msmod_samples_2stagesim.csv", header=TRUE)#if not loaded in already
head(mod.samples_2stage)
#mod.samples_4stage<-read.csv("msmod_samples_4stagesim.csv", header=TRUE)#if not loaded in already
head(mod.samples_4stage)
colnames(mod.samples_4stage)
colnames(mod.samples_2stage)
#ok, now try makeing a figure to compare:
#First, survival
survs<-c(sV,sV,sF,sF)
quartz(height=6, width=10)
par(mfcol=c(2,4))
for (i in 1:4){
hist(mod.samples_2stage[,i], main=paste("2stage",substr(colnames(mod.samples_2stage)[i],6,9)), col="gray")
abline(v=mean(mod.samples_2stage[,i]), lty=2, lwd=2)
abline(v=survs[i], col="red", lwd=2)

hist(mod.samples_4stage[,i], main=paste("4stage",substr(colnames(mod.samples_4stage)[i],6,9)), col="gray")
abline(v=mean(mod.samples_4stage[,i]), lty=2, lwd=2)
abline(v=survs[i], col="red", lwd=2)
}
#Now, Flowering prob:
quartz(height=6, width=10)
par(mfcol=c(2,4))
flows<-c(fV,fV,fF,fF)
for (i in 1:4){
  hist(mod.samples_2stage[,i+4], main=paste("2stage",substr(colnames(mod.samples_2stage)[i+4],6,9)), col="gray")
  abline(v=mean(mod.samples_2stage[,i+4]), lty=2, lwd=2)
  abline(v=flows[i], col="red", lwd=2)
  
  hist(mod.samples_4stage[,i+8], main=paste("4stage",substr(colnames(mod.samples_4stage)[i+8],6,9)), col="gray")
  abline(v=mean(mod.samples_4stage[,i+8]), lty=2, lwd=2)
  abline(v=flows[i], col="red", lwd=2)
}
#now dormancuy (=1-p for 2-stage model, d for 4-stage)
quartz(height=6, width=10)
par(mfcol=c(2,4))
dorms<-c(dV,dV,dF,dF)
for (i in 1:4){
  hist(1-mod.samples_2stage[,i+8], main=paste("2stage, 1-",substr(colnames(mod.samples_2stage)[i+8],6,9)), col="gray")
  abline(v=mean(1-mod.samples_2stage[,i+8]), lty=2, lwd=2)
  abline(v=dorms[i], col="red", lwd=2)
  
  hist(mod.samples_4stage[,i+16], main=paste("4stage",substr(colnames(mod.samples_4stage)[i+16],6,9)), col="gray")
  abline(v=mean(mod.samples_4stage[,i+16]), lty=2, lwd=2)
  abline(v=dorms[i], col="red", lwd=2)
}
#no check the dormant stages (just for 4-stage model)
quartz(height=6, width=7)
par(mfcol=c(2,2))
dsurvs<-c(sUV,sUV,sUF,sUF)

for (i in 1:4){
  hist(mod.samples_4stage[,i+4], main=paste("4stage",substr(colnames(mod.samples_4stage)[i+4],6,9)), col="gray")
  abline(v=mean(mod.samples_4stage[,i+4]), lty=2, lwd=2)
  abline(v=dsurvs[i], col="red", lwd=2)
}
#Now, Flowering prob:
quartz(height=6, width=7)
par(mfcol=c(2,4))
dflows<-c(fUV,fUV,fUF,fUF)
for (i in 1:2){
  hist(mod.samples_4stage[,i+12], main=paste("4stage",substr(colnames(mod.samples_4stage)[i+12],6,9)), col="gray")
  abline(v=mean(mod.samples_4stage[,i+12]), lty=2, lwd=2)
  abline(v=dflows[i], col="red", lwd=2)
}
#now dormancuy (=1-p for 2-stage model, d for 4-stage)
quartz(height=6, width=7)
par(mfcol=c(2,2))
ddorms<-c(dUV,dUV,dUF,dUF)
for (i in 1:4){
  hist(mod.samples_4stage[,i+20], main=paste("4stage",substr(colnames(mod.samples_4stage)[i+20],6,9)), col="gray")
  abline(v=mean(mod.samples_4stage[,i+20]), lty=2, lwd=2)
  abline(v=ddorms[i], col="red", lwd=2)
}