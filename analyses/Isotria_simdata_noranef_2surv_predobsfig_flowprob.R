#Isotria multistate model
#Simulate data with known parameters to test the ability of my model to recover them
#Do this for a range of values for each vital rate, so that we can assess the model's ability to recover true values
#Save all estimates, posterior, and model output for later use
#Modified May 16, 2017
#Modified again Sept 11, 2017 to run through different flowering probabilities, keeping survivial and dormancy constant

rm(list=ls()) 

options(stringsAsFactors=FALSE)
 
library(rjags)
library(jagsUI)
library(lattice)
library(coda)
library(boot)
#start with simple dataset with no variation between groups x and y, and with no annual variation in vital rates

#Range of values to use for each vital rate:
#after july discussion with andy, we do not want to use low estimates of survival- has trouble estimating with low values
#instead, use literature to get reasonable range of estimates for annual survival. 
#0.76-1.0 for cleistes bifaria (Gregg and Kery 2005, Kery et al 2005)
#0.96, 0.98 in orchis purpurea (JAcquemyn et al 2010)
#0.90-0.95 for Platanthera macrophylla & P. orbiculata (Cleavitt & Berry 2016)
#hutchings 2010 ophrys sphegodes has very low survival(0.2, half life=2.25years, but some live>20 years)

#since there are many, many potential combinations of vital rates to test, i will start with varying s, f, and d at 4 different values each. i will assume that flowering plants and vegetative plants have different, but correlated vital rates (i.e. when a given vital rate is high for flowering plants, it is also high  for vetetative plants)
#for now, assume that uf and uv plants have same d and f probs
#svvals<-c(0.7, 0.8, 0.9)
fvvals<-c(0.25,0.9)
#dvvals<-c(0.2, 0.4, 0.6, 0.8)
#sfvals<-c(0.7, 0.8, 0.9)
ffvals<-c(0.25,0.9)
#dfvals<-c(0.1, 0.3, 0.5, 0.7)
fufvals<-c(0.15,0.8)
#dufvals<-c(0.5, 0.6, 0.7, 0.8)
fuvvals<-c(0.15,0.8)
#duvvals<-c(0.5, 0.6, 0.7, 0.8)
#first, define mean survival, flowering probability, and dormancy/emergence probability, 
#as well as the number of occassions, states, observations, and marked individuals
#for now, andy said to vary survival only- use realistic values for f and d

#Create dataframes with true values and estimated values, to be used to make a figure of estimated vs true values for paper
sV.df<-data.frame(sV_true=numeric(),
                   sV_est=numeric(),
                   sV_est_2.5=numeric(),
                   sV_est_97.5=numeric(),
                   stringsAsFactors = FALSE)
sF.df<-data.frame(sF_true=numeric(),
                  sF_est=numeric(),
                  sF_est_2.5=numeric(),
                  sF_est_97.5=numeric(),
                  stringsAsFactors = FALSE)
fV.df<-data.frame(fV_true=numeric(),
                  fV_est=numeric(),
                  fV_est_2.5=numeric(),
                  fV_est_97.5=numeric(),
                  stringsAsFactors = FALSE)
fF.df<-data.frame(fF_true=numeric(),
                  fF_est=numeric(),
                  fF_est_2.5=numeric(),
                  fF_est_97.5=numeric(),
                  stringsAsFactors = FALSE)
fUV.df<-data.frame(fUV_true=numeric(),
                   fUV_est=numeric(),
                   fUV_est_2.5=numeric(),
                   fUV_est_97.5=numeric(),
                   stringsAsFactors = FALSE)
fUF.df<-data.frame(fUF_true=numeric(),
                   fUF_est=numeric(),
                   fUF_est_2.5=numeric(),
                   fUF_est_97.5=numeric(),
                   stringsAsFactors = FALSE)
dV.df<-data.frame(dV_true=numeric(),
                  dV_est=numeric(),
                  dV_est_2.5=numeric(),
                  dV_est_97.5=numeric(),
                  stringsAsFactors = FALSE)
dF.df<-data.frame(dF_true=numeric(),
                  dF_est=numeric(),
                  dF_est_05=numeric(),
                  dF_est_97.5=numeric(),
                  stringsAsFactors = FALSE)
dUV.df<-data.frame(dUV_true=numeric(),
                   dUV_est=numeric(),
                   dUV_est_2.5=numeric(),
                   dUV_est_97.5=numeric(),
                   stringsAsFactors = FALSE)
dUF.df<-data.frame(dUF_true=numeric(),
                   dUF_est=numeric(),
                   dUF_est_2.5=numeric(),
                   dUF_est_97.5=numeric(),
                   stringsAsFactors = FALSE)

#for(a in 1:length(svvals)){
for (b in 1:length(fvvals)){
#for(c in 1:length(dvvals)){
  for(j in 1:20){#do each simulation 20 times to check model. eventually will want to do this 100 times
      #sV<-svvals[a]#survival=phi
      #sF<-sfvals[a]
      #sUV<-svvals[a]
      #sUF<-sfvals[a]
      fV<-fvvals[b]#flowering probability
      fF<- ffvals[b]
      fUV<-fuvvals[b]
      fUF<-fufvals[b]
      #dV<-dvvals[c]#probability of dormancy
      #dF<-dfvals[c]
      #dUV<-duvvals[c]
      #dUF<-dufvals[c]    
      sV<-0.7#survival=phi
      sF<-0.7
      sUV<-0.7
      sUF<-0.7
      #fV<-0.3
      #fF<- 0.5
      #fUV<-0.2
      #fUF<-0.4
      dV<-0.5
      dF<-0.5
      dUV<-0.5
      dUF<-0.5
    
n.occasions<-31#31 year dataset
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
  #CH: capture-histories to be used
  #CH.TRUE: capture histories with perfect observation
}
sim<-simul.ms(STATE,OBS,marked)
CH<-sim$CH
#Replace first observations that are 0 with NAs
CH[which(CH[,1]==0),1]<-NA
#compute vector with occasion of first capture
get.first<-function(x) min(which(x!=0))
#CH2<-rbind(CH,CH)
CH2<-CH
f<-apply(CH2,1,get.first)
#recode CH matrix since 0 is not allowed in WinBUGS
rCH<-CH2
rCH[rCH==0]<-3
rCH[rCH==4]<-3
n=dim(rCH)[1]#number plants in sulated dataset
###now use my 4-stage model to estimate parameters of these data
#group<-c(rep(1,n/2))#1 group only for simulations

#Fit a 4-stage multistate model in JAGS, with no random effects of year, and shared survival for V-UV and F-UF plants####
sink("ms-4stages_2surv.jags")
cat("
    
    model {
    # -------------------------------------------------
    # Parameters:
    # sV: survival probability for nonreproductive plants (=sUV: survival probability for dormant plants that were vegetative above ground for this model)
    # sF: survival probability for reproductive plants(=sUF: survival probability for dormant plants that were reproductive above ground) for this model
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
    logit(sV[i,t]) <- eta.sV[t]
    logit(sF[i,t]) <- eta.sF[t]
    #logit(sUV[i,t]) <- eta.sUV[t]
    #logit(sUF[i,t]) <- eta.sUF[t]
    logit(fV[i,t]) <- eta.fV[t]
    logit(fF[i,t]) <- eta.fF[t]
    logit(fUV[i,t]) <- eta.fUV[t]
    logit(fUF[i,t]) <- eta.fUF[t]
    logit(dV[i,t]) <- eta.dV[t]
    logit(dF[i,t]) <- eta.dF[t]
    logit(dUV[i,t]) <- eta.dUV[t]
    logit(dUF[i,t]) <- eta.dUF[t]
    }#t
    }#i
    
    for (t in 1:(n.occasions-1)){
    eta.sV[t]  <-mu.sV
    eta.sF[t]  <-mu.sF
    #eta.sUV[t]  <-mu.sUV
    #eta.sUF[t]  <-mu.sUF
    eta.fV[t]  <-mu.fV
    eta.fF[t]  <-mu.fF
    eta.fUV[t]  <-mu.fUV
    eta.fUF[t]  <-mu.fUF
    eta.dV[t]  <-mu.dV
    eta.dF[t]  <-mu.dF
    eta.dUV[t]  <-mu.dUV
    eta.dUF[t]  <-mu.dUF
    
    }#t
    mean.sV~dunif(0,1)# Priors for mean survival 
    mu.sV<-log(mean.sV/(1-mean.sV)) 
    mean.sF~dunif(0,1)# Priors for mean survival 
    mu.sF<-log(mean.sF/(1-mean.sF))  
    #mean.sUV~dunif(0,1)# Priors for mean group-specific survival 
    #mu.sUV<-log(mean.sUV/(1-mean.sUV))  
    #mean.sUF~dunif(0,1)# Priors for mean group-specific survival 
    #mu.sUF<-log(mean.sUF/(1-mean.sUF))  
    mean.fV~dunif(0,1)# Priors for mean group-specific prob of rep for veg plants
    mu.fV<-log(mean.fV/(1-mean.fV))       
    mean.fF~dunif(0,1)# Priors for mean group-specific prob of rep for rep plants
    mu.fF<-log(mean.fF/(1-mean.fF))    
    mean.fUV~dunif(0,1)# Priors for mean group-specific prob of rep for dorm (prev veg) plants
    mu.fUV<-log(mean.fUV/(1-mean.fUV)) 
    mean.fUF~dunif(0,1)# Priors for mean group-specific prob of rep for dorm (prev rep) plants
    mu.fUF<-log(mean.fUF/(1-mean.fUF)) 
    mean.dV~dunif(0,1)# Priors for mean group-specific dormancy for veg plants
    mu.dV<-log(mean.dV/(1-mean.dV))       
    mean.dF~dunif(0,1)# Priors for mean group-specific dormancy for rep plants
    mu.dF<-log(mean.dF/(1-mean.dF))    
    mean.dUV~dunif(0,1)# Priors for mean group-specific dormancy for dorm (prev veg) plants
    mu.dUV<-log(mean.dUV/(1-mean.dUV))
    mean.dUF~dunif(0,1)# Priors for mean group-specific dormancy for dorm (prev rep) plants
    mu.dUF<-log(mean.dUF/(1-mean.dUV)) 

 
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
  for(i in 1:dim(ch)[1]){#warning for 163(ok), 400(ok), 404(ok), 424(ok), 456(ok),458 (ok), 479(ok), 509(ok), 513(ok), 528(ok), 529(ok), 551(ok), 553(ok), are others ok?
    z<-ch[i,]
    v <- which(z>=3)# occurences that are unknown states get 3s by default (already coded)
    w<-which(z<3)
    x<-which(z==2)
    z[-v] <- NA#observed states get NA
    if(length(x)==0){z[w]<-NA
    } else if
    (length(which(z==3))> 0 & length(x)>0 & !('2' %in% ch[i,which(z==3)-1])){z[w]<-NA
    } else if
    (length(which(ch[i,]>=3))==0){z[w]<-NA
    } else if (length(which(z==3))> 0 & length(x)>0 & '2' %in% ch[i,which(z==3)-1]){
      if(length(x)>=1){
        temp<-ch[i,]
        temp[which(is.na(temp))]<-0
        for (j in 1:sum(rle(temp)$values == 2)){
          if(is.na(rle(temp)$values[(which(rle(temp)$values==2)+1)[j]]==3)){next}
          if(rle(temp)$values[(which(rle(temp)$values==2)+1)[j]]==3){
         fourstart<-sum(rle(temp)$lengths[1:which(rle(temp)$values==2)[j]])+1
          threestofours<-which(rle(temp)$values==3)
          threestofours<-threestofours[threestofours>which(rle(temp)$values==2)[j]]
          print(threestofours)
          fourend<-sum(rle(temp)$lengths[1:threestofours])
          z[fourstart:fourend]<-4
          }
          else {next}
          }#j
      }
    }
    zstart<-rbind(zstart,z)
  }
      return(zstart)
}

#for (i in 1:length(sim_ds){
  
zst=ms.init.z(rCH,f)
# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = known.state.ms(rCH))
 

 

# Initial values
inits<-function(){list(z = zst)}
# Parameters monitored
parameters <- c("mean.sV","mean.sF","mean.fV","mean.fF","mean.fUV","mean.fUF","mean.dV","mean.dF","mean.dUV","mean.dUF")

# MCMC settings
ni <- 10000
nt <- 10
nb <- 2000
#na<-2000
nc <- 3
# Call JAGS from R #
#complex model
ms4.rf <- jags(jags.data, inits, parameters, "ms-4stages_2surv.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin=nb, parallel=T)#n.burnin = nb,

print(ms4.rf, digits=3)
#plot(ms4.rf)

###save all samples
name<-paste("fVfF", fV,"fUVfUF",fUV,j,".csv", sep="_")
path="C:/Users/aettinger/Documents/isotria/simdat_modsamples/fprobs"
file<-file.path(path,name)
mod.samples_4stage<- as.data.frame(do.call("rbind", ms4.rf$samples))
write.csv(mod.samples_4stage,paste(file),row.names=T)

##Table of model summary, with betas
mod.sum<-ms4.rf$summary
namesum<-paste("modsum",name,sep="_")
pathsum="C:/Users/aettinger/Documents/isotria/simdat_modsums/fprobs"
filesum<-file.path(pathsum,namesum)
write.csv(mod.sum,paste(filesum),row.names=T)
#ok, now try making a figure to compare:
#First, survival
survs<-c(sV,sF)
x11(height=6, width=12)
par(mfcol=c(2,5))
for (i in 1:2){
  hist(mod.samples_4stage[,i], main=paste(substr(colnames(mod.samples_4stage)[i],6,7)), col="gray")
  abline(v=mean(mod.samples_4stage[,i]), lty=2, lwd=2)
  abline(v=survs[i], col="red", lwd=1)
  
 if(i==1) {mtext(paste("sv=", sV,"sf=",sF,":",j),side=3, line=3, adj=8)}
}

#Now, Flowering prob:
flows<-c(fV,fF,fUV,fUF)
for (i in 1:4){
  hist(mod.samples_4stage[,i+2], main=paste(substr(colnames(mod.samples_4stage)[i+2],6,8)), col="gray")
  abline(v=mean(mod.samples_4stage[,i+2]), lty=2, lwd=2)
  abline(v=flows[i], col="red", lwd=2)
}
#now dormancy (=1-p for 2-stage model, d for 4-stage)

dorms<-c(dV,dF,dUV,dUF)
for (i in 1:4){
  hist(mod.samples_4stage[,i+6], main=paste(substr(colnames(mod.samples_4stage)[i+6],6,8)), col="gray")
  abline(v=mean(mod.samples_4stage[,i+6]), lty=2, lwd=2)
  abline(v=dorms[i], col="red", lwd=2)
  
}
#fill in dataframes for this iteration
#for now, to create and test the code: 
#ms4.rf<-read.csv("simdat_modsums/modsum_sv_0.3_sf_0.5_1_.csv", header=T)
#row.names(ms4.rf)<-ms4.rf[,1]
#filesum<-file.path(pathsum,namesum)
#write.csv(mod.sum,paste(filesum),row.names=T)
#ms4.rf<-ms4.rf[,-1]
sV.df[j,1]<-sV
sV.df[j,2]<-ms4.rf$summary[grep("mean.sV",rownames(ms4.rf$summary)),grep("mean",colnames(ms4.rf$summary))]
#sV.df[j,2]<-ms4.rf$mean[1]
sV.df[j,3]<-ms4.rf$summary[grep("mean.sV",rownames(ms4.rf$summary)),grep("2.5",colnames(ms4.rf$summary))]#sV_est_2.5=credible interval
#sV.df[j,3]<-ms4.rf[1,grep("2.5",colnames(ms4.rf))]
sV.df[j,4]<-ms4.rf$summary[grep("mean.sV",rownames(ms4.rf$summary)),grep("97.5",colnames(ms4.rf$summary))]#sV_est_97.5=credible interval                  
#sV.df[j,4]<-ms4.rf[1,grep("97.5",colnames(ms4.rf))]

sF.df[j,1]<-sF
sF.df[j,2]<-ms4.rf$summary[grep("mean.sF",rownames(ms4.rf$summary)),grep("mean",colnames(ms4.rf$summary))]#sF_est
#sF.df[j,2]<-ms4.rf$mean[2]
sF.df[j,3]<-ms4.rf$summary[grep("mean.sF",rownames(ms4.rf$summary)),grep("2.5",colnames(ms4.rf$summary))]#sF_est_2.5=credible interval
#sF.df[j,3]<-ms4.rf[2,grep("2.5",colnames(ms4.rf))]
sF.df[j,4]<-ms4.rf$summary[grep("mean.sF",rownames(ms4.rf$summary)),grep("97.5",colnames(ms4.rf$summary))]#sF_est_97.5=credible interval                  
#sF.df[j,4]<-ms4.rf[2,grep("97.5",colnames(ms4.rf))]

fV.df[j,1]<-fvvals[b]
fV.df[j,2]<-ms4.rf$summary[grep("mean.fV",rownames(ms4.rf$summary)),grep("mean",colnames(ms4.rf$summary))]
#fV.df[j,2]<-ms4.rf$mean[3]
fV.df[j,3]<-ms4.rf$summary[grep("mean.fV",rownames(ms4.rf$summary)),grep("2.5",colnames(ms4.rf$summary))]#fV_est_2.5=credible interval
#fV.df[j,3]<-ms4.rf[3,grep("2.5",colnames(ms4.rf))]
fV.df[j,4]<-ms4.rf$summary[grep("mean.fV",rownames(ms4.rf$summary)),grep("97.5",colnames(ms4.rf$summary))]#fV_est_97.5=credible interval  
#fV.df[j,4]<-ms4.rf[3,grep("97.5",colnames(ms4.rf))]

fF.df[j,1]<-ffvals[b]
fF.df[j,2]<-ms4.rf$summary[grep("mean.fF",rownames(ms4.rf$summary)),grep("mean",colnames(ms4.rf$summary))]#fF_est
#fF.df[j,2]<-ms4.rf$mean[4]
fF.df[j,3]<-ms4.rf$summary[grep("mean.fF",rownames(ms4.rf$summary)),grep("2.5",colnames(ms4.rf$summary))]#fF_est_2.5=credible interval
#fF.df[j,3]<-ms4.rf[4,grep("2.5",colnames(ms4.rf))]
fF.df[j,4]<-ms4.rf$summary[grep("mean.fF",rownames(ms4.rf$summary)),grep("97.5",colnames(ms4.rf$summary))]#fF_est_97.5=credible interval     
#fF.df[j,4]<-ms4.rf[4,grep("97.5",colnames(ms4.rf))]

fUV.df[j,1]<-fuvvals[b]
fUV.df[j,2]<-ms4.rf$summary[grep("mean.fUV",rownames(ms4.rf$summary)),grep("mean",colnames(ms4.rf$summary))]
#fUV.df[j,2]<-ms4.rf$mean[5]
fUV.df[j,3]<-ms4.rf$summary[grep("mean.fUV",rownames(ms4.rf$summary)),grep("2.5",colnames(ms4.rf$summary))]#fUV_est_2.5=credible interval
#fUV.df[j,3]<-ms4.rf[5,grep("2.5",colnames(ms4.rf))]
fUV.df[j,4]<-ms4.rf$summary[grep("mean.fUV",rownames(ms4.rf$summary)),grep("97.5",colnames(ms4.rf$summary))]#fUV_est_97.5=credible interval                  
#fUV.df[j,4]<-ms4.rf[5,grep("97.5",colnames(ms4.rf))]

fUF.df[j,1]<-fufvals[b]

fUF.df[j,2]<-ms4.rf$summary[grep("mean.fUF",rownames(ms4.rf$summary)),grep("mean",colnames(ms4.rf$summary))]#fUF_est
#fUF.df[j,2]<-ms4.rf$mean[6]
fUF.df[j,3]<-ms4.rf$summary[grep("mean.fUF",rownames(ms4.rf$summary)),grep("2.5",colnames(ms4.rf$summary))]#fUF_est_2.5=credible interval
#fUF.df[j,3]<-ms4.rf[6,grep("2.5",colnames(ms4.rf))]

fUF.df[j,4]<-ms4.rf$summary[grep("mean.fUF",rownames(ms4.rf$summary)),grep("97.5",colnames(ms4.rf$summary))]#fUF_est_97.5=credible interval  
#fUF.df[j,4]<-ms4.rf[6,grep("97.5",colnames(ms4.rf))]
dV.df[j,1]<-dV
dV.df[j,2]<-ms4.rf$summary[grep("mean.dV",rownames(ms4.rf$summary)),grep("mean",colnames(ms4.rf$summary))]
#dV.df[j,2]<-ms4.rf$mean[7]
dV.df[j,3]<-ms4.rf$summary[grep("mean.dV",rownames(ms4.rf$summary)),grep("2.5",colnames(ms4.rf$summary))]#dV_est_2.5=credible interval
#dV.df[j,3]<-ms4.rf[7,grep("2.5",colnames(ms4.rf))]
dV.df[j,4]<-ms4.rf$summary[grep("mean.dV",rownames(ms4.rf$summary)),grep("97.5",colnames(ms4.rf$summary))]#dV_est_97.5=credible interval                  
#dV.df[j,4]<-ms4.rf[7,grep("97.5",colnames(ms4.rf))]
dF.df[j,1]<-dF
dF.df[j,2]<-ms4.rf$summary[grep("mean.dF",rownames(ms4.rf$summary)),grep("mean",colnames(ms4.rf$summary))]#dF_est
#dF.df[j,2]<-ms4.rf$mean[8]
dF.df[j,3]<-ms4.rf$summary[grep("mean.dF",rownames(ms4.rf$summary)),grep("2.5",colnames(ms4.rf$summary))]#dF_est_2.5=credible interval
#dF.df[j,3]<-ms4.rf[8,grep("2.5",colnames(ms4.rf))]
dF.df[j,4]<-ms4.rf$summary[grep("mean.dF",rownames(ms4.rf$summary)),grep("97.5",colnames(ms4.rf$summary))]#dF_est_97.5=credible interval     
#dF.df[j,4]<-ms4.rf[8,grep("97.5",colnames(ms4.rf))]
dUV.df[j,1]<-dUV
dUV.df[j,2]<-ms4.rf$summary[grep("mean.dUV",rownames(ms4.rf$summary)),grep("mean",colnames(ms4.rf$summary))]
#dUV.df[j,2]<-ms4.rf$mean[9]
dUV.df[j,3]<-ms4.rf$summary[grep("mean.dUV",rownames(ms4.rf$summary)),grep("2.5",colnames(ms4.rf$summary))]#dUV_est_2.5=credible interval
#dUV.df[j,3]<-ms4.rf[9,grep("2.5",colnames(ms4.rf))]
dUV.df[j,4]<-ms4.rf$summary[grep("mean.dUV",rownames(ms4.rf$summary)),grep("97.5",colnames(ms4.rf$summary))]#dUV_est_97.5=credible interval                  
#dUV.df[j,4]<-ms4.rf[9,grep("97.5",colnames(ms4.rf))]
dUF.df[j,1]<-dUF
dUF.df[j,2]<-ms4.rf$summary[grep("mean.dUF",rownames(ms4.rf$summary)),grep("mean",colnames(ms4.rf$summary))]#dUF_est
#dUF.df[j,2]<-ms4.rf$mean[10]
dUF.df[j,3]<-ms4.rf$summary[grep("mean.dUF",rownames(ms4.rf$summary)),grep("2.5",colnames(ms4.rf$summary))]#dUF_est_2.5=credible interval
#dUF.df[j,3]<-ms4.rf[10,grep("2.5",colnames(ms4.rf))]
dUF.df[j,4]<-ms4.rf$summary[grep("mean.dUF",rownames(ms4.rf$summary)),grep("97.5",colnames(ms4.rf$summary))]#dUF_est_97.5=credible interval  
#dUF.df[j,4]<-ms4.rf[10,grep("97.5",colnames(ms4.rf))]
}#j
#    }#c
  }#b
# }#a


###If code stops running before full set of simulations has been done, then use the following:


pathsum="C:/Users/aettinger/Documents/isotria/simdat_modsums/fprobs"
files<-list.files(path = pathsum)
for(j in 1:length(files)){
namesum<-files[j]
filesum<-file.path(pathsum,namesum)
ms4.rf<-read.csv(filesum)
row.names(ms4.rf)<-ms4.rf[,1]
ms4.rf<-ms4.rf[,-1]
#b<-which(fvvals==as.numeric(substr(namesum,13,15)))
sV.df[j,1]<-sV
sV.df[j,2]<-ms4.rf$mean[1]
sV.df[j,3]<-ms4.rf[1,grep("2.5",colnames(ms4.rf))]
sV.df[j,4]<-ms4.rf[1,grep("97.5",colnames(ms4.rf))]

sF.df[j,1]<-sF
sF.df[j,2]<-ms4.rf$mean[2]
sF.df[j,3]<-ms4.rf[2,grep("2.5",colnames(ms4.rf))]
sF.df[j,4]<-ms4.rf[2,grep("97.5",colnames(ms4.rf))]
fV.df[j,1]<-as.numeric(substr(namesum,13,15))
print(fV.df[j,1]);
fV.df[j,2]<-ms4.rf$mean[3]
fV.df[j,3]<-ms4.rf[3,grep("2.5",colnames(ms4.rf))]
fV.df[j,4]<-ms4.rf[3,grep("97.5",colnames(ms4.rf))]
fF.df[j,1]<-as.numeric(substr(namesum,13,15))
print(fF.df[j,1])
fF.df[j,2]<-ms4.rf$mean[4]
fF.df[j,3]<-ms4.rf[4,grep("2.5",colnames(ms4.rf))]
fF.df[j,4]<-ms4.rf[4,grep("97.5",colnames(ms4.rf))]
fUV.df[j,1]<-as.numeric(substr(namesum,24,26))
fUV.df[j,2]<-ms4.rf$mean[5]
print(fUV.df[j,1])
fUV.df[j,3]<-ms4.rf[5,grep("2.5",colnames(ms4.rf))]
fUV.df[j,4]<-ms4.rf[5,grep("97.5",colnames(ms4.rf))]
fUF.df[j,1]<-as.numeric(substr(namesum,24,26))
print(fUF.df[j,1])
fUF.df[j,2]<-ms4.rf$mean[6]
fUF.df[j,3]<-ms4.rf[6,grep("2.5",colnames(ms4.rf))]
fUF.df[j,4]<-ms4.rf[6,grep("97.5",colnames(ms4.rf))]

dV.df[j,1]<-dV
dV.df[j,2]<-ms4.rf$mean[7]
dV.df[j,3]<-ms4.rf[7,grep("2.5",colnames(ms4.rf))]
dV.df[j,4]<-ms4.rf[7,grep("97.5",colnames(ms4.rf))]
dF.df[j,1]<-dF
dF.df[j,2]<-ms4.rf$mean[8]
dF.df[j,3]<-ms4.rf[8,grep("2.5",colnames(ms4.rf))]
dF.df[j,4]<-ms4.rf[8,grep("97.5",colnames(ms4.rf))]
dUV.df[j,1]<-dUV
dUV.df[j,2]<-ms4.rf$mean[9]
dUV.df[j,3]<-ms4.rf[9,grep("2.5",colnames(ms4.rf))]
dUV.df[j,4]<-ms4.rf[9,grep("97.5",colnames(ms4.rf))]
dUF.df[j,1]<-dUF
dUF.df[j,2]<-ms4.rf$mean[10]
dUF.df[j,3]<-ms4.rf[10,grep("2.5",colnames(ms4.rf))]
dUF.df[j,4]<-ms4.rf[10,grep("97.5",colnames(ms4.rf))]
}#j

#Now make a figure of estimate vs true SURVIVAL:
sV.df[,1]<-as.numeric(sV.df[,1])
sV.df[,2]<-as.numeric(sV.df[,2])
x11(height=6, width=10)
par(mfrow=c(1,2))
plot(sV.df[,1],sV.df[,2], pch=16, col="gray", xlim=c(0,1), ylim=c(0,1), xlab="True value", ylab="Estimate", main="sV (2 Surv. Mod)")
abline(a=0, b=1, lty=1 )
for(i in 1:dim(sV.df)[1]){
  arrows(sV.df[i,1],sV.df[i,3],sV.df[i,1],sV.df[i,4], length=0.01, angle=90, code=3)
}
points(sV.df[,1],sV.df[,2], pch=21, bg="gray")

sF.df[,1]<-as.numeric(sF.df[,1])
sF.df[,2]<-as.numeric(sF.df[,2])
plot(sF.df[,1],sF.df[,2], pch=16, col="gray", xlim=c(0,1), ylim=c(0,1), xlab="True value", ylab="Estimate", main="sF (2 Surv. Mod)")
abline(a=0, b=1, lty=1 )
for(i in 1:dim(sF.df)[1]){
  arrows(sF.df[i,1],sF.df[i,3],sF.df[i,1],sF.df[i,4], length=0.01, angle=90, code=3)
}
points(sF.df[,1],sF.df[,2], pch=21, bg="gray")

#Now make a figure of estimate vs true FLOWERING:
fV.df[,1]<-as.numeric(fV.df[,1])
fV.df[,2]<-as.numeric(fV.df[,2])
x11(height=6, width=15)
par(mfrow=c(1,4))
plot(fV.df[,1],fV.df[,2], pch=16, col="gray", xlim=c(0,1), ylim=c(0,1), xlab="True value", ylab="Estimate", main="fV (2 Surv. Mod)")
abline(a=0, b=1, lty=1 )
for(i in 1:dim(fV.df)[1]){
  arrows(fV.df[i,1],fV.df[i,3],fV.df[i,1],fV.df[i,4], length=0.01, angle=90, code=3)
}
points(fV.df[,1],fV.df[,2], pch=21, bg="gray")

fF.df[,1]<-as.numeric(fF.df[,1])
fF.df[,2]<-as.numeric(fF.df[,2])
plot(fF.df[,1],fF.df[,2], pch=16, col="gray", xlim=c(0,1), ylim=c(0,1), xlab="True value", ylab="Estimate", main="fF (2 Surv. Mod)")
abline(a=0, b=1, lty=1 )
for(i in 1:dim(fF.df)[1]){
  arrows(fF.df[i,1],fF.df[i,3],fF.df[i,1],fF.df[i,4], length=0.01, angle=90, code=3)
}
points(fF.df[,1],fF.df[,2], pch=21, bg="gray")
#fUV
plot(fUV.df[,1],fUV.df[,2], pch=16, col="gray", xlim=c(0,1), ylim=c(0,1), xlab="True value", ylab="Estimate", main="fUV (2 Surv. Mod)")
abline(a=0, b=1, lty=1 )
for(i in 1:dim(fUV.df)[1]){
  arrows(fUV.df[i,1],fUV.df[i,3],fUV.df[i,1],fUV.df[i,4], length=0.01, angle=90, code=3)
}
points(fUV.df[,1],fUV.df[,2], pch=21, bg="gray")
#fUF

fUF.df[,1]<-as.numeric(fUF.df[,1])
fUF.df[,2]<-as.numeric(fUF.df[,2])
plot(fUF.df[,1],fUF.df[,2], pch=16, col="gray", xlim=c(0,1), ylim=c(0,1), xlab="True value", ylab="Estimate", main="fUF (2 Surv. Mod)")
abline(a=0, b=1, lty=1 )
for(i in 1:dim(fUF.df)[1]){
  arrows(fUF.df[i,1],fUF.df[i,3],fUF.df[i,1],fUF.df[i,4], length=0.01, angle=90, code=3)
}
points(fUF.df[,1],fUF.df[,2], pch=21, bg="gray")

#Now make a figure of estimate vs true dormancy:
dV.df[,1]<-as.numeric(dV.df[,1])
dV.df[,2]<-as.numeric(dV.df[,2])
x11(height=6, width=10)
par(mfrow=c(1,2))
plot(dV.df[,1],dV.df[,2], pch=21, bg="gray", xlim=c(0,1), ylim=c(0,1), xlab="True value", ylab="Estimate", main="dV (2 Surv. Mod)")
abline(a=0, b=1, lty=1 )
for(i in 1:dim(dV.df)[1]){
  arrows(dV.df[i,1],dV.df[i,3],dV.df[i,1],dV.df[i,4], length=0.01, angle=90, code=3)
}
points(dV.df[,1],dV.df[,2], pch=21, bg="gray")


dF.df[,1]<-as.numeric(dF.df[,1])
dF.df[,2]<-as.numeric(dF.df[,2])

plot(dF.df[,1],dF.df[,2], pch=21, bg="gray", xlim=c(0,1), ylim=c(0,1), xlab="True value", ylab="Estimate", main="dF (2 Surv. Mod)")
abline(a=0, b=1, lty=1 )
for(i in 1:dim(dF.df)[1]){
  arrows(dF.df[i,1],dF.df[i,3],dF.df[i,1],dF.df[i,4], length=0.01, angle=90, code=3)
}
points(dF.df[,1],dF.df[,2], pch=21, bg="gray")
