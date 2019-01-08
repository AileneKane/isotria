#Isotria multistate model
#Simulate data with known parameters to test the ability of the two-stage model to recover them
#Do this for a range of values for each vital rate, so that we can assess the model's ability to recover true values
#Save all estimates, posterior, and model output for later use
#Modified from 4-stage code June 14, 2018 
#Modified July 2018 to run through different survival probabilities, keeping survivial and dormancy constant

rm(list=ls()) 
setwd(setwd("~/Documents/GitHub/isotria"))
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
#svvals<-c(0.6,0.7, 0.8, 0.9)
#fvvals<-c(0.1,0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9)
#fvvals<-c(0.5, 0.6, 0.7, 0.8, 0.9)
#fuvvals<-c(0.5, 0.6, 0.7, 0.8, 0.9)
dvvals<-c(0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5, 0.55,0.6, 0.65, 0.7,0.75, 0.8, 0.85, 0.9,0.95)
#sfvals<-c(0.6, 0.7, 0.8, 0.9)
#ffvals<-c(0.5, 0.6, 0.7, 0.8, 0.9)
#fufvals<-c(0.5, 0.6, 0.7, 0.8, 0.9)

dfvals<-c(0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5, 0.55,0.6, 0.65, 0.7,0.75, 0.8, 0.85, 0.9,0.95)
#first, define mean survival, flowering probability, and dormancy/emergence probability, 
#as well as the number of occassions, states, observations, and marked individuals
duvvals<-c(0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5, 0.55,0.6, 0.65, 0.7,0.75, 0.8, 0.85, 0.9,0.95)
dufvals<-c(0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5, 0.55,0.6, 0.65, 0.7,0.75, 0.8, 0.85, 0.9,0.95)

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

#for(a in 1:length(svvals)){
#for (b in 1:length(fvvals)){
  for(c in 1:length(dvvals)){
  for(j in 1:3){#do each simulation 20 times to check model. eventually will want to do this 100 times
    #sV<-svvals[a]#survival=phi
    #sF<-sfvals[a]
    
    #fV<-fvvals[b]#flowering probability
    #fF<- ffvals[b]
    #fUV<-fuvvals[b]
   # fUF<-fufvals[b]
    dV<-dvvals[c]#probability of dormancy
    dF<-dfvals[c]
    dUV<-duvvals[c]
    dUF<-dufvals[c]    
    #sV<-0.7#survival=phi
    #sF<-0.7
   sV<-0.7
    sF<-0.7
    fV<-0.5
    fF<- 0.7
    fUV<-0.5
    fUF<-0.7
    #dV<-0.5
    #dF<-0.5
    #dUV<-0.5
    #dUF<-0.5
    
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
          sV*(1-fUV)*(1-dUV),sV*fUV*(1-dUV),sV*dUV,0,1-sV,
          sF*(1-fUF)*(1-dUF),sF*fUF*(1-dUF),0,sF*dUF,1-sF,
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
    n=dim(rCH)[1]#number plants in simulated dataset
    ###now use my 2-stage model to estimate parameters of these data

#Fit a 2-stage multistate model in JAGS, with no random effects of year####
sink("ms-2stages.jags")
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
    # dV: dormancy probability for nonreproductive plants
    # dF: dormancy  probability for reproductive plants
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
    
    for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
    logit(sV[i,t]) <- eta.sV[t]
    logit(sF[i,t]) <- eta.sF[t]
    logit(fV[i,t]) <- eta.fV[t]
    logit(fF[i,t]) <- eta.fF[t]
    logit(dV[i,t]) <- eta.dV[t]
    logit(dF[i,t]) <- eta.dF[t]
  }#t
}#i

    for (t in 1:(n.occasions-1)){
    eta.sV[t]  <-mu.sV
    eta.sF[t]  <-mu.sF
    eta.fV[t]  <-mu.fV
    eta.fF[t]  <-mu.fF
    eta.dV[t]  <-mu.dV
    eta.dF[t]  <-mu.dF
  }#t
    mean.sV~dunif(0,1)# Priors for mean survival 
    mu.sV<-log(mean.sV/(1-mean.sV)) 
    mean.sF~dunif(0,1)# Priors for mean survival 
    mu.sF<-log(mean.sF/(1-mean.sF))  
     mean.fV~dunif(0,1)# Priors for mean group-specific prob of rep for veg plants
    mu.fV<-log(mean.fV/(1-mean.fV))       
    mean.fF~dunif(0,1)# Priors for mean group-specific prob of rep for rep plants
    mu.fF<-log(mean.fF/(1-mean.fF))    
    mean.dV~dunif(0,1)# Priors for mean group-specific dormancy for veg plants
    mu.dV<-log(mean.dV/(1-mean.dV))       
    mean.dF~dunif(0,1)# Priors for mean group-specific dormancy for rep plants
    mu.dF<-log(mean.dF/(1-mean.dF))    
    
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
    po[1,i,t,1] <- 1-dV[i,t]
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- dV[i,t]
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- 1-dF[i,t]
    po[2,i,t,3] <- dF[i,t]
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

zst=ms.init.z(rCH,f)

# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = known.state.ms(rCH))

# Initial values
#inits3a<-function(){list(mean.sV=c(.5,.5), mean.sF=c(.5,.5),mean.fV=c(.5,.5),mean.fF=c(.5,.5),mean.dV=c(.5,.5),mean.dF=c(.5,.5),sigma.sF=c(1,1),sigma.fV=c(1,1),sigma.fF=c(1,1),sigma.pV=c(1,1),sigma.pF=c(1,1),z = zst1)}#Error in jags.model(file = model.file, data = data, inits = inits, n.chains = n.chains,  : #Error in node y[1,3]
inits<-function(){list(z = zst)}

# Parameters monitored
parameters <- c("mean.sV","mean.sF","mean.fV","mean.fF","mean.dV","mean.dF")

# MCMC settings
ni <- 10000
nt <- 10
nb <- 2000
nc <- 2

# Call JAGS from R #
#complex model
ms2stage <- jags(jags.data, inits, parameters, "ms-2stages.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T)
print(ms2stage, digits=3)


##Table of model summary, with betas
mod.sum<-ms2stage$summary
name<-paste("dV", dV,"dF",dF,j,".csv", sep="_")

namesum<-paste("modsum",name,sep="_")
pathsum="/Users/aileneettinger/Documents/GitHub/isotria/analyses/simdat_modsums/2stage/dprobs"
filesum<-file.path(pathsum,namesum)
write.csv(mod.sum,paste(filesum),row.names=T)
}#j
    }#c
#  }#b
 #}#a


###If code stops running before full set of simulations has been done, then use the following:


pathsum="Users/aettinger/Documents/isotria/simdat_modsums/2stage/dprobs"
files<-list.files(path = pathsum)
for(j in 1:length(files)){
namesum<-files[j]
filesum<-file.path(pathsum,namesum)
ms2stage<-read.csv(filesum)
row.names(ms2stage)<-ms2stage[,1]
ms2stage<-ms2stage[,-1]
#b<-which(fvvals==as.numeric(substr(namesum,13,15)))
sV.df[j,1]<-as.numeric(substr(namesum,11,13))
sV.df[j,2]<-ms2stage$mean[1]
sV.df[j,3]<-ms2stage[1,grep("2.5",colnames(ms2stage))]
sV.df[j,4]<-ms2stage[1,grep("97.5",colnames(ms2stage))]

sF.df[j,1]<-as.numeric(substr(namesum,18,20))
sF.df[j,2]<-ms2stage$mean[2]
sF.df[j,3]<-ms2stage[2,grep("2.5",colnames(ms2stage))]
sF.df[j,4]<-ms2stage[2,grep("97.5",colnames(ms2stage))]
fV.df[j,1]<-fV
print(sV.df[j,1]);
fV.df[j,2]<-ms2stage$mean[3]
fV.df[j,3]<-ms2stage[3,grep("2.5",colnames(ms2stage))]
fV.df[j,4]<-ms2stage[3,grep("97.5",colnames(ms2stage))]
fF.df[j,1]<-fF
print(sF.df[j,1])
fF.df[j,2]<-ms2stage$mean[4]
fF.df[j,3]<-ms2stage[4,grep("2.5",colnames(ms2stage))]
fF.df[j,4]<-ms2stage[4,grep("97.5",colnames(ms2stage))]

dV.df[j,1]<-dV
dV.df[j,2]<-ms2stage$mean[5]
dV.df[j,3]<-ms2stage[5,grep("2.5",colnames(ms2stage))]
dV.df[j,4]<-ms2stage[5,grep("97.5",colnames(ms2stage))]
dF.df[j,1]<-dF
dF.df[j,2]<-ms2stage$mean[6]
dF.df[j,3]<-ms2stage[6,grep("2.5",colnames(ms2stage))]
dF.df[j,4]<-ms2stage[6,grep("97.5",colnames(ms2stage))]
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
x11(height=6, width=10)
par(mfrow=c(1,2))
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
