#Isotria multistate model
#Simulate data with known parameters to test the ability of my model to recover them
#start with simple dataset with no variation between groups x and y
#first, define mean survival, transitions, and emergence probability, 
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

