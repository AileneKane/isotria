#Analysis of Bill's isotria data
#Estimate vital rates separately for groups x and y
#these vital rates will be used to fit a GAM with environmental data
#setwd("~/isotria") #at usgs
setwd("~/git/isotria/analyses")
rm(list=ls()) 
options(stringsAsFactors=FALSE)
library(lme4)
#update.packages()
isoinds<-read.csv("Isotria_Stage3_2016.csv", header=T)#
isoinds<-isoinds[-which(isoinds$UniqueID=="X-01-244"),]#individual has to be removed because of errors in its monitoring
isoinds[which(isoinds$UniqueID=="X-03-236" & isoinds$Year==2002),]$TotNoStems<-0#typo in "totnostems" column- says 1, but plant was dormant
isoinds[which(isoinds$UniqueID=="Y-07-063" & isoinds$Year==2001),]$TotNoStems<-0#typo in "totnostems" column- says 1, but plant was dormant

#Add column for emergent/nonemergent
isoinds$Emerg<-NA
isoinds[which(isoinds$TotNoStems>0),]$Emerg=1
isoinds[which(isoinds$TotNoStems==0),]$Emerg=0
#Kirsi used the following emergent stages:smallveg, largeveg, flowering (i will include arrested in this- check with her if that's what she did) and fruiting, plus 2 dormant stages (dormant stage for previously small vegetative plants and dormant stage for everything else)
#i will follow these stages, which requires merging large and small fruiting plants in the stage3 column and mergeing large arrested/flowering with small arrested/flowering
#add column for Stage3, which is like Stage2 but merges  sm Arrested with sm Flowering abd Lg arrested withLg Flowering
isoinds$StageKirsi<-NA
isoinds[which(isoinds$Stage3=="ArFlLg"|isoinds$Stage3=="ArFlSm"),]$StageKirsi="fl"
isoinds[which(isoinds$Stage3=="FruitLg"|isoinds$Stage3=="FruitSm"),]$StageKirsi="fr"
#isoinds[which(isoinds$Stage3=="VegSm"),]$StageKirsi="sv"
#isoinds[which(isoinds$Stage3=="VegLg"),]$StageKirsi="lv"
isoinds[which(isoinds$Stage3=="VegLg"|isoinds$Stage3=="VegSm"),]$StageKirsi="v"#because for just groups x and y, there are not enough instantces of lv to estimate vital rates
isoinds[which(isoinds$Stage3=="D"),]$StageKirsi="d"
unique(isoinds$StageKirsi)
#need to use previous state to get 2 dormant stages
past = isoinds[isoinds$Year < 2015,]; past = past[order(past$UniqueID),]; past = past[order(past$Year),]
dim(past)# 32135 (dim(isoinds) is 32558 rows)
past$YearPrev=past$Year+1
head(past)
colnames(past)[5]<-"YearCurrent"
#Merge present and past to get fates
isoinds_prev<-merge(isoinds,past,by.x=c("UniqueID","Year","Group",  "Block", "Plant","YearFound"),by.y=c("UniqueID","YearPrev","Group",  "Block", "Plant","YearFound"),all.x = TRUE)#this includes all plants, without matches in past, as well
#delete unncessary columns (we only want stage)
isoinds_prev2=subset(isoinds_prev,select=c("UniqueID","Year", "Group",  "Block", "Plant",  "YearFound","Emerg.x","Emerg.y","StageKirsi.x","StageKirsi.y"))
#I want to add a column with a column to identify if a plant is dormant AND the last observed state was small vegetative
isoinds_prev2$dormstage=NA
isoinds_prev2[which(isoinds_prev2$StageKirsi.x=="d"),]$dormstage="d"
isoinds_prev2[which(isoinds_prev2$StageKirsi.y=="v" & isoinds_prev2$StageKirsi.x=="d"),]$dormstage="dv"
#isoinds_prev2[which(isoinds_prev2$StageKirsi.y=="lv" & isoinds_prev2$StageKirsi.x=="d"),]$dormstage="d"
isoinds_prev2[which(isoinds_prev2$StageKirsi.y=="fl" & isoinds_prev2$StageKirsi.x=="d"),]$dormstage="d"
isoinds_prev2[which(isoinds_prev2$StageKirsi.y=="fr" & isoinds_prev2$StageKirsi.x=="d"),]$dormstage="d"
unique(isoinds_prev2$dormstage)
for(i in 2:dim(isoinds_prev2)[1]){
  ind<-isoinds_prev2[i,which(colnames(isoinds_prev2)=="UniqueID")]
  currentstage<-isoinds_prev2[i,which(colnames(isoinds_prev2)=="StageKirsi.x")]
  previousstage<-isoinds_prev2[i-1,which(colnames(isoinds_prev2)=="StageKirsi.x")]
  previousdormstage<-isoinds_prev2[i-1,which(colnames(isoinds_prev2)=="dormstage")]
  currentdormstage<-isoinds_prev2[i,which(colnames(isoinds_prev2)=="dormstage")]
  if(currentstage!="d"){next}
  if(is.na(previousdormstage)) {if (currentstage=="d" & previousstage=="fl"|previousstage=="fr"){next}
    else if (currentstage=="d" & previousstage=="v"){isoinds_prev2[i,which(colnames(isoinds_prev2)=="dormstage")]<-"dv"}}
  else if(!is.na(previousdormstage)){
    if(isoinds_prev2$UniqueID[i]==ind & currentstage=="d" & previousdormstage=="dv"){isoinds_prev2[i,which(colnames(isoinds_prev2)=="dormstage")]<-"dv"}
    if(isoinds_prev2$UniqueID[i]==ind & currentstage=="d" & previousdormstage=="d"){isoinds_prev2[i,which(colnames(isoinds_prev2)=="dormstage")]<-"d"}
    if(currentstage=="d" & isoinds_prev2[i-1,]$StageKirsi.x=="v"){isoinds_prev2[i,which(colnames(isoinds_prev2)=="dormstage")]<-"dv"}
    if(currentstage=="d" & !is.na(currentdormstage)|previousdormstage!="dv"|isoinds_prev2[i-1,]$StageKirsi.x!="v"){next}
    if(isoinds_prev2$UniqueID[i]==ind & currentstage=="d" & previousdormstage=="dv"){isoinds_prev2[i,which(colnames(isoinds_prev2)=="dormstage")]<-"dv"}
  }
    else if (isoinds_prev2$UniqueID[i]==ind & currentstage=="d" & previousdormstage=="d"){isoinds_prev2[i,which(colnames(isoinds_prev2)=="dormstage")]<-"d"}
}
head(isoinds_prev2)
#rename columns
colnames(isoinds_prev2)[7:10]<-c("Emerg_cur","Emerg_prev","Stage_cur","Stage_prev")
#now want to get fates of each plant
present = isoinds_prev2[isoinds_prev2$Year <2015 & isoinds_prev2$Year >1983, ]; present = present[order(present$UniqueID),]; present = present[order(present$Year),]
dim(present)# 31837, 11 columns
future = isoinds_prev2[isoinds_prev2$Year > 1984, ]; future = future[order(future$UniqueID),]; future = future[order(future$Year),]
dim(future)#31751 rows, 11 columns
head(present);head(future)
colnames(future)
future$YearFate=future$Year-1
#Merge present and future to get fates
isofate<-merge(present,future,by.x=c("UniqueID","Year","Group",  "Block", "Plant","YearFound"),by.y=c("UniqueID","YearFate","Group",  "Block", "Plant","YearFound"),all = TRUE)#this includes all plants, without matches in past/future, as well
isofate2<-merge(present,future,by.x=c("UniqueID","Year","Group",  "Block", "Plant","YearFound"),by.y=c("UniqueID","YearFate","Group",  "Block", "Plant","YearFound"))#this only includes plants for which there is fate data that matches
dim(isofate);dim(isofate2)
head(isofate2)
colnames(isofate)
isofate=isofate[,-12]#remove year column
isofate2=isofate2[,-12]
#add column names
colnames(isofate)=c("UniqueID","Year", "Group",  "Block", "Plant",  "YearFound","Emerg_cur","Emerg_prev","Stage_cur","Stage_prev","dormstage","FutEmerg","Emerg_REMOVE","Fate","Stage_cur2","dormstage_fut")
colnames(isofate2)=c("UniqueID","Year", "Group",  "Block", "Plant",  "YearFound","Emerg_cur","Emerg_prev","Stage_cur","Stage_prev","dormstage","FutEmerg","Emerg_REMOVE","Fate","Stage_cur2","dormstage_fut")

isofate<-subset(isofate,select=c("UniqueID","Year", "Group",  "Block", "Plant",  "YearFound","Emerg_cur","Emerg_prev","Stage_cur","dormstage","Stage_prev","FutEmerg","Fate","dormstage_fut"))
isofate2<-subset(isofate2,select=c("UniqueID","Year", "Group",  "Block", "Plant",  "YearFound","Emerg_cur","Emerg_prev","Stage_cur","dormstage","Stage_prev","FutEmerg","Fate","dormstage_fut"))
#I think I will mostly use isofate2, which only includes matches, so stick with that one for rest of formatting
#Add column for survival
#To do this,loop through each individual for each row, if dormant for less than 5 years, then alive. if not, then dead

#Include both dormant stages in Stage_cur
isofate2[which(isofate2$dormstage=="dv"),]$Stage_cur<-"dv"
isofate2$Stage_cur<-as.factor(isofate2$Stage_cur)
#Add column for Large (vs small)
isofate2$tov=0
isofate2[which(isofate2$Fate=="v"),]$tov=1
#Add column for Flowering(includes arrested) and Fruiting
isofate2$tofl=0
isofate2[which(isofate2$Fate=="fl"|isofate2$Fate=="fr"),]$tofl=1
#Add column for Fruiting
isofate2$tofr=0
isofate2[which(isofate2$Fate=="fr"),]$tofr=1

##Select out just groups X and Y for this analysis
groupX=isofate2[isofate2$Group=="X",]
groupY=isofate2[isofate2$Group=="Y",]
#some individuals have error in "year found column (X-07-235 is only one in groups x or y). for this individual, it was not actually observed in 2007 (emergence=0), so remnove this row
groupX<-groupX[-which(groupX$UniqueID=="X-07-235" & groupX$Year==2007),]
head(groupX)
###Fit binomial models for each transition to get vital rate estimates
#First, a model for emergence by stage
#observed emergence by stage
numsdsaliveX<-tapply(groupX$FutEmerg,groupX$Stage_cur, sum)
numsdsX<-tapply(groupX$FutEmerg,groupX$Stage_cur, length)
propemergX<-numsdsaliveX/numsdsX
propemergX#
numsdsaliveY<-tapply(groupY$FutEmerg,groupY$Stage_cur, sum)
numsdsY<-tapply(groupY$FutEmerg,groupY$Stage_cur, length)
propemergY<-numsdsaliveY/numsdsY
propemergY#
#Model for prob of emergence (need to condition this on survival, but have to model this first for dormant stages)
xemerg.mod <- glmer(FutEmerg ~ -1 +Stage_cur + (-1+Stage_cur|Year), data=groupX, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),family=binomial(link=logit))#
yemerg.mod <- glmer(FutEmerg ~ -1 +Stage_cur + (-1+Stage_cur|Year), data=groupY, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),family=binomial(link=logit))#
summary(xemerg.mod)#midentifiable with just vegetative state (with lv and sv, not identifiable)
B.xemerg.yr=coef(xemerg.mod)#
xprobemerg = exp(B.xemerg.yr$Year)/(1+exp(B.xemerg.yr$Year))#back-transformed= these are the probabilities of emergence, 
colnames(xprobemerg)
write.csv(xprobemerg,"xprobemerg.csv")
summary(yemerg.mod)#summary(emerg.mod)
B.yemerg.yr=coef(yemerg.mod)#
yprobemerg = exp(B.yemerg.yr$Year)/(1+exp(B.yemerg.yr$Year))#back-transformed= these are the probabilities of emergence, 
colnames(yprobemerg)
write.csv(yprobemerg,"yprobemerg.csv")
#Model transitions across life stages
#2) for emergent plants, probability of flowering (which includes arrested) 
emergX=groupX[groupX$FutEmerg==1,]
emergY=groupY[groupY$FutEmerg==1,]
emergX$Stage_cur<-factor(emergX$Stage_cur)
emergY$Stage_cur<-factor(emergY$Stage_cur)
emergX$Year<-factor(emergX$Year)
emergY$Year<-factor(emergY$Year)
xflow.mod <- glmer(tofl ~ -1 +Stage_cur + (-1+Stage_cur|Year), data=emergX, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),family=binomial(link=logit))#probability of flowering conditioned on emergence
yflow.mod <- glmer(tofl ~ -1 +Stage_cur + (-1+Stage_cur|Year), data=emergY, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),family=binomial(link=logit))#probability of flowering conditioned on emergence
B.xfl.yr=coef(xflow.mod)#
xprobfl = exp(B.xfl.yr$Year)/(1+exp(B.xfl.yr$Year))#back-transformed= these are the probabilities of flowering, 
B.yfl.yr=coef(yflow.mod)#
yprobfl = exp(B.yfl.yr$Year)/(1+exp(B.yfl.yr$Year))#back-transformed= these are the probabilities of flowering, 
write.csv(xprobfl,"xprobfl.csv")
write.csv(yprobfl,"yprobfl.csv")

#3) for flowering plants, probability of fruiting:too little data for the below!
#select out only plants that flowered
flowX=emergX[emergX$tofl==1,]
flowY=emergY[emergY$tofl==1,]
flowX$Stage_cur<-factor(flowX$Stage_cur)
flowY$Stage_cur<-factor(flowY$Stage_cur)
flowX$Year<-factor(flowX$Year)
flowY$Year<-factor(flowY$Year)
xfr.mod <- glmer(tofr ~ -1 +Stage_cur + (-1+Stage_cur|Year), data=flowX, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),family=binomial(link=logit))#probability of fruiting, conditioned on flowering
yfr.mod <- glmer(tofr ~ -1 +Stage_cur + (-1+Stage_cur|Year), data=flowY, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),family=binomial(link=logit))#probability of fruiting, conditioned on flowering
#convergence error- i think because there is compelte linear separate (all d plants both flower and fruit; no dv plants fruit for group Y; )
#for group X, all dv plants fo to fr and 
#table(flowY$Stage_cur,flowY$Fate)
#table(flowX$Stage_cur,flowX$Fate)
B.xfr.yr=coef(xfr.mod)#
xprobfr = exp(B.xfr.yr$Year)/(1+exp(B.xfr.yr$Year))#back-transformed= these are the probabilities of fruiting, 
write.csv(xprobfr,"xprobfr.csv")
B.yfr.yr=coef(yfr.mod)#
yprobfr = exp(B.yfr.yr$Year)/(1+exp(B.yfr.yr$Year))#back-transformed= these are the probabilities of fruiting, 
write.csv(yprobfr,"yprobfr.csv")
