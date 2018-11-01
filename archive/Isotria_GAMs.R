#Following Elizabeth's code to fit a GAM to isotria vital rate data, to see if climate data is associated with any vital rates
#Feb 2017
#By Ailene Ettinger
rm(list=ls()) 

library(mgcv)
library(boot)
library(reshape)
#setwd("C:\\Users\\ecrone01\\Box Sync\\Box Sync\\recovered files\\Ailene Isotria")
setwd("~/git/isotria/analyses")
# read in climate data
ClimDat = read.csv("isotria_clim.csv", header = T)#cleaned/summarized climatedata
# read in isotria data
FlwrDat = read.csv("probfr.csv")
colnames(FlwrDat)[1]<-"year"
# both files need year in lower case as a column for matching data sets
ClimDat=subset(ClimDat,year>1980 & year<2015)
colnames(ClimDat)
pairs(ClimDat[,3:9])
# function from Teller 2016 to aggregate data by month-my data are already in this format
aggLags<-function(datC, fun, meas, seg){
  pars <- as.list(match.call()[-1])
  measure<-datC[,as.character(pars$meas)]
  segment<-datC[,as.character(pars$seg)]
  agg<-with(datC,ftable(tapply(measure, list(datC$year, segment), fun)))
  agg<-as.data.frame(agg)
  names(agg)<-c("year",pars$seg,paste(pars$meas, pars$fun,"Now", sep=""))
  
  return(list(agg,paste(pars$fun)))
}

# function from Teller 2016 to create lagged variables
appLags<-function(agDatC,nUnits,name,fun){
  storeNames<-names(agDatC)
  for (i in 1:nUnits){
    agDatC<-as.data.frame(agDatC)
    waved<-agDatC[[3]]
    newWave<-c(rep(NA, times=i),waved) 
    newThing<-newWave[1:nrow(agDatC)]
    agDatC<-cbind(agDatC, newThing)
  }
  names(agDatC)<-c(storeNames,paste(name,fun,1:nUnits, sep=""))
  return(agDatC)
}

makeDat = function(useclim, usevr, monthStart){
  tempAggWT<-aggLags(datC=useclim, fun=mean, meas=tmax, seg=month)#aggregating temperature by mean
  subAgg<-tempAggWT[[1]]
  subAgg$month<-as.numeric(as.character(subAgg$month))
  subAgg<-subAgg[with(subAgg,order(year,month)),]

  tempWaveWT<-appLags(agDatC=subAgg, nUnits=nUnits, name=paste("month","temp", sep=""), fun=tempAggWT[[2]])
  tempWaveWT<-na.omit(tempWaveWT)

  monthStart=monthStart
  subsWT<-subset(tempWaveWT,tempWaveWT$month==monthStart) #choose start month
  subsWT<-subsWT[,-which(names(subsWT)=="month")] #Get rid of "month" column
  subsWT$year<-as.numeric(as.character(subsWT$year)) #make "year" readable to merge

  #reshape to show temp over the last nUnits
  labs<-c("year","t.00", paste("t.0",1:9, sep=""),
        paste("t.",10:nUnits, sep=""))
  names(subsWT) <- labs

  # merge climate and isotria data, following code from Teller 2016
  datag<-merge(usevr, subsWT, by="year")
  datag$year<-as.factor(datag$year);

  ## define and enter covariates the way the fitting function wants
  tvars <- which(substr(names(datag),1,2)=="t."); 
  datag$tcovar <- as.matrix(datag[,tvars]) 

  lags <- matrix(0,nrow(datag),length(tvars)); 
  for(i in 1:ncol(lags)) lags[,i]=i; 
  datag$lags=as.matrix(lags); 
  datag$yr = strtoi(datag$year)

  return(datag)
}

nUnits = 36
strt = 5
dat1 = makeDat(useclim = ClimDat,usevr = FlwrDat, monthStart = strt)
#str(dat1)
#head(dat1)
dat2<-cbind(dat1[2:7],dat1$year,dat1[8:47])
dat3<-reshape(dat2,varying=list(names(dat2)[1:6]),direction = "long", v.names = c("probfruit"), times = c(names(dat2)[1:6]))
dat4<-dat3[,-44]
colnames(dat4)[42]<-"Stage"
colnames(dat4)[1]<-"year"

quartz(height=5,width=8)
par(mfrow=c(2,3))
# and now, fit a GAM for dormant(prev veg) plants
#try comparing model fit with stages and without stages
#Model that includes stage:
mod1 <- gam(logit(probfruit) ~ s(lags, by=tcovar, bs="cs") + s(Stage, bs = "cs", by =tcovar) + te(lags, Stage, by =tcovar, bs = rep("cs",2)),data=dat4, method="GCV.Cp",gamma=1.2, family = gaussian)

#No Stage-level differences in effect of lag
  mod0 <- gam(logit(probfruit) ~ Stage + s(lags, by=tcovar, bs="cs") + 
                s(Doy, bs = "cc", by = Loc, k = 5, m = 1) + 
                s(Tod, bs = "cc", k = 5) + 
                s(Tod, bs = "cc", by = Loc, k = 5, m = 1),
              data = dat4, method = "ML")




gam0a = gam(logit(probfruit) ~ s(lags, by=tcovar, bs="cs"),
            data=dat4, method="GCV.Cp",gamma=1.2, family = gaussian) 
print(summary(gam0a))
gam0a = gam(logit(probfruit) ~ Stage+ s(lags, by=tcovar, bs="cs"),
            data=dat4, method="GCV.Cp",gamma=1.2, family = gaussian) 
print(summary(gam0b))

AIC(gam0a,gam0b)


layout(matrix(1:2, ncol = 2))
plot(gam0b)
layout(1)
acf(resid(gam0b), lag.max = 36, main = "ACF")
pacf(resid(gam0b), lag.max = 36, main = "ACF")


gam1b = gam(logit(Stage_curd) ~ as.numeric(year)+s(lags, by=tcovar, bs="cs"),
            data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
print(summary(gam1b))
AIC(gam1b,gam0b)
layout(matrix(1:2, ncol = 2))
plot(gam1b)
layout(1)

plot(gam0b,main=paste("stage=",substr(summary(gam0b)$formula,16,17)[2],",nmos=",nUnits,"strt=",strt))
plot(gam1b,main=paste("stage=",substr(summary(gam0b)$formula,16,17)[2],",nmos=",nUnits,"strt=",strt))

# and now, fit a GAM for dormant(prev veg) plants
gam0b = gam(logit(Stage_curdsv) ~ s(lags, by=tcovar, bs="cs"),
            data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
print(summary(gam0b))
plot(gam0b,main=paste("stage=",substr(summary(gam0b)$formula,16,17)[2],",nmos=",nUnits,"strt=",strt))
# and now, fit a GAM for small vegetative plants
gam0b = gam(logit(Stage_cursv) ~ s(lags, by=tcovar, bs="cs"),
            data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
print(summary(gam0b))
plot(gam0b,main=paste("stage=",substr(summary(gam0b)$formula,16,17)[2],",nmos=",nUnits,"strt=",strt))
# and now, fit a GAM for large vegetative plants
gam0b = gam(logit(Stage_curlv) ~ s(lags, by=tcovar, bs="cs"),
            data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
print(summary(gam0b))
plot(gam0b,main=paste("stage=",substr(summary(gam0b)$formula,16,17)[2],",nmos=",nUnits,"strt=",strt))
# and now, fit a GAM for flowering plants
gam0b = gam(logit(Stage_curfl) ~ s(lags, by=tcovar, bs="cs"),
            data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
print(summary(gam0b))
plot(gam0b,main=paste("stage=",substr(summary(gam0b)$formula,16,17)[2],",nmos=",nUnits,"strt=",strt))
# and now, fit a GAM for fruiting plants
gam0b = gam(logit(Stage_curfr) ~ s(lags, by=tcovar, bs="cs"),
            data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
print(summary(gam0b))
plot(gam0b,main=paste("stage=",substr(summary(gam0b)$formula,16,17)[2],",nmos=",nUnits,"strt=",strt))

#Write a loop to fit models with different climate variables and different life stages to compare AICs
#want to save lag used, AIC, Dev explained & climate variab
stages<-c("Stage_curd","Stage_curdsv","Stage_curfl","Stage_curfr","Stage_curlv","Stage_cursv") 
lags_months<-c(24,36)

modresults.lags<-c()
for (j in 1:length(lags_months)){
  nUnits=lags_months[j]
  dat1 = makeDat(useclim = ClimDat,usevr = FlwrDat, monthStart = strt)
  modresults<-c()
for (i in 1:length(stages)){
  resp=logit(dat1[,i+1])
  dat1 = makeDat(useclim = ClimDat,usevr = FlwrDat, monthStart = strt)
  mod1=gam(resp ~ s(lags, by=tcovar, bs="cs"),
           data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
  mod2=gam(resp ~ as.numeric(year)+s(lags, by=tcovar, bs="cs"),
           data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
  modsums<-c(stages[i],lags_months[j],strt,"tmin",summary(mod1)$s.pv,summary(mod1)$r.sq,summary(mod1)$dev.expl,AIC(mod1),"year+tmin",summary(mod2)$s.pv,summary(mod2)$r.sq,summary(mod2)$dev.expl,AIC(mod2))
modresults=rbind(modresults,modsums)
}
modresults.lags=rbind(modresults.lags,modresults)
}
  colnames(modresults.lags)<-c("stage","lag_mos","start.mos","mod1","mod1.pval","mod1.rsq","mod1.dev.expl","mod1.aic","mod2","mod2.pval","mod2.rsq","mod2.dev.expl","mod1.aic")

  mod.results.lags_tmin<-modresults.lags

  #rerun code at top of script with tmax-i need to redo this so i don't have to do so much by hand!!!
modresults.lags<-c()
  for (j in 1:length(lags_months)){
    nUnits=lags_months[j]
    dat1 = makeDat(useclim = ClimDat,usevr = FlwrDat, monthStart = strt)
    modresults<-c()
    for (i in 1:length(stages)){
      resp=logit(dat1[,i+1])
      dat1 = makeDat(useclim = ClimDat,usevr = FlwrDat, monthStart = strt)
      mod1=gam(resp ~ s(lags, by=tcovar, bs="cs"),
               data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
      mod2=gam(resp ~ as.numeric(year)+s(lags, by=tcovar, bs="cs"),
               data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
      modsums<-c(stages[i],lags_months[j],strt,"tmax",summary(mod1)$s.pv,summary(mod1)$r.sq,summary(mod1)$dev.expl,AIC(mod1),"year+tmax",summary(mod2)$s.pv,summary(mod2)$r.sq,summary(mod2)$dev.expl,AIC(mod2))
      modresults=rbind(modresults,modsums)
    }
    modresults.lags=rbind(modresults.lags,modresults)
  }

colnames(modresults.lags)<-c("stage","lag_mos","start.mos","mod1","mod1.pval","mod1.rsq","mod1.dev.expl","mod1.aic","mod2","mod2.pval","mod2.rsq","mod2.dev.expl","mod1.aic")
mod.results.lags_tmax<-modresults.lags

#rerun code at top of script with prcp-i need to redo this so i don't have to do so much by hand!!!
modresults.lags<-c()
for (j in 1:length(lags_months)){
  nUnits=lags_months[j]
  dat1 = makeDat(useclim = ClimDat,usevr = FlwrDat, monthStart = strt)
  
  modresults<-c()
  for (i in 1:length(stages)){
    resp=logit(dat1[,i+1])
    dat1 = makeDat(useclim = ClimDat,usevr = FlwrDat, monthStart = strt)
    mod1=gam(resp ~ s(lags, by=tcovar, bs="cs"),
             data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
    mod2=gam(resp ~ as.numeric(year)+s(lags, by=tcovar, bs="cs"),
             data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
    modsums<-c(stages[i],lags_months[j],strt,"prcp",summary(mod1)$s.pv,summary(mod1)$r.sq,summary(mod1)$dev.expl,AIC(mod1),"year+prcp",summary(mod2)$s.pv,summary(mod2)$r.sq,summary(mod2)$dev.expl,AIC(mod2))
    modresults=rbind(modresults,modsums)
  }
  modresults.lags=rbind(modresults.lags,modresults)
}

colnames(modresults.lags)<-c("stage","lag_mos","start.mos","mod1","mod1.pval","mod1.rsq","mod1.dev.expl","mod1.aic","mod2","mod2.pval","mod2.rsq","mod2.dev.expl","mod1.aic")
mod.results.lags_prcp<-modresults.lags

#rerun code at top of script with snow-i need to redo this so i don't have to do so much by hand!!!
modresults.lags<-c()
for (j in 1:length(lags_months)){
  nUnits=lags_months[j]
  dat1 = makeDat(useclim = ClimDat,usevr = FlwrDat, monthStart = strt)
  modresults<-c()
  for (i in 1:length(stages)){
    resp=logit(dat1[,i+1])
    mod1=gam(resp ~ s(lags, by=tcovar, bs="cs"),
             data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
    mod2=gam(resp ~ as.numeric(year)+s(lags, by=tcovar, bs="cs"),
             data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
    modsums<-c(stages[i],lags_months[j],strt,"snow",summary(mod1)$s.pv,summary(mod1)$r.sq,summary(mod1)$dev.expl,AIC(mod1),"year+snow",summary(mod2)$s.pv,summary(mod2)$r.sq,summary(mod2)$dev.expl,AIC(mod2))
    modresults=rbind(modresults,modsums)
  }
  modresults.lags=rbind(modresults.lags,modresults)
}
colnames(modresults.lags)<-c("stage","lag_mos","start.mos","mod1","mod1.pval","mod1.rsq","mod1.dev.expl","mod1.aic","mod2","mod2.pval","mod2.rsq","mod2.dev.expl","mod1.aic")
mod.results.lags_snow<-modresults.lags

mod.results.lags_all<-as.data.frame(rbind(mod.results.lags_tmin,mod.results.lags_tmax,mod.results.lags_prcp,mod.results.lags_snow))
mod.results.lags_all <- mod.results.lags_all[order(mod.results.lags_all$stage,mod.results.lags_all$mod1.aic),]#order by stage, then AIC

write.csv(mod.results.lags_all,"prob_emerg.gammods.csv", row.names = FALSE)

###Probability of fruiting
stages<-c("Stage_curd","Stage_curdsv","Stage_curfl","Stage_curfr","Stage_curlv","Stage_cursv") 
lags_months<-c(24,36)

modresults.lags<-c()
for (j in 1:length(lags_months)){
  nUnits=lags_months[j]
  dat1 = makeDat(useclim = ClimDat,usevr = FlwrDat, monthStart = strt)
  modresults<-c()
  for (i in 1:length(stages)){
    resp=logit(dat1[,i+1])
    dat1 = makeDat(useclim = ClimDat,usevr = FlwrDat, monthStart = strt)
    mod1=gam(resp ~ s(lags, by=tcovar, bs="cs"),
             data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
    mod2=gam(resp ~ as.numeric(year)+s(lags, by=tcovar, bs="cs"),
             data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
    modsums<-c(stages[i],lags_months[j],strt,"tmin",summary(mod1)$s.pv,summary(mod1)$r.sq,summary(mod1)$dev.expl,AIC(mod1),"year+tmin",summary(mod2)$s.pv,summary(mod2)$r.sq,summary(mod2)$dev.expl,AIC(mod2))
    modresults=rbind(modresults,modsums)
  }
  modresults.lags=rbind(modresults.lags,modresults)
}
colnames(modresults.lags)<-c("stage","lag_mos","start.mos","mod1","mod1.pval","mod1.rsq","mod1.dev.expl","mod1.aic","mod2","mod2.pval","mod2.rsq","mod2.dev.expl","mod1.aic")

mod.results.lags_tmin<-modresults.lags

#rerun code at top of script with tmax-i need to redo this so i don't have to do so much by hand!!!
modresults.lags<-c()
for (j in 1:length(lags_months)){
  nUnits=lags_months[j]
  dat1 = makeDat(useclim = ClimDat,usevr = FlwrDat, monthStart = strt)
  modresults<-c()
  for (i in 1:length(stages)){
    resp=logit(dat1[,i+1])
    dat1 = makeDat(useclim = ClimDat,usevr = FlwrDat, monthStart = strt)
    mod1=gam(resp ~ s(lags, by=tcovar, bs="cs"),
             data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
    mod2=gam(resp ~ as.numeric(year)+s(lags, by=tcovar, bs="cs"),
             data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
    modsums<-c(stages[i],lags_months[j],strt,"tmax",summary(mod1)$s.pv,summary(mod1)$r.sq,summary(mod1)$dev.expl,AIC(mod1),"year+tmax",summary(mod2)$s.pv,summary(mod2)$r.sq,summary(mod2)$dev.expl,AIC(mod2))
    modresults=rbind(modresults,modsums)
  }
  modresults.lags=rbind(modresults.lags,modresults)
}

colnames(modresults.lags)<-c("stage","lag_mos","start.mos","mod1","mod1.pval","mod1.rsq","mod1.dev.expl","mod1.aic","mod2","mod2.pval","mod2.rsq","mod2.dev.expl","mod1.aic")
mod.results.lags_tmax<-modresults.lags

#rerun code at top of script with prcp-i need to redo this so i don't have to do so much by hand!!!
modresults.lags<-c()
for (j in 1:length(lags_months)){
  nUnits=lags_months[j]
  dat1 = makeDat(useclim = ClimDat,usevr = FlwrDat, monthStart = strt)
  
  modresults<-c()
  for (i in 1:length(stages)){
    resp=logit(dat1[,i+1])
    dat1 = makeDat(useclim = ClimDat,usevr = FlwrDat, monthStart = strt)
    mod1=gam(resp ~ s(lags, by=tcovar, bs="cs"),
             data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
    mod2=gam(resp ~ as.numeric(year)+s(lags, by=tcovar, bs="cs"),
             data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
    modsums<-c(stages[i],lags_months[j],strt,"prcp",summary(mod1)$s.pv,summary(mod1)$r.sq,summary(mod1)$dev.expl,AIC(mod1),"year+prcp",summary(mod2)$s.pv,summary(mod2)$r.sq,summary(mod2)$dev.expl,AIC(mod2))
    modresults=rbind(modresults,modsums)
  }
  modresults.lags=rbind(modresults.lags,modresults)
}

colnames(modresults.lags)<-c("stage","lag_mos","start.mos","mod1","mod1.pval","mod1.rsq","mod1.dev.expl","mod1.aic","mod2","mod2.pval","mod2.rsq","mod2.dev.expl","mod1.aic")
mod.results.lags_prcp<-modresults.lags

#rerun code at top of script with snow-i need to redo this so i don't have to do so much by hand!!!
modresults.lags<-c()
for (j in 1:length(lags_months)){
  nUnits=lags_months[j]
  dat1 = makeDat(useclim = ClimDat,usevr = FlwrDat, monthStart = strt)
  modresults<-c()
  for (i in 1:length(stages)){
    resp=logit(dat1[,i+1])
    mod1=gam(resp ~ s(lags, by=tcovar, bs="cs"),
             data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
    mod2=gam(resp ~ as.numeric(year)+s(lags, by=tcovar, bs="cs"),
             data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
    modsums<-c(stages[i],lags_months[j],strt,"snow",summary(mod1)$s.pv,summary(mod1)$r.sq,summary(mod1)$dev.expl,AIC(mod1),"year+snow",summary(mod2)$s.pv,summary(mod2)$r.sq,summary(mod2)$dev.expl,AIC(mod2))
    modresults=rbind(modresults,modsums)
  }
  modresults.lags=rbind(modresults.lags,modresults)
}
colnames(modresults.lags)<-c("stage","lag_mos","start.mos","mod1","mod1.pval","mod1.rsq","mod1.dev.expl","mod1.aic","mod2","mod2.pval","mod2.rsq","mod2.dev.expl","mod1.aic")
mod.results.lags_snow<-modresults.lags

mod.results.lags_all<-as.data.frame(rbind(mod.results.lags_tmin,mod.results.lags_tmax,mod.results.lags_prcp,mod.results.lags_snow))
mod.results.lags_all <- mod.results.lags_all[order(mod.results.lags_all$stage,mod.results.lags_all$mod1.aic),]#order by stage, then AIC

write.csv(mod.results.lags_all,"prob_fruit.gammods.csv", row.names = FALSE)




#now do same thing for flow and fru prob.
#ProbFlow- these are not significantly affected by any climate variable
modresults.lags<-c()
for (j in 1:length(lags_months)){
  nUnits=lags_months[j]
  dat1 = makeDat(useclim = ClimDat,usevr = FlwrDat, monthStart = strt)
  modresults<-c()
  for (i in 1:length(stages)){
    resp=logit(dat1[,i+1])
    dat1 = makeDat(useclim = ClimDat,usevr = FlwrDat, monthStart = strt)
    mod1=gam(resp ~ s(lags, by=tcovar, bs="cs"),
             data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
    mod2=gam(resp ~ as.numeric(year)+s(lags, by=tcovar, bs="cs"),
             data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
    modsums<-c(stages[i],lags_months[j],strt,"tmin",summary(mod1)$s.pv,summary(mod1)$r.sq,summary(mod1)$dev.expl,AIC(mod1),"year+tmin",summary(mod2)$s.pv,summary(mod2)$r.sq,summary(mod2)$dev.expl,AIC(mod2))
    modresults=rbind(modresults,modsums)
  }
  modresults.lags=rbind(modresults.lags,modresults)
}
colnames(modresults.lags)<-c("stage","lag_mos","start.mos","mod1","mod1.pval","mod1.rsq","mod1.dev.expl","mod1.aic","mod2","mod2.pval","mod2.rsq","mod2.dev.expl","mod1.aic")

mod.results.lags_tmin<-modresults.lags

#rerun code at top of script with tmax-i need to redo this so i don't have to do so much by hand!!!
modresults.lags<-c()
for (j in 1:length(lags_months)){
  nUnits=lags_months[j]
  dat1 = makeDat(useclim = ClimDat,usevr = FlwrDat, monthStart = strt)
  modresults<-c()
  for (i in 1:length(stages)){
    resp=logit(dat1[,i+1])
    dat1 = makeDat(useclim = ClimDat,usevr = FlwrDat, monthStart = strt)
    mod1=gam(resp ~ s(lags, by=tcovar, bs="cs"),
             data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
    mod2=gam(resp ~ as.numeric(year)+s(lags, by=tcovar, bs="cs"),
             data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
    modsums<-c(stages[i],lags_months[j],strt,"tmax",summary(mod1)$s.pv,summary(mod1)$r.sq,summary(mod1)$dev.expl,AIC(mod1),"year+tmax",summary(mod2)$s.pv,summary(mod2)$r.sq,summary(mod2)$dev.expl,AIC(mod2))
    modresults=rbind(modresults,modsums)
  }
  modresults.lags=rbind(modresults.lags,modresults)
}

colnames(modresults.lags)<-c("stage","lag_mos","start.mos","mod1","mod1.pval","mod1.rsq","mod1.dev.expl","mod1.aic","mod2","mod2.pval","mod2.rsq","mod2.dev.expl","mod1.aic")
mod.results.lags_tmax<-modresults.lags

#rerun code at top of script with prcp-i need to redo this so i don't have to do so much by hand!!!
modresults.lags<-c()
for (j in 1:length(lags_months)){
  nUnits=lags_months[j]
  dat1 = makeDat(useclim = ClimDat,usevr = FlwrDat, monthStart = strt)
  
  modresults<-c()
  for (i in 1:length(stages)){
    resp=logit(dat1[,i+1])
    dat1 = makeDat(useclim = ClimDat,usevr = FlwrDat, monthStart = strt)
    mod1=gam(resp ~ s(lags, by=tcovar, bs="cs"),
             data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
    mod2=gam(resp ~ as.numeric(year)+s(lags, by=tcovar, bs="cs"),
             data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
    modsums<-c(stages[i],lags_months[j],strt,"prcp",summary(mod1)$s.pv,summary(mod1)$r.sq,summary(mod1)$dev.expl,AIC(mod1),"year+prcp",summary(mod2)$s.pv,summary(mod2)$r.sq,summary(mod2)$dev.expl,AIC(mod2))
    modresults=rbind(modresults,modsums)
  }
  modresults.lags=rbind(modresults.lags,modresults)
}

colnames(modresults.lags)<-c("stage","lag_mos","start.mos","mod1","mod1.pval","mod1.rsq","mod1.dev.expl","mod1.aic","mod2","mod2.pval","mod2.rsq","mod2.dev.expl","mod1.aic")
mod.results.lags_prcp<-modresults.lags

#rerun code at top of script with snow-i need to redo this so i don't have to do so much by hand!!!
modresults.lags<-c()
for (j in 1:length(lags_months)){
  nUnits=lags_months[j]
  dat1 = makeDat(useclim = ClimDat,usevr = FlwrDat, monthStart = strt)
  modresults<-c()
  for (i in 1:length(stages)){
    resp=logit(dat1[,i+1])
    mod1=gam(resp ~ s(lags, by=tcovar, bs="cs"),
             data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
    mod2=gam(resp ~ as.numeric(year)+s(lags, by=tcovar, bs="cs"),
             data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
    modsums<-c(stages[i],lags_months[j],strt,"snow",summary(mod1)$s.pv,summary(mod1)$r.sq,summary(mod1)$dev.expl,AIC(mod1),"year+snow",summary(mod2)$s.pv,summary(mod2)$r.sq,summary(mod2)$dev.expl,AIC(mod2))
    modresults=rbind(modresults,modsums)
  }
  modresults.lags=rbind(modresults.lags,modresults)
}
colnames(modresults.lags)<-c("stage","lag_mos","start.mos","mod1","mod1.pval","mod1.rsq","mod1.dev.expl","mod1.aic","mod2","mod2.pval","mod2.rsq","mod2.dev.expl","mod2.aic")
mod.results.lags_snow<-modresults.lags

mod.results.lags_all<-as.data.frame(rbind(mod.results.lags_tmin,mod.results.lags_tmax,mod.results.lags_prcp,mod.results.lags_snow))
mod.results.lags_all <- mod.results.lags_all[order(mod.results.lags_all$stage,mod.results.lags_all$mod1.aic),]#order by stage, then AIC

write.csv(mod.results.lags_all,"prob_flow.gammods.csv", row.names = FALSE)
