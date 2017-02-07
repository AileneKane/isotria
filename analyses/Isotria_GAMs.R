#Following Elizabeth's code to fit a GAM to isotria vital rate data, to see if climate data is associated with any vital rates
#Feb 2017
#By Ailene Ettinger
library(mgcv)
library(boot)
#setwd("C:\\Users\\ecrone01\\Box Sync\\Box Sync\\recovered files\\Ailene Isotria")
setwd("~/git/isotria/analyses")
#rm(list=ls()) 
# read in climate data
ClimDat = read.csv("isotria_clim.csv", header = T)#cleaned/summarized climatedata
# read in isotria data
FlwrDat = read.csv("probfl.csv")
colnames(FlwrDat)[1]<-"year"
# both files need year in lower case as a column for matching data sets
#ClimDat=subset(ClimDat,year>1982 & year<2015)

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
  tempAggWT<-aggLags(datC=useclim, fun=sum, meas=tmin, seg=month)#aggregating temperature by mean
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

#Loop through with different start months 1:12, just to see if that matters:

nUnits = 36
strt = 5
dat1 = makeDat(useclim = ClimDat,usevr = FlwrDat, monthStart = strt)
#str(dat1)
#head(dat1)
quartz(height=5,width=8)
par(mfrow=c(2,3))
# and now, fit a GAM for dormant(prev veg) plants
gam0b = gam(logit(Stage_curd) ~ s(lags, by=tcovar, bs="cs"),
            data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
print(summary(gam0b))
plot(gam0b,main="Dormant,nU=36,strt=5")
# and now, fit a GAM for dormant(prev veg) plants
gam0b = gam(logit(Stage_curdsv) ~ s(lags, by=tcovar, bs="cs"),
            data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
print(summary(gam0b))
plot(gam0b,main="Dormant (prev. veg),nU=36,strt=5")
# and now, fit a GAM for small vegetative plants
gam0b = gam(logit(Stage_cursv) ~ s(lags, by=tcovar, bs="cs"),
            data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
print(summary(gam0b))
plot(gam0b,main="Small Vegetative,nU=36,strt=5")
# and now, fit a GAM forlargevegetative plants
gam0b = gam(logit(Stage_curlv) ~ s(lags, by=tcovar, bs="cs"),
            data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
print(summary(gam0b))
plot(gam0b,main="Large Vegetative,nU=36,strt=5")
# and now, fit a GAM for flowering plants
gam0b = gam(logit(Stage_curfl) ~ s(lags, by=tcovar, bs="cs"),
            data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
print(summary(gam0b))
plot(gam0b,main="Flowering,nU=36,strt=5")
# and now, fit a GAM for fruiting plants
gam0b = gam(logit(Stage_curfr) ~ s(lags, by=tcovar, bs="cs"),
            data=dat1, method="GCV.Cp",gamma=1.2, family = gaussian) 
print(summary(gam0b))
plot(gam0b,main="Fruiting,nU=36,strt=5")
