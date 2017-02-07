#Cleaning and getting the climate data in good shape for isotria
#want to get climate variables by month and year, removing duplicate years and months
rm(list=ls()) 
options(stringsAsFactors=FALSE)
allclim = read.csv("../data/climate_isotria_NH_1976_2016.csv", header = T)#These data are form multiple stations
clim=allclim[which(substr(allclim$NAME,1,8)=="LAKEPORT"),]
head(clim);dim(clim)
clim$year<-substr(clim$DATE,1,4)
clim$month<-substr(clim$DATE,6,7)
#Select the climate variables of interest:
clim2<-subset(clim,select = c(NAME,DATE,year,month,TMAX,TMIN,HTDD,CLDD,PRCP,SNOW,DSND))
tmax<-aggregate(x=subset(clim2, select=c("TMAX")), by=list(clim2$year,clim2$month), FUN=mean,na.rm=T)
tmin<-aggregate(x=subset(clim2, select=c("TMIN")), by=list(clim2$year,clim2$month), FUN=mean,na.rm=T)
htdd<-aggregate(x=subset(clim2, select=c("HTDD")), by=list(clim2$year,clim2$month), FUN=mean,na.rm=T)
cldd<-aggregate(x=subset(clim2, select=c("CLDD")), by=list(clim2$year,clim2$month), FUN=mean,na.rm=T)
prcp<-aggregate(x=subset(clim2, select=c("PRCP")), by=list(clim2$year,clim2$month), FUN=mean,na.rm=T)
snow<-aggregate(x=subset(clim2, select=c("SNOW")), by=list(clim2$year,clim2$month), FUN=mean,na.rm=T)
dsnd<-aggregate(x=subset(clim2, select=c("DSND")), by=list(clim2$year,clim2$month), FUN=mean,na.rm=T)
tmax[tmax$TMAX=="NaN",]$TMAX<-NA
tmin[tmin$TMIN=="NaN",]$TMIN<-NA
htdd[htdd$HTDD=="NaN",]$HTDD<-NA
cldd[cldd$CLDD=="NaN",]$CLDD<-NA
prcp[prcp$PRCP=="NaN",]$PRCP<-NA
snow[snow$SNOW=="NaN",]$SNOW<-NA
dsnd[dsnd$DSND=="NaN",]$DSND<-NA

##For now, replace missing data with monthly mean across all years. talk to Elizabeth about other options
clim2a<-cbind(tmax,tmin$TMIN,prcp$PRCP,snow$SNOW,dsnd$DSND,htdd$HTDD,cldd$CLDD)
colnames(clim2a)<-substr(colnames(clim2a),1,4)
colnames(clim2a)[1:3]<-c("year","month","tmax")
#for some reason, september 1982 is missing completely! add in missing data
clim3<-rbind(clim2a,c("1982","09",NA,NA,NA,NA,NA,NA,NA))
clim3$month<-as.character(clim3$month)
clim3$tmax <-as.numeric(clim3$tmax)
clim3$tmin <-as.numeric(clim3$tmin)
clim3$prcp <-as.numeric(clim3$prcp)
clim3$snow <-as.numeric(clim3$snow)
clim3$dsnd <-as.numeric(clim3$dsnd)
clim3$htdd <-as.numeric(clim3$htdd)
clim3$cldd <-as.numeric(clim3$cldd)
mos<-unique(clim3$month)
for(i in 1:12){
  print(i)
  if(length(clim3[which(is.na(clim3$tmax) & clim3$month==mos[i]),]$tmax)!=0){clim3[which(is.na(clim3$tmax) & clim3$month==mos[i]),]$tmax<-mean(as.numeric(clim3[clim3$month==mos[i],]$tmax), na.rm=T)}
  print(clim3[which(is.na(clim3$tmax) & clim3$month==mos[i]),])
  if(length(clim3[which(is.na(clim3$tmin) & clim3$month==mos[i]),]$tmin)!=0){clim3[which(is.na(clim3$tmin) & clim3$month==mos[i]),]$tmin<-mean(as.numeric(clim3[clim3$month==mos[i],]$tmin), na.rm=T)}
  if(length(clim3[which(is.na(clim3$htdd) & clim3$month==mos[i]),]$htdd)!=0){clim3[which(is.na(clim3$htdd) & clim3$month==mos[i]),]$htdd<-mean(as.numeric(clim3[clim3$month==mos[i],]$htdd), na.rm=T)}
  if(length(clim3[which(is.na(clim3$cldd) & clim3$month==mos[i]),]$htdd)!=0){clim3[which(is.na(clim3$cldd) & clim3$month==mos[i]),]$cldd<-mean(as.numeric(clim3[clim3$month==mos[i],]$cldd), na.rm=T)}
  if(length(clim3[which(is.na(clim3$prcp) & clim3$month==mos[i]),]$prcp)!=0){clim3[which(is.na(clim3$prcp) & clim3$month==mos[i]),]$prcp<-mean(as.numeric(clim3[clim3$month==mos[i],]$prcp), na.rm=T)}
  if(length(clim3[which(is.na(clim3$snow) & clim3$month==mos[i]),]$snow)!=0){clim3[which(is.na(clim3$snow) & clim3$month==mos[i]),]$snow<-mean(as.numeric(clim3[clim3$month==mos[i],]$snow), na.rm=T)}
  if(length(clim3[which(is.na(clim3$dsnd) & clim3$month==mos[i]),]$dsnd)!=0){clim3[which(is.na(clim3$dsnd) & clim3$month==mos[i]),]$dsnd<-mean(as.numeric(clim3[clim3$month==mos[i],]$dsnd), na.rm=T)}
}
clim4 <- clim3[order(clim3$year),]#order by year
write.csv(clim4,"isotria_clim.csv", row.names = FALSE)
dim(clim4)
