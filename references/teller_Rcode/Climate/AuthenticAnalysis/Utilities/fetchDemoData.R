fetchDemoData<-function(doSpp,nUnits,WEEK){
  
  #doSpp=doSpp; nUnits=nUnits; WEEK=WEEKset;
  
  alphaList <- c(0.2208319, 0.5672076, 0.5160774, 0.3499988);
  midRings <- c(seq(1,19,2),seq(22.5,47.5,5),seq(55,145,10));
  
  #Simulation Parameters 
  monthStart=6;   
    
  #Climate Data
  DETREND=FALSE; SCALE=FALSE; 
        
  ##############################################################
  #load climate data
  ##############################################################
  climDfile=paste(dataDir,"/climateData/interpClim.csv",sep="")
  iClim=read.csv(file=climDfile)
  iClim<-iClim[,-1]
  
  #assign week of year
  iClim$week<-NA
  yearvec<-unique(iClim$year)
  for (i in 1:length(yearvec)){
    yeari<-yearvec[i]
    subs<-subset(iClim,iClim$year==yeari)
    week<-rep(1:53, each=7)
    iClim$week[iClim$year==yeari]<-week[1:nrow(subs)]  
  }
  
  ## Change units: Ppt in mm, temperatures in Celsius. 
  iClim$Ppt<-iClim$Ppt/0.039370
  iClim$Tmax<-(5/9)*(iClim$Tmax-32);
  iClim$Tmin<-(5/9)*(iClim$Tmin-32);
  iClim$Tmid<-(iClim$Tmin+iClim$Tmax)/2; 
  
  #plot(iClim$Tmid~iClim$week, col=iClim$year, xlab="Week of Year", ylab="Mean Temp.")
  #plot(iClim$Ppt~iClim$week, col=iClim$year, xlab="Week of Year", ylab="Mean Ppt.")  
  
  if(DETREND) {
    require(mgcv); par(mfrow=c(2,2))
    iClim$DOY<- 30*(iClim$month-1)+iClim$day
    
    fit=gam(Ppt~s(DOY,bs="cc"),gamma=1.4,data=iClim,method="GCV.Cp")
    iClim$Ppt <- fit$residuals; plot(fit,main="Precip");
    
    fit=gam(Tmax~s(DOY,bs="cc"),gamma=1.4,data=iClim,method="GCV.Cp")
    iClim$Tmax <- fit$residuals; plot(fit,main="Tmax");
    
    fit=gam(Tmin~s(DOY,bs="cc"),gamma=1.4,data=iClim,method="GCV.Cp")
    iClim$Tmin <- fit$residuals; plot(fit,main="Tmin");
    
    fit=gam(Tmid~s(DOY,bs="cc"),gamma=1.4,data=iClim,method="GCV.Cp")
    iClim$Tmid <- fit$residuals; plot(fit,main="Tmax");
  }
  
  if(SCALE) {
    iClim$Ppt<-scale(iClim$Ppt,center=FALSE,scale=TRUE); 
    iClim$Tmax<-scale(iClim$Tmax,center=FALSE,scale=TRUE); 
    iClim$Tmin<-scale(iClim$Tmin,center=FALSE,scale=TRUE); 
    iClim$Tmid<-scale(iClim$Tmid,center=FALSE,scale=TRUE);  
  }
  
  #######################################################################
  #Generate lags of climate signals
  #######################################################################
  aggLags<-function(datC, fun, meas, seg){
    pars <- as.list(match.call()[-1])
    measure<-datC[,as.character(pars$meas)]
    segment<-datC[,as.character(pars$seg)]
    agg<-with(datC,ftable(tapply(measure, list(datC$year, segment), fun)))
    agg<-as.data.frame(agg)
    names(agg)<-c("year",pars$seg,paste(pars$meas, pars$fun,"Now", sep=""))
    
    return(list(agg,paste(pars$fun)))
  }
  
  appLags<-function(agDatC,nUnits,name, fun){
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
  
  if(WEEK==F){
    #make monthly averages 
    #iClim<-iClim
    tempAggWT<-aggLags(datC=iClim, fun=mean, meas=Tmid, seg=month)
    subAgg<-tempAggWT[[1]]
    subAgg$month<-as.numeric(as.character(subAgg$month))
    subAgg<-subAgg[with(subAgg,order(year,month)),]
    
    tempWaveWT<-appLags(agDatC=subAgg, nUnits=nUnits, name=paste("month","Tmid", sep=""), fun=tempAggWT[[2]])
    tempWaveWT<-na.omit(tempWaveWT)
    
    monthStart=monthStart
    subsWT<-subset(tempWaveWT,tempWaveWT$month==monthStart) #choose start month
    subsWT<-subsWT[,-which(names(subsWT)=="month")] #Get rid of "month" column
    subsWT$year<-as.numeric(as.character(subsWT$year)) #make "year" readable to merge
    
    #reshape to show temp over the last nUnits
    labs<-c("year","t.00", paste("t.0",1:9, sep=""),
            paste("t.",10:nUnits, sep=""))
    names(subsWT) <- labs
    
    tempAggWP<-aggLags(datC=iClim, fun=mean, meas=Ppt, seg=month)
    subAgg<-tempAggWP[[1]]
    subAgg$month<-as.numeric(as.character(subAgg$month))
    subAgg<-subAgg[with(subAgg,order(year,month)),]
    
    tempWaveWP<-appLags(agDatC=subAgg, nUnits=nUnits, name=paste("month","Ppt", sep=""), fun=tempAggWP[[2]])
    tempWaveWP<-na.omit(tempWaveWP)
    subsWP<-subset(tempWaveWP,tempWaveWP$month==monthStart) #choose a start month
    subsWP<-subsWP[,-which(names(subsWP)=="month")] #remove month column
    subsWP$year<-as.numeric(as.character(subsWP$year)) #make "year" readable to merge
    
    labs<-c("year","p.00",paste("p.0",1:9, sep=""),
            paste("p.",10:nUnits, sep=""))
    names(subsWP) <- labs
    
    climDat<-cbind(subsWP,subsWT)
    
  }
  
  if(WEEK==T){
    #make monthly averages 
    tempAggWT<-aggLags(datC=iClim, fun=mean, meas=Tmid, seg=week)
    subAgg<-tempAggWT[[1]]
    subAgg$week<-as.numeric(as.character(subAgg$week))
    subAgg<-subAgg[with(subAgg,order(year,week)),]
    
    tempWaveWT<-appLags(agDatC=subAgg, nUnits=nUnits, name=paste("week","Tmid", sep=""), fun=tempAggWT[[2]])
    tempWaveWT<-na.omit(tempWaveWT)
    
    Start=monthStart*4
    subsWT<-subset(tempWaveWT,tempWaveWT$week==Start) #choose start month
    subsWT<-subsWT[,-which(names(subsWT)=="week")] #Get rid of "week" column
    subsWT$year<-as.numeric(as.character(subsWT$year)) #make "year" readable to merge
    
    #reshape to show temp over the last nUnits
    labs<-c("year","t.00", paste("t.0",1:9, sep=""),
            paste("t.",10:nUnits, sep=""))
    names(subsWT) <- labs
    
    tempAggWP<-aggLags(datC=iClim, fun=mean, meas=Ppt, seg=week)
    subAgg<-tempAggWP[[1]]
    subAgg$week<-as.numeric(as.character(subAgg$week))
    subAgg<-subAgg[with(subAgg,order(year,week)),]
    
    tempWaveWP<-appLags(agDatC=subAgg, nUnits=nUnits, name=paste("week","Ppt", sep=""), fun=tempAggWP[[2]])
    tempWaveWP<-na.omit(tempWaveWP)
    subsWP<-subset(tempWaveWP,tempWaveWP$week==Start) #choose a start month
    subsWP<-subsWP[,-which(names(subsWP)=="week")] #remove week column
    subsWP$year<-as.numeric(as.character(subsWP$year)) #make "year" readable to merge
    
    labs<-c("year","p.00", paste("p.0",1:9, sep=""),
            paste("p.",10:nUnits, sep=""))
    names(subsWP) <- labs
    
    climDat<-cbind(subsWP,subsWT)
    
  }
  
  ##############################################################
  #load demographic data
  ##############################################################
  # load growth data
  nonCompLength.g=5; # Number of columns in SppData that are not measures of competitors
  groDfile=paste(dataDir,"/speciesdata/",doSpp,"/growDnoNA.csv",sep="")
  groD=read.csv(file=groDfile)
  D=groD[groD$allEdge==0,]; 
  #D$year=as.factor(D$year)
  D$logarea.t0=log(D$area.t0)
  D$logarea.t1=log(D$area.t1)
  D$quad=as.character(D$quad)
  D=D[order(D$X),]
  
  ################# drop seedlings 
  e <- (D$logarea.t0>0)
  D <- D[e,];
  
  ##########################################################
  # Read in data on neighbors 
  ##########################################################
  ringD <- read.csv(paste(dataDir,"/speciesdata/",doSpp,"/",doSpp,"_nbhood_rings.csv",sep=""))
  ringD$year<-as.factor(ringD$year)
  
  # merge D with ringD (D contains fewer rows)
  D<-merge(D,ringD,by.x=c("quad","year","trackID"),by.y=c("quad","year","genetID"))
  D=D[order(D$X),]
  rm(ringD)
  row.names(D) <- NULL  
  
  ## pull out annulus data for the focal species  
  sppCols=which(substr(names(D),1,4)==doSpp); 
  sppData=cbind(D$logarea.t1,as.factor(D$quad),D$year,D$logarea.t0,as.factor(D$Group),as.matrix(D[,sppCols]))
  colnames(sppData)<-c("logarea.t1","quad","year","logarea.t0","Group", colnames(sppData)[(1+nonCompLength.g):length(colnames(sppData))])
  
  intraD<-sppData  # focal spp and intraspecific neighbors
  intraD<-data.frame(intraD)
  
  sppNum <- which(sppList==paste(doSpp)); 
  alpha <- alphaList[sppNum]; 
  dist_wts <- exp(-alpha*midRings); 
  
  C <- as.matrix(intraD[,grep(paste(doSpp),names(intraD))]) #matrix of conspecific areas in the annuli 
  W <- C%*%dist_wts; 
  dataG<-cbind(intraD[,c("logarea.t1","year","logarea.t0","Group")], W)
  
  climDat<-climDat[,-grep("year",names(climDat))[2]]
  dataG$year<-as.numeric(paste(19,dataG$year,sep=""))
  
  datag<-merge(dataG, climDat, by="year")
  datag$year<-as.factor(datag$year);
  datag$Group <- as.factor(datag$Group);
  
  ## define and enter covariates the way the fitting function wants
  pvars <- which(substr(names(datag),1,2)=="p."); 
  datag$pcovar <- as.matrix(datag[,pvars])
  
  tvars <- which(substr(names(datag),1,2)=="t."); 
  datag$tcovar <- as.matrix(datag[,tvars]) 
  
  lags <- matrix(0,nrow(datag),length(tvars)); 
  for(i in 1:ncol(lags)) lags[,i]=i; 
  datag$lags=as.matrix(lags); 
  
  data<-datag
  
  return(data)
}

