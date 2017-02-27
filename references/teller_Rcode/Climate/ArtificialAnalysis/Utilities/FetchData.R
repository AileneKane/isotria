FetchData<-function(paramList){
  
  #Directories 
  root=ifelse(.Platform$OS.type=="windows","c:/repos","~/repos"); # modify as needed
  workDir=paste(root,"/drivers/Manuscripts/Sheffield/Rcode",sep=""); 
  dataDir=paste(workDir,"/Data",sep="");
  utilDir=paste(workDir,"/Climate/ArtificialAnalysis/Utilities/",sep="")
  simDatDir=paste(workDir,"/Climate/ArtificialAnalysis/SimData/",sep="")
  
  sppList=sort(c("ARTR","HECO"))
  alphaList <- c(0.2208319, 0.5672076, 0.5160774, 0.3499988);
  midRings <- c(seq(1,19,2),seq(22.5,47.5,5),seq(55,145,10));
  
  #Simulation Parameters 
  monthStart=paramList$monthStart;  nYrs=paramList$nYrs; 
  nMonths=paramList$nMonths; WEEK=paramList$WEEK;
  nFyears=paramList$nFyears; nIndiv=paramList$nIndiv;
  
  #Climate Data
  DETREND=paramList$DETREND; SCALE=paramList$SCALE; 
  noiseSD_logP=paramList$noiseSD_logP; noiseSD_T= paramList$noiseSD_T; 
  tauHalf.T = paramList$tauHalf.T; tauHalf.P = paramList$tauHalf.P; rPT=paramList$rPT;
  
  #Demographic Data
  doSpp=paramList$doSpp;
    
  ## Constants that affect simulated response variables
  nEffect=paramList$nEffect; sigmaFrac=paramList$sigmaFrac; 
  fPrecip=paramList$fPrecip; varFrac=paramList$varFrac;
  
  
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
  
  #############################################################
  #Generate new years of climate data to extend series
  #############################################################
  fClim <- iClim;  
  fClim$year <- fClim$year-min(fClim$year)+max(iClim$year)+1; 
  fClim <- rbind(iClim,fClim); 
  
  ### AR(1) coefficients to give desired correlation decay 
  a.T <- 0.5^(1/tauHalf.T);
  a.P <- 0.5^(1/tauHalf.P);
  
  ### Compute the innovation correlation corresponding to rPT
  rho <- rPT* (1-a.T*a.P)/sqrt( (1-a.T^2)*(1-a.P^2) ) 
  if(abs(rho)>1) stop("Error: impossible P-T correlation requested") 
  
  ### generate the AR(1) time series x(t) and y(t) 
  nBurn=25*365; 
  nt <- nrow(fClim) + nBurn; 
  Xt = Yt = rep(0,nt); 
  Z = mvrnorm(nt,mu=c(0,0),Sigma=matrix(c(1,rho,rho,1),nrow=2,byrow=TRUE)); 
  sigma.1 = sqrt(1-a.T^2); sigma.2 = sqrt(1-a.P^2)
  for(j in 2:nt) {
    Xt[j]=a.T*Xt[j-1]+sigma.1*Z[j-1,1];
    Yt[j]=a.P*Yt[j-1]+sigma.2*Z[j-1,2];
  }
  Xt=Xt[-(1:nBurn)]; Yt=Yt[-(1:nBurn)];   
  
  #### make up fake Tmid: trend plus AR(1) 
  DOY = 30.37*(iClim$month-1)+iClim$day; 
  TmidFit=gam(iClim$Tmid~s(DOY),data=iClim);  
  fClim$Tmid <- rep(TmidFit$fitted,2) + noiseSD_T*Xt; 
  
  ### probability of rain as a function of Day of Year 
  fit=gam( (Ppt>0)~s(DOY),data=iClim,family="binomial"); 
  px=1:365; py=predict(fit,newdata=data.frame(DOY=px),type="response"); 
  pRain=py[round(DOY)]; 
  
  ### log of Rain when it happens
  PosPpt = iClim$Ppt; e = PosPpt>0; 
  y=PosPpt[e]; x=DOY[e]; 
  fit=gam(log(y)~s(x));  
  px=1:365; py=predict(fit,newdata=data.frame(x=px),type="response"); 
  meanLogRain = py[round(DOY)] 
  
  fakeRain = exp(rep(meanLogRain,2) + noiseSD_logP*Yt); 
  doesRain = as.numeric(runif(length(fClim$Ppt))<rep(pRain,2)) 
  fClim$Ppt = fakeRain*doesRain; 
    
  
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
    #iClim<-fClim
    tempAggWT<-aggLags(datC=fClim, fun=mean, meas=Tmid, seg=month)
    subAgg<-tempAggWT[[1]]
    subAgg$month<-as.numeric(as.character(subAgg$month))
    subAgg<-subAgg[with(subAgg,order(year,month)),]
    
    tempWaveWT<-appLags(agDatC=subAgg, nUnits=nMonths, name=paste("month","Tmid", sep=""), fun=tempAggWT[[2]])
    tempWaveWT<-na.omit(tempWaveWT)
    
    monthStart=monthStart
    subsWT<-subset(tempWaveWT,tempWaveWT$month==monthStart) #choose start month
    subsWT<-subsWT[,-which(names(subsWT)=="month")] #Get rid of "month" column
    subsWT$year<-as.numeric(as.character(subsWT$year)) #make "year" readable to merge
    
    #reshape to show temp over the last 24 months
    labs<-c("year","t.00", paste("t.0",1:9, sep=""),
            paste("t.",10:nMonths, sep=""))
    names(subsWT) <- labs
    
    tempAggWP<-aggLags(datC=fClim, fun=mean, meas=Ppt, seg=month)
    subAgg<-tempAggWP[[1]]
    subAgg$month<-as.numeric(as.character(subAgg$month))
    subAgg<-subAgg[with(subAgg,order(year,month)),]
    
    tempWaveWP<-appLags(agDatC=subAgg, nUnits=nMonths, name=paste("month","Ppt", sep=""), fun=tempAggWP[[2]])
    tempWaveWP<-na.omit(tempWaveWP)
    subsWP<-subset(tempWaveWP,tempWaveWP$month==monthStart) #choose a start month
    subsWP<-subsWP[,-which(names(subsWP)=="month")] #remove month column
    subsWP$year<-as.numeric(as.character(subsWP$year)) #make "year" readable to merge
    
    labs<-c("year","p.00",paste("p.0",1:9, sep=""),
            paste("p.",10:nMonths, sep=""))
    names(subsWP) <- labs
    
    climDat<-cbind(subsWP,subsWT)
    
  }
  
  if(WEEK==T){
    #make monthly averages 
    tempAggWT<-aggLags(datC=fClim, fun=mean, meas=Tmid, seg=week)
    subAgg<-tempAggWT[[1]]
    subAgg$week<-as.numeric(as.character(subAgg$week))
    subAgg<-subAgg[with(subAgg,order(year,week)),]
    
    tempWaveWT<-appLags(agDatC=subAgg, nUnits=52*(nYrs), name=paste("week","Tmid", sep=""), fun=tempAggWT[[2]])
    tempWaveWT<-na.omit(tempWaveWT)
    
    Start=monthStart*4
    subsWT<-subset(tempWaveWT,tempWaveWT$week==Start) #choose start month
    subsWT<-subsWT[,-which(names(subsWT)=="week")] #Get rid of "week" column
    subsWT$year<-as.numeric(as.character(subsWT$year)) #make "year" readable to merge
    
    #reshape to show temp over the last 24 months
    labs<-c("year","t.00", paste("t.0",1:9, sep=""),
            paste("t.",10:(52*nYrs), sep=""))
    names(subsWT) <- labs
    
    tempAggWP<-aggLags(datC=fClim, fun=mean, meas=Ppt, seg=week)
    subAgg<-tempAggWP[[1]]
    subAgg$week<-as.numeric(as.character(subAgg$week))
    subAgg<-subAgg[with(subAgg,order(year,week)),]
    
    tempWaveWP<-appLags(agDatC=subAgg, nUnits=52*nYrs, name=paste("week","Ppt", sep=""), fun=tempAggWP[[2]])
    tempWaveWP<-na.omit(tempWaveWP)
    subsWP<-subset(tempWaveWP,tempWaveWP$week==Start) #choose a start month
    subsWP<-subsWP[,-which(names(subsWP)=="week")] #remove week column
    subsWP$year<-as.numeric(as.character(subsWP$year)) #make "year" readable to merge
    
    labs<-c("year","p.00", paste("p.0",1:9, sep=""),
            paste("p.",10:(52*nYrs), sep=""))
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
  
  dataG$year<-as.numeric(paste(19,dataG$year,sep=""))
  
  #sample a number of individuals
  dataG<-dataG[sample(1:nrow(dataG),size=nIndiv, replace=F),]
  
  #sample a number of climate years
  sampY<-sample(unique(climDat$year), size=nFyears, replace=F)
  
  #assign new year names to data
  datag<-dataG
  for (i in 1:nFyears){
    datagSb<-dataG
    datagSb$year<-sampY[i]
    datag<-if(i==1)(datagSb)else{rbind(datag,datagSb)}}
  
  climSub<-subset(climDat,climDat$year %in% sampY) #OR sample a sequence of years of size nFyears?
  
  datag<-merge(datag, climSub)
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
  
  ##########################################################
  #Augment demographic response to climate
  ##########################################################
  datag_all<-datag
  nlags = ncol(datag_all$pcovar);  
  fit.noClimate <- lm(logarea.t1 ~ logarea.t0, data=datag_all); 
  
  if(WEEK==F){
    # make the betas for precip 
    beta_year = rep(0,12); beta_year[2:5] = c(0.5,1,1,0.5); 
    beta_true_p = rep(beta_year,round(nlags/12)+1); 
    beta_true_p =beta_true_p[1:nlags]; 
    beta_true_p[(nEffect+1):nlags]=0; 
    
    # make the betas for temp 
    beta_year = rep(0,12); beta_year[6:9] = -c(0.75,1,1,0.75); 
    beta_true_t = rep(beta_year,round(nlags/12)+1); 
    beta_true_t = beta_true_t[1:nlags]; 
    beta_true_t[(nEffect+1):nlags]=0; 
    
  }else{
    
    beta_year = rep(0,52); 
    beta_year[5:20] = c(seq(0,1, length.out=6),rep(1,4),seq(1,0,length.out=6)); 
    beta_true_p = rep(beta_year,round(nlags/52)+1); 
    beta_true_p =beta_true_p[1:nlags]; 
    beta_true_p[(156*(nEffect/nMonths)+1):nlags]=0; 
    
    # make the betas for temp 
    beta_year = rep(0,52); 
    #beta_year[21:36] = - c(seq(0.75,1,length.out=4),rep(1,8),seq(1,.75,length.out=4)); 
    beta_year[21:36] = -c(seq(0, 1, length.out=6),rep(1,4),seq(1,0, length.out=6)); 
    beta_true_t = rep(beta_year,round(nlags/52)+1); 
    beta_true_t = beta_true_t[1:nlags]; 
    beta_true_t[(156*(nEffect/nMonths)+1):nlags]=0;
    
  }
  
  # make the response 
  response_p = (datag_all$pcovar)%*%beta_true_p*exp(-datag_all$logarea.t0) 
  response_t = (datag_all$tcovar)%*%beta_true_t*exp(-datag_all$logarea.t0)  
  
  # scale betas to get desired target fraction of response variance 
  targetVar = varFrac*var(datag_all$logarea.t1-datag_all$logarea.t0); 
  
  targetVar_p = fPrecip*targetVar;
  targetVar_t = (1-fPrecip)*targetVar;
  beta_true_t = beta_true_t*sqrt(targetVar_t/var(response_t)); 
  beta_true_p = beta_true_p*sqrt(targetVar_p/var(response_p)); 
  
  response_p = (datag_all$pcovar)%*%beta_true_p*exp(-datag_all$logarea.t0)  
  response_t = (datag_all$tcovar)%*%beta_true_t*exp(-datag_all$logarea.t0)  
  response = response_p + response_t;
  
  datag_Myclim <- datag_all; 
  datag_Myclim$logarea.t1 <- predict(fit.noClimate,data=datag_all) + 
    sigmaFrac*sample(fit.noClimate$residuals) + response
  
  p.beta <- beta_true_p; t.beta <- beta_true_t; 
  
  datag<-datag_Myclim 
  datag$Group <- as.factor(datag$Group)
  datag$year<-as.factor(datag$year);
  datag<-datag[,-grep("year.1",names(datag))]
  
  data<-list(datag=datag,p.beta=beta_true_p,t.beta=beta_true_t,response_p=response_p,response_t=response_t)
  
  return(data)
}

