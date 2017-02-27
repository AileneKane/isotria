AnalysisFun<-function(it,consList){

  params<-consList[[1]]
  ListID<-consList[[2]]
  strt<-consList[[3]]
  
  root=ifelse(.Platform$OS.type=="windows","c:/repos","~/repos"); # modify as needed
  workDir=paste(root,"/drivers/Manuscripts/Sheffield/Rcode",sep=""); 
  dataDir=paste(workDir,"/Data",sep="");
  utilDir=paste(workDir,"/Climate/ArtificialAnalysis/Utilities/",sep="")
  simDatDir=paste(workDir,"/Climate/ArtificialAnalysis/SimData/",sep="")

  set.seed(it);   
  
  paramList<-list(monthStart=6,nYrs=3,
                  nMonths=params[it,"nMonths"],WEEK=params[it,"WEEK"],
                  nFyears=params[it,"nFyears"],nIndiv=params[it,"nIndiv"],
                  DETREND=FALSE,SCALE=FALSE,
                  noiseSD_logP=params[it,"noiseSD_logP"],noiseSD_T=params[it,"noiseSD_T"],
                  tauHalf.T=params[it,"tauHalf.T"],tauHalf.P=params[it,"tauHalf.P"],
                  rPT=params[it,"rPT"],doSpp=params[it,"doSpp"],
                  nEffect=params[it,"nEffect"],sigmaFrac=params[it,"sigmaFrac"],
                  fPrecip=params[it,"fPrecip"],varFrac=params[it,"varFrac"])
  
  
  #####################################
  #Import artficial data
  #####################################
  source(paste(utilDir,"FetchData.R", sep=""))
  data<-FetchData(paramList=paramList)
  
  #############################################
  ## Linear Routine
  #############################################
  avYrT1<-rowSums(data$datag$tcovar[,1:12])/12
  avYrT2<-rowSums(data$datag$tcovar[,12:24])/12
  avYrT3<-rowSums(data$datag$tcovar[,25:36])/12
  avYrP1<-rowSums(data$datag$pcovar[,1:12])/12
  avYrP2<-rowSums(data$datag$pcovar[,12:24])/12
  avYrP3<-rowSums(data$datag$pcovar[,25:36])/12
  
  LinRes1<-lm(logarea.t1~logarea.t0+avYrT1+avYrT2+avYrT3+avYrP1+avYrP2+avYrP3,data=data$datag)
  
  #############################################
  ## Flm Routine
  #############################################
  source(paste(utilDir,"FuncFLM.R", sep=""))
  nUnits = ifelse(WEEK,params[it,"nMonths"]*4,params[it,"nMonths"]); 
  FLMRes<-FuncFLM(data=data$datag,nUnits,
                              WEEK=params[it,"WEEK"], SENS=T)
  
  #############################################
  ## RF Routine
  #############################################
  y <- data$datag$logarea.t1; 
  evars <- which(substr(names(data$datag),2,2)==".")
  xvars <- c(evars); 
  x <- as.data.frame(data$datag[,xvars]);
  x$la<-data$datag$logarea.t0
  
  source(paste(utilDir,"FuncRF.R", sep=""))
  RFRes<-FuncRF(x=x,y=y)
  
  #############################################
  # Lasso routine
  #############################################
  y <- data$datag$logarea.t1
  evars <- which(substr(names(data$datag),2,2)==".")
  zvars <- which(names(data$datag)=="logarea.t0"); 
  xvars <- c(zvars,evars); 
  x <- as.data.frame(data$datag[,xvars]); 
  
  source(paste(utilDir,"FuncLas.R", sep=""))
  LARes<-FuncLas(data=data$datag,x=x,y=y, SENS=T)
      
  #############################################
  # RMSE 
  #############################################
  predLin<-predict(LinRes1, type="response")
  RMSE.Lin<-sqrt( sum( (data$datag$logarea.t1 - predLin)^2 , na.rm = TRUE ) / nrow(data$datag) )
  predGam<-as.vector(predict(FLMRes$gam.ridge, type="response"))
  RMSE.gam<-sqrt( sum( (data$datag$logarea.t1 - predGam)^2 , na.rm = TRUE ) / nrow(data$datag) )
  predRF<-predict(RFRes$rFfit, type="response")
  RMSE.RF<-sqrt( sum( (data$datag$logarea.t1 - predRF)^2 , na.rm = TRUE ) / nrow(data$datag) )
  predLA<-predict(LARes$cv1g, newx=as.matrix(x),type="response")
  RMSE.LA<-sqrt( sum( (data$datag$logarea.t1 - predLA)^2 , na.rm = TRUE ) / nrow(data$datag) )
  
  
  resRMSE<-data.frame(RMSE.Lin, RMSE.gam, RMSE.LA, RMSE.RF)
  
  if(it==1){write.table(resRMSE,paste(simDatDir,ListID,"_RMSE.txt",sep=""),
                        na="NA",col.names=T, row.names=F)}else{
                          write.table(resRMSE,paste(simDatDir,ListID,"_RMSE.txt",sep=""),append=T,
                                      na="NA",col.names=F, row.names=F)}


  
  #############################################
  # Storage 
  #############################################
  #source(paste(utilDir,"StoreSimRes.R", sep="")) #Storage for plotting
  #StoreSimRes(ListID=ListID,x=x,y=y,data=data,FLMRes=FLMRes,RFRes=RFRes,LARes=LARes,
  #            nMonths=params[it,"nMonths"], storeX=F,storeX_sum=F,storeY=F,it=it,
  #            params=params)
  
  print(list(paste("iteration",it),strt-proc.time()));
  
  }

