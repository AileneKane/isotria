StoreSimRes<-function(ListID,x,y,data,FLMRes,RFRes,LARes,nMonths,storeX=F,storeX_sum,storeY,it, params){
  
  #Directories 
  root=ifelse(.Platform$OS.type=="windows","c:/repos","~/repos"); # modify as needed
  workDir=paste(root,"/drivers/Manuscripts/Sheffield/Rcode",sep=""); 
  dataDir=paste(workDir,"/Data",sep="");
  utilDir=paste(workDir,"/Climate/ArtificialAnalysis/Utilities/",sep="")
  simDatDir=paste(workDir,"/Climate/ArtificialAnalysis/SimData/",sep="")
  
  lags<-seq(0:(length(RFRes$SensP)-1))
  
  x<-x[,grep("t|p",names(x))]
    
  #store predictors (makes a very large file)
  if(storeX==T){
    storX<-merge(x,params[it,])
    storX$it<-it
    if(it==1){write.table(storX,paste(simDatDir,ListID,"_storX.txt",sep=""),
                          na="NA",col.names=T, row.names=F)}else{
                            write.table(storX,paste(simDatDir,ListID,"_storX.txt",sep=""),append=T,
                                        na="NA",col.names=F, row.names=F)}
  }
  
  #store summary of predictors; mean of lag correlations with all other lags
  if(storeX_sum==T){       
    
    cors.t<-cor(as.matrix(x))[,(nMonths+2):ncol(x)]
    mean.cor.t<-apply(cors.t,2,mean, na.rm=T)
    
    cors.p<-cor(as.matrix(x))[,1:(nMonths+1)]
    mean.cor.p<-apply(cors.p,2,mean, na.rm=T)
    
    storX_sum<-c(mean.cor.p,mean.cor.t)
    corty<-rep(c("PPT","TMP"), each=length(mean.cor.p))
    storX_sum<-data.frame(corty,storX_sum)
    storX_sum$lagno<-rep(lags,times=2)
    storX_sum$it<-it
    
    stX<-merge(storX_sum, params[it,])
    
    if(it==1){write.table(stX,paste(simDatDir,ListID,"_storX_sum.txt",sep=""),
                          na="NA",col.names=T, row.names=F)}else{
                            write.table(stX,paste(simDatDir,ListID,"_storX_sum.txt",sep=""),append=T,
                                        na="NA",col.names=F, row.names=F)}
    
    
  }
  
  #store the responses being analyzed
  if(storeY==T){
    storY<-merge(y,params[it,],it)
    if(it==1){write.table(storY,paste(simDatDir,ListID,"_storY.txt",sep=""),
                          na="NA",col.names=T, row.names=F)}else{
                            write.table(storY,paste(simDatDir,ListID,"_storY.txt",sep=""),append=T,
                                        na="NA",col.names=F, row.names=F)}
  }
  
  #store estimated coeffients (c) or probed responses (s)
 
  Ppt_lag<-data.frame(lags,FLMRes$beta_precip, FLMRes$SensP, RFRes$SensP, LARes$SensP, LARes$pptCoef)
  names(Ppt_lag)<-c("Lag","FLMc","FLMs","RFs","LAs","LAc")
  Ppt_lag$Pred<-"PPT"
  Ppt_lag$it<-it
  
  Temp_lag<-data.frame(lags,FLMRes$beta_temp, FLMRes$SensT, RFRes$SensT, LARes$SensT, LARes$tmpCoef)
  names(Temp_lag)<-c("Lag","FLMc","FLMs","RFs","LAs", "LAc")
  Temp_lag$Pred<-"TMP"
  Temp_lag$it<-it
  
  All_lag<-rbind(Ppt_lag,Temp_lag)
  All_lag<-merge(All_lag,params[it,])
  
  if(it==1){write.table(All_lag,paste(simDatDir,ListID,"_All_lag.txt",sep=""),
                        na="NA",col.names=T, row.names=F)}else{
                          write.table(All_lag,paste(simDatDir,ListID,"_All_lag.txt",sep=""),append=T,
                                      na="NA",col.names=F, row.names=F)}
  #store beta_true
  pBetaTrue<-data$p.beta
  tBetaTrue<-data$t.beta
  
  beta_true<-data.frame(lags,pBetaTrue,tBetaTrue)
  beta_true<-merge(beta_true,params[it,])
  beta_true$it<-it
  
  if(it==1){write.table(beta_true,paste(simDatDir,ListID,"_betaTrue.txt",sep=""),
                        na="NA",col.names=T, row.names=F)}else{
                          write.table(beta_true,paste(simDatDir,ListID,"_betaTrue.txt",sep=""),append=T,
                                      na="NA",col.names=F, row.names=F)}
  
}