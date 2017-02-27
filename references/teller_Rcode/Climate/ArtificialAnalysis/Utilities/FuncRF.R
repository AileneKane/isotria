#This function runs RF based on x (predictors matrix), and y (the response)
#This script returns only probed values

FuncRF<-function(x,y){
  
  rFtune=tuneRF(x=x, y=y, mtryStart=3, ntreeTry=250, stepFactor=1.5, improve=0.05,
                trace=TRUE, plot=FALSE, doBest=FALSE)
  jt=which(rFtune[,2]==min(rFtune[,2]))
  mtry=rFtune[jt,1]; 
  rFfit=randomForest(x=x,y=y,mtry=mtry,ntree=500)
  #rFfitted = predict(rFfit,newdata=x); 
  
  p.orig<-mean(predict(rFfit,newdata=x))
  
  pDat<-x[,grep("p\\.",names(x))]; epsp<-.05*mean(apply(pDat,2,sd));
  tDat<-x[,grep("t\\.",names(x))]; epst<-.05*mean(apply(tDat,2,sd));
  
  sens=numeric(ncol(x)-1); newDat<-x;
  
  for (ii in 1:length(sens)){
    newDat<-x; 
    eps<-if(ii %in% grep("p\\.",names(newDat))){epsp}else{epst}; 
    newDat[,ii]<-newDat[,ii]+eps;
    p.up<-predict(rFfit,newdata=newDat)
    sens[ii]=(mean(p.up)-p.orig)/eps;
  }
  
  RFrun<-data.frame(cbind(names(x)[-length(names(x))],sens))
  names(RFrun)<-c("cov","sens")
  RFrun$sens<-as.numeric(as.character(RFrun$sens))
  
  SensP<-RFrun$sens[grep("p\\.",RFrun$cov)]
  SensT<-RFrun$sens[grep("t\\.",RFrun$cov)]
  
  RFlist<-list(SensP,SensT)
  names(RFlist)<-c("SensP","SensT")
  
  return(RFlist)}

