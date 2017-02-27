loadSurvivalData=function(dataDir,doSpp) {
######################################################
# Read in the data on Size, Quadrat, and Survival 
######################################################
  nonCompLength.s=5 #Number of columns in SppData that are not measures of competitors 
  survDfile=paste(dataDir,doSpp,"/survD.csv",sep="")
  survD=read.csv(file=survDfile)
  D=survD[survD$allEdge==0,];
  D$year=as.factor(D$year)
  D$logarea=log(D$area)
  D$quad=as.character(D$quad)
  D1=D=D[order(D$X),];
  
##########################################################
# Read in data on neighbors 
##########################################################
  ringD <- read.csv(paste(dataDir,doSpp,"/",doSpp,"_nbhood_rings.csv",sep=""))
  ringD$year<-as.factor(ringD$year)

  # merge D with ringD (D contains fewer rows)
  D<-merge(D,ringD,by.x=c("quad","year","trackID"),by.y=c("quad","year","genetID"))
  D=D[order(D$X),]
  rm(ringD)
  row.names(D) <- NULL  

## pull out annulus data for the focal species  
  sppCols=which(substr(names(D),1,4)==doSpp); 
  sppData=cbind(D$survives,as.factor(D$quad),D$year,D$logarea,as.factor(D$Group),as.matrix(D[,sppCols]));  
                        ### Change this so sppData includes the response 
  colnames(sppData)<-c("survives","quad","year","logarea","Group", colnames(sppData)[(1+nonCompLength.s):length(colnames(sppData))])
  return(list(D=D,sppData=sppData,nonCompLength.s=nonCompLength.s))
}

loadGrowthData=function(dataDir,doSPP) {
##########################################################
# Read in the data on Size, Quadrat, and Growth 
##########################################################
  nonCompLength.g=5; # Number of columns in SppData that are not measures of competitors
  groDfile=paste(dataDir,doSpp,"/growDnoNA.csv",sep="")
  groD=read.csv(file=groDfile)
  D=groD[groD$allEdge==0,]; 
  D$year=as.factor(D$year)
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
  ringD <- read.csv(paste(dataDir,doSpp,"/",doSpp,"_nbhood_rings.csv",sep=""))
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
  return(list(D=D,sppData=sppData,nonCompLength.g=nonCompLength.g))
}

