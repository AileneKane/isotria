rm(list=ls(all=TRUE));
graphics.off();

## Working directory must be set here, so the source()'s below run
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/repos"); # modify as needed

workDir=paste(root,"/drivers/Manuscripts/Sheffield/Rcode",sep=""); 
simDatDir=paste(workDir,"/Climate/ArtificialAnalysis/SimData/",sep="");
outDir=paste(workDir,"/Climate/ArtificialAnalysis/Output/",sep="");

ListID<-"MS_cEffort" #"MS_cIndiv";

params<-read.table(paste(simDatDir,ListID,"_params.txt",sep=""),head=T);
params$it<-1:nrow(params);

fits<-read.table(paste(simDatDir,ListID,"_All_lag.txt",sep=""),head=T);

betas<-read.table(paste(simDatDir,ListID,"_betaTrue.txt",sep=""),head=T);
betas$Lag<-betas$lags;

######################################################################
#Figure 1: One set of parameters viewed up close
######################################################################
windows()
colCF<-c("#000000","#FF0000", "#BF3EFF", "#009E73") 
par(mfrow=c(2,2), mar=c(5,5,1,2))

WEEKset=FALSE
lagTy<-if(WEEKset==TRUE){"Weekly"}else{"Monthly"}
dn<-if(WEEKset==TRUE){-.05}else{-.1}
up<-if(WEEKset==TRUE){.2}else{.4}

pfits<-subset(fits, fits$nFyears==50 & fits$varFrac==0.5 & 
                fits$rPT==0 & fits$WEEK==WEEKset)
pbetas<-subset(betas, betas$nFyears==50 & betas$varFrac==0.5 & 
                 betas$rPT==0 & betas$WEEK==WEEKset)

m.beta<-as.vector(with(pbetas,{tapply(pBetaTrue,list(Lag), mean, na.rm=T)}))
plot(m.beta~c(0:(length(m.beta)-1)), col=colCF[1], type="l", lwd=2, cex.lab=2, cex.axis=1.5,
     ylab="Coefs for PPT", xlab=paste(lagTy,"lags",sep=" "), frame.plot=FALSE, 
     ylim=c(min(m.beta)-.05,max(m.beta)+.05))

mtext(paste("(a)",sep=""),side=3,outer=F, adj=0, padj=0)
xW<-if(WEEKset==TRUE){110}else{25}
yW<-if(WEEKset==TRUE){.13}else{.35}
legend(x=xW,y=yW,lty=c(1,1,1,1), lwd=c(2,2,2,2), bty="n", cex=1.5,
       col=c(colCF[1],colCF[2],colCF[3],colCF[4]),
       c("True","FLM","RF","Lasso"))

m.LAc<-as.vector(with(pfits[pfits$Pred=="PPT",],{tapply(LAc,list(Lag), mean, na.rm=T)}))
s.LAc<-sqrt(as.vector(with(pfits[pfits$Pred=="PPT",],{tapply(LAc,list(Lag), var, na.rm=T)}))/5)
lines(m.LAc~c(0:(length(m.beta)-1)), col=colCF[4], type="l", lwd=2)
polygon(c(c(0:(length(m.beta)-1)),rev(c(0:(length(m.beta)-1)))),c(m.LAc-s.LAc,rev(m.LAc+s.LAc)),
        lty = 2, border=colCF[4], col=paste(colCF[4],50, sep=""))

m.FLMc<-as.vector(with(pfits[pfits$Pred=="PPT",],{tapply(FLMc,list(Lag), mean, na.rm=T)}))
s.FLMc<-sqrt(as.vector(with(pfits[pfits$Pred=="PPT",],{tapply(FLMc,list(Lag), var, na.rm=T)}))/5)
lines(m.FLMc~c(0:(length(m.beta)-1)), col=colCF[2], type="l", lwd=2)
polygon(c(c(0:(length(m.beta)-1)),rev(c(0:(length(m.beta)-1)))),c(m.FLMc-s.FLMc,rev(m.FLMc+s.FLMc)),
        lty = 2, border=colCF[2], col=paste(colCF[2],50, sep=""))
lines(m.beta~c(0:(length(m.beta)-1)), col=colCF[1], type="l", lwd=2)
abline(h=0, lty=2)

m.beta<-as.vector(with(pbetas,{tapply(tBetaTrue,list(Lag), mean, na.rm=T)}))
plot(m.beta~c(0:(length(m.beta)-1)), col=colCF[1], type="l", lwd=2, cex.lab=2,cex.axis=1.5,
     ylab="Coefs for TMP", xlab=paste(lagTy,"lags",sep=" "), frame.plot=FALSE, 
     ylim=c(min(m.beta)-.03,max(m.beta)+.001))
mtext(paste("(b)",sep=""),side=3,outer=F, adj=0, padj=0)
m.LAc<-as.vector(with(pfits[pfits$Pred=="TMP",],{tapply(LAc,list(Lag), mean, na.rm=T)}))
s.LAc<-sqrt(as.vector(with(pfits[pfits$Pred=="TMP",],{tapply(LAc,list(Lag), var, na.rm=T)}))/5)
lines(m.LAc~c(0:(length(m.beta)-1)), col=colCF[4], type="l", lwd=2)
polygon(c(c(0:(length(m.beta)-1)),rev(c(0:(length(m.beta)-1)))),c(m.LAc-s.LAc,rev(m.LAc+s.LAc)),
        lty = 2, border=colCF[4], col=paste(colCF[4],50, sep=""))

m.FLMc<-as.vector(with(pfits[pfits$Pred=="TMP",],{tapply(FLMc,list(Lag), mean, na.rm=T)}))
s.FLMc<-sqrt(as.vector(with(pfits[pfits$Pred=="TMP",],{tapply(FLMc,list(Lag), var, na.rm=T)}))/5)
lines(m.FLMc~c(0:(length(m.beta)-1)), col=colCF[2], type="l", lwd=2) 
polygon(c(c(0:(length(m.beta)-1)),rev(c(0:(length(m.beta)-1)))),c(m.FLMc-s.FLMc,rev(m.FLMc+s.FLMc)),
        lty = 2, border=colCF[2], col=paste(colCF[2],50, sep=""))
lines(m.beta~c(0:(length(m.beta)-1)), col=colCF[1], type="l", lwd=2)
abline(h=0, lty=2)

m.RFs<-as.vector(with(pfits[pfits$Pred=="PPT",],{tapply(RFs,list(Lag), mean, na.rm=T)}))
s.RFs<-sqrt(as.vector(with(pfits[pfits$Pred=="PPT",],{tapply(RFs,list(Lag), var, na.rm=T)}))/5)
plot(m.RFs~c(0:(length(m.beta)-1)), col=colCF[3], type="l", lwd=2, ylim=c(-.05,up), cex.lab=2,cex.axis=2,
     ylab="Sens for PPT", xlab=paste(lagTy,"lags",sep=" "), frame.plot=FALSE)
mtext(paste("(c)",sep=""),side=3,outer=F, adj=0, padj=0)
polygon(c(c(0:(length(m.beta)-1)),rev(c(0:(length(m.beta)-1)))),c(m.RFs-s.RFs,rev(m.RFs+s.RFs)),
        lty = 2, border=colCF[3], col=paste(colCF[3],50, sep=""))

m.LAs<-as.vector(with(pfits[pfits$Pred=="PPT",],{tapply(LAs,list(Lag), mean, na.rm=T)}))
s.LAs<-sqrt(as.vector(with(pfits[pfits$Pred=="PPT",],{tapply(LAs,list(Lag), var, na.rm=T)}))/5)
points(m.LAs~c(0:(length(m.beta)-1)), col=colCF[4], type="l", lwd=2)
polygon(c(c(0:(length(m.beta)-1)),rev(c(0:(length(m.beta)-1)))),c(m.LAs-s.LAs,rev(m.LAs+s.LAs)),
        lty = 2, border=colCF[4], col=paste(colCF[4],50, sep=""))

m.FLMs<-as.vector(with(pfits[pfits$Pred=="PPT",],{tapply(FLMs,list(Lag), mean, na.rm=T)}))
s.FLMs<-sqrt(as.vector(with(pfits[pfits$Pred=="PPT",],{tapply(FLMs,list(Lag), var, na.rm=T)}))/5)
points(m.FLMs~c(0:(length(m.beta)-1)),col=colCF[2], type="l",)
polygon(c(c(0:(length(m.beta)-1)),rev(c(0:(length(m.beta)-1)))),c(m.FLMs-s.FLMs,rev(m.FLMs+s.FLMs)),
        lty = 2, border=colCF[2], col=paste(colCF[2],50, sep=""))
abline(h=0, lty=2)

m.RFs<-as.vector(with(pfits[pfits$Pred=="TMP",],{tapply(RFs,list(Lag), mean, na.rm=T)}))
s.RFs<-sqrt(as.vector(with(pfits[pfits$Pred=="TMP",],{tapply(RFs,list(Lag), var, na.rm=T)}))/5)
plot(m.RFs~c(0:(length(m.beta)-1)), col=colCF[3], type="l", ylim=c(dn,.01), lwd=2,cex.lab=2,cex.axis=2,
     ylab="Sens for TMP", xlab=paste(lagTy,"lags",sep=" "), frame.plot=FALSE)
mtext(paste("(d)",sep=""),side=3,outer=F, adj=0, padj=0)
polygon(c(c(0:(length(m.beta)-1)),rev(c(0:(length(m.beta)-1)))),c(m.RFs-s.RFs,rev(m.RFs+s.RFs)),
        lty = 2, border=colCF[3], col=paste(colCF[3],50, sep=""))

m.LAs<-as.vector(with(pfits[pfits$Pred=="TMP",],{tapply(LAs,list(Lag), mean, na.rm=T)}))
s.LAs<-sqrt(as.vector(with(pfits[pfits$Pred=="TMP",],{tapply(LAs,list(Lag), var, na.rm=T)}))/5)
points(m.LAs~c(0:(length(m.beta)-1)), col=colCF[4], type="l", lwd=2)
polygon(c(c(0:(length(m.beta)-1)),rev(c(0:(length(m.beta)-1)))),c(m.LAs-s.LAs,rev(m.LAs+s.LAs)),
        lty = 2, border=colCF[4], col=paste(colCF[4],50, sep=""))

m.FLMs<-as.vector(with(pfits[pfits$Pred=="TMP",],{tapply(FLMs,list(Lag), mean, na.rm=T)}))
s.FLMs<-sqrt(as.vector(with(pfits[pfits$Pred=="TMP",],{tapply(FLMs,list(Lag), var, na.rm=T)}))/5)
points(m.FLMs~c(0:(length(m.beta)-1)), col=colCF[2], type="l")
polygon(c(c(0:(length(m.beta)-1)),rev(c(0:(length(m.beta)-1)))),c(m.FLMs-s.FLMs,rev(m.FLMs+s.FLMs)),
        lty = 2, border=colCF[2], col=paste(colCF[2],50, sep=""))
abline(h=0, lty=2)

#dev.copy2pdf(width=12, height=7, file=paste(outDir,"CL1m.pdf", sep=""))


##########################################################
# Summarize fits of each iteration with one cor.test and one error
# Merge with parameters from original sim
##########################################################
fAb<-merge(fits,betas)
fAb<-fAb[order(fAb$it, fAb$Pred, fAb$Lag),]
rownames(fAb) <- 1:nrow(fAb)

itVals<-unique(fAb$it)
Preds<-unique(fAb$Pred)
RPT<-unique(fAb$rPT)
WEEKset<-c(TRUE,FALSE)
modelty<-c("FLMc","FLMs","RFs","LAc","LAs")
sul<-params[,c("it","rPT","WEEK")]
sul2<-sul[rep(1:nrow(sul),times = length(modelty)),]
sul2$modTy <- rep(modelty, each=nrow(sul))
sl<-sul2[rep(1:nrow(sul2),times = length(Preds)),]
sl$Pred <- rep(Preds, each=nrow(sul2))

newCor <- function(x,y) {ifelse(sd(x)*sd(y)>0,cor.test(x,y)$estimate,0)}
newErr<- function(x,y) {sqrt(sum((x-y)^2))/sqrt(sum((x)^2))}

for (i in 1:nrow(sl)){
  subfAb<-subset(fAb, it==sl[i,"it"] & Pred==sl[i,"Pred"] & rPT==sl[i,"rPT"] & WEEK==sl[i,"WEEK"])
  #choice<-as.character(sl[i,"modTy"])
  
  sl$cor[i]<-if(sl[i,"Pred"]=="PPT"){
      newCor(subfAb[,paste(sl[i,"modTy"])],subfAb$pBetaTrue)
    }else{newCor(subfAb[,paste(sl[i,"modTy"])],subfAb$tBetaTrue)}
  
  sl$err[i]<-if(sl[i,"Pred"]=="PPT"){
    bound<-if(sl[i,"WEEK"]==FALSE){22}else{90}
    newErr(subfAb$pBetaTrue[1:bound],subfAb[1:bound,paste(sl[i,"modTy"])])
  }else{newErr(subfAb$tBetaTrue[1:bound],subfAb[1:bound,paste(sl[i,"modTy"])] )
  }
}

slm<-merge(sl,params[,-c(11,12)], by="it")
plot(cor~nFyears, col=as.factor(modTy), data=slm)

modelty<-c("FLMc","FLMs","RFs","LAc","LAs")

for (i in 1:length(modelty)){
  sub<-subset(slm, slm$modTy==paste(modelty[i]))
  
  OI<-with(sub,as.data.frame(ftable(tapply(cor,list(Pred,varFrac,nFyears,rPT,WEEK), mean, na.rm=T))))
  names(OI)<-c("Pred","varFrac","nFyears","RPT","WEEK","mean.Cor")
  
  var.Cor<-with(sub,as.data.frame(ftable(tapply(cor,list(Pred,varFrac,nFyears,rPT,WEEK), var, na.rm=T))))[,6]
  len.Cor<-with(sub,as.data.frame(ftable(tapply(cor,list(Pred,varFrac,nFyears,rPT,WEEK), length))))[,6]
  OI$se.Cor<-sqrt(var.Cor/len.Cor)
  
  OI$mean.err<-with(sub,as.data.frame(ftable(tapply(err,list(Pred,varFrac,nFyears,rPT, WEEK), mean, na.rm=T))))[,6]
  var.err<-with(sub,as.data.frame(ftable(tapply(err,list(Pred,varFrac,nFyears,rPT,WEEK), var, na.rm=T))))[,6]
  len.err<-with(sub,as.data.frame(ftable(tapply(err,list(Pred,varFrac,nFyears,rPT,WEEK), length))))[,6]
  OI$se.err<-sqrt(var.err/len.err)
  
  OI$modTy<-as.factor(as.character(modelty[i]))
  
  OIL<-if(i==1){OI}else{rbind(OIL,OI)}
}

OIL$nFyears<-as.numeric(as.character(OIL$nFyears))-3
OIL$varFrac<-as.numeric(as.character(OIL$varFrac))
OIL$RPT<-as.numeric(as.character(OIL$RPT))

vF<-unique(OIL$varFrac)
pred<-unique(OIL$Pred)


#############################################################
#Summary Plots (2 and 3)
#############################################################
windows()
colCF<-c("#FF0000", "#BF3EFF", "#009E73") 
rPTsub=0; WEEKset=FALSE

par(mfrow=c(2,2), mar=c(5,5,1,2))
for (j in 1:2){
for (jj in 1:2){    
  s.OIL<-subset(OIL, OIL$Pred==pred[jj] & OIL$varFrac==vF[j] & OIL$RPT==rPTsub & OIL$WEEK==WEEKset &
                  #OIL$modTy %in% c("FLMc","LAc"))
                  OIL$modTy %in% c("FLMs","RFs","LAs"))
  modTy<-unique(s.OIL$modTy)
  s.OIL$nFyears<-s.OIL$nFyears+3
  plot(NA, ylim=c(0,1), xlim=c(10,50), #main=paste(pred[jj],"; vF=", vF[j],"; rPT=", rPTsub, sep=""), 
       cex.lab=1.5, ylab="Correlation", xlab="N Years", frame.plot=FALSE)
  
  lv<-if(j==1){c(1,2)}else{c(3,4)}; lettervec<-c("a","b","c","d")
  mtext(paste("(",lettervec[lv[jj]],")",sep=""),side=3,outer=F, adj=0, padj=0)
  
  for (i in 1:3){
    with(s.OIL[s.OIL$modTy==modTy[i],],{
      polygon(c(nFyears+3,rev(nFyears)),c(mean.Cor-se.Cor,rev(mean.Cor+se.Cor)),
              lty = 2, border=colCF[i])})
  } #, col=paste(colCF[i],50, sep="")
  
  for (i in 1:3){
    points(mean.Cor~nFyears, data=s.OIL[s.OIL$modTy==modTy[i],], 
           col=colCF[i], type="o", pch=i+14, lwd=2, cex=1.5)
  }
}}
  
  legend(x=30, y=.3,lty=c(1,1,1), lwd=c(2,2,2), bty="n", cex=1.25, pch=c(15,16,17),
         col=c(colCF[1],colCF[2],colCF[3]),
         c("FLM","RF","Lasso"))

#dev.copy2pdf(width = 8, height = 8, file=paste(outDir,"CL2m_cEff.pdf", sep=""))


#########################################################################
#RPT plot
#########################################################################
windows()
par(mfrow=c(2,2), mar=c(5,5,1,2))
for (j in 1:2){
  for (jj in 1:2){    
  s.OIL<-subset(OIL, OIL$Pred==pred[jj] & OIL$varFrac==vF[j] 
                  & OIL$nFyears==47 & OIL$RPT %in% c(0,0.5) & #OIL$WEEK==FALSE &
                    #OIL$modTy %in% c("FLMc","LAc"))
                    OIL$modTy %in% c("FLMs","RFs","LAs"))

  s.OIL<-s.OIL[order(s.OIL$modTy,s.OIL$WEEK,s.OIL$RPT),]
  s.OIL$Fac<-seq(1:nrow(s.OIL))
  
  xv<-s.OIL$Fac
  yv<-s.OIL$mean.Cor
  z<-s.OIL$se.Cor
  g=(max(xv)-min(xv))/50
  
  plot(yv~xv, ylim=c(0,1), xaxt="n", ylab="Correlation", cex.lab=1.5, xlab=NA, pch=rep(c(15,16,17),each=4),
       col=rep(c(colCF[1],colCF[2],colCF[3]), each=4), frame.plot=FALSE, cex=1.25)
  axis(side=1, labels=c("FLM", "RF", "Lasso"), at=c(2.5,6.5,10.5), 
       tick=FALSE, line=2, font=2, cex.axis=1.5, cex=1.5)
  axis(side=1, labels=rep(c("M","W"),times=3), 
       at=c(1.5,3.5,5.5,7.5, 9.5,11.5), tick=FALSE, line=0)
  axis(side=1, labels=rep(c("0.0","0.5"),times=6), 
       at=seq(1,12,by=1), tick=FALSE, line=-1)
  
  for(i in 1:length(xv)){
    lines(c(xv[i],xv[i]),c(yv[i]+z[i],yv[i]-z[i]),col=rep(c(colCF[1],colCF[2],colCF[3]), each=4)[i])
    lines(c(xv[i]-g,xv[i]+g),c(yv[i]+z[i],yv[i]+z[i]),col=rep(c(colCF[1],colCF[2],colCF[3]), each=4)[i])
    lines(c(xv[i]-g,xv[i]+g),c(yv[i]-z[i],yv[i]-z[i]),col=rep(c(colCF[1],colCF[2],colCF[3]), each=4)[i])
    }

  abline(v=c(4.5,8.5), lty=2)  
  abline(h=0)
  
  lv<-if(j==1){c(1,2)}else{c(3,4)}; lettervec<-c("a","b","c","d")
  mtext(paste("(",lettervec[lv[jj]],")",sep=""),side=3,outer=F, adj=0, padj=0)
  
}}

#dev.copy2pdf(width = 10, height = 8, file=paste(outDir,"CL3_cEff.pdf", sep=""))


