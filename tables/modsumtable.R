#Isotria
#Turn jags model output into pretty table for isotria manuscript
#By Ailene
#Started May 8, 2018


options(stringsAsFactors = FALSE)

setwd("~/git/isotria/tables")
modsums<-read.csv("isotria4stagemodsum_sept2017.csv", header=TRUE)
row.names(modsums)<-NULL
colnames(modsums)[1]<-"param"
head(modsums)
sumtab <- subset(modsums, select=c(param,X50.,X2.5.,X97.5.))
sumtab<-sumtab[-which(substr(sumtab$param,1,4)=="mean"),]
betas<-sumtab[which(substr(sumtab$param,1,4)=="beta"),]
sigmas<-sumtab[which(substr(sumtab$param,1,4)=="sigm"),]
mus<-sumtab[which(substr(sumtab$param,1,2)=="mu"),]
mus$vitrate<-substr(mus$param,4,6)
mus$vitrate[which(substr(mus$vitrate,3,3)=="[")]<-substr(mus$vitrate[which(substr(mus$vitrate,3,3)=="[")],1,2)
betas$vitrate<-substr(betas$param,6,9)
betas$vitrate[which(substr(betas$vitrate,3,3)=="[")]<-substr(betas$vitrate[which(substr(betas$vitrate,3,3)=="[")],1,2)
betas$vitrate[which(substr(betas$vitrate,4,4)=="[")]<-substr(betas$vitrate[which(substr(betas$vitrate,4,4)=="[")],1,3)
betas$vitrate[which(substr(betas$vitrate,1,1)==".")]<-substr(betas$vitrate[which(substr(betas$vitrate,1,1)==".")],2,4)
sigmas$vitrate<-substr(sigmas$param,7,9)
sigmas$vitrate[which(substr(sigmas$vitrate,3,3)=="2")]<-substr(sigmas$vitrate[which(substr(sigmas$vitrate,3,3)=="2")],1,2)

#sumtab<-sumtab[-which(substr(sumtab$param,1,4)=="beta"),]
#sumtab<-sumtab[-which(substr(sumtab$param,1,4)=="sigm"),]
#sumtab<-sumtab[-which(substr(sumtab$param,1,4)=="devi"),]
#sumtab$vitrate<-substr(sumtab$param,1,3)

sumtab2<-as.data.frame(rbind(mus,sigmas,betas))
sumtab2$group<-sub(".*\\[(.*)\\].*", "\\1", sumtab2$param, perl=TRUE) 
sumtab2<-sumtab2[order(sumtab2$group,sumtab2$vitrate, decreasing=TRUE),]
colnames(sumtab2)[2]<-c("median")
sumtab2$median<-round(sumtab2$median, digits=2)
sumtab2$X2.5.<-round(sumtab2$X2.5., digits=2)
sumtab2$X97.5.<-round(sumtab2$X97.5., digits=2)

write.csv(sumtab2,"output/isotria_modsum.csv", row.names = FALSE)
