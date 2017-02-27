#This function runs LASSO based on x (predictors matrix), and y (the response)
#It also needs a vector of years from the orginal dataframe (not x)
#This script returns both coefficients and probed values

FuncLas<-function(data,x,y,SENS=TRUE){
	
   #data=data$datag;x=x;y=y; SENS=TRUE;	
	
# subset data frame to remove year and response as 'covariates' 
tryBIGg = glmnet(as.matrix(x),y, standardize=TRUE, alpha=1)
# plot(tryBIGg, sub="growth")

# Cross validate LASSO, according to year
foldid<-as.numeric(as.factor(data$year))
cv1g=cv.glmnet(as.matrix(x),y,
               foldid=foldid, standardize=TRUE, alpha=1)

cv1g.coefs<-as.matrix(coef(cv1g, s = "lambda.min"))
cv1g.coefs<-data.frame(cv1g.coefs)
cv1g.coefs$names<-row.names(cv1g.coefs)
names(cv1g.coefs)<-c("cv.coefs","names")

pptCoef<-cv1g.coefs$cv.coefs[grep("p\\.",cv1g.coefs$names)] #exclude intercept
tmpCoef<-cv1g.coefs$cv.coefs[grep("t\\.",cv1g.coefs$names)] #exclude intercept

LAlist<-list(pptCoef,tmpCoef)
names(LAlist)<-c("pptCoef","tmpCoef")

if(SENS){
  
  p.orig<-mean(as.vector(predict(cv1g, as.matrix(x), type="response",s=cv1g$lambda.min)))
    
    pDat<-x[,grep("p\\.",names(x))]; epsp<-.05*mean(apply(pDat,2,sd));
    tDat<-x[,grep("t\\.",names(x))]; epst<-.05*mean(apply(tDat,2,sd));
    
    sens=numeric(ncol(x)-1); newDat<-x;
    
    for (ii in 1:length(sens)){
      newDat<-x; 
      eps<-if(ii %in% grep("p\\.",names(newDat))){epsp}else{epst}; 
      newDat[,ii+1]<-newDat[,ii+1]+eps;
      p.up<-as.vector(predict(cv1g, as.matrix(newDat), type="response",s=cv1g$lambda.min))
      sens[ii]=(mean(p.up)-p.orig)/eps;
    }
    
    LArun<-data.frame(cbind(names(x)[-1],sens))
    names(LArun)<-c("cov","sens")
    
    SensP<-as.numeric(as.character(LArun$sens[grep("p\\.",LArun$cov)]))
    SensT<-as.numeric(as.character(LArun$sens[grep("t\\.",LArun$cov)]))
    
    LAlist<-list(pptCoef,tmpCoef,SensP,SensT)
    names(LAlist)<-c("pptCoef","tmpCoef","SensP","SensT")
  }

return(LAlist)}

