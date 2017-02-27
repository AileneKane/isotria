## First source DistanceSplineSubs.R, then this file.   
## This modifies some of the functions so that the fitted 
## survival and growth models include a year-dependent intercept. 

###########################################
# Function to make model matrix for years
###########################################
mk_Ymat <- function(year) {
	A = model.matrix(~y, data=data.frame(y=factor(year)))
	return(A[,-1])	
}


#######################################
# Omega Construction function         #                                      
#######################################
mk_omega<-function(par,sppData.s,sppData.g,B1,nonCompLength.s,nonCompLength.g,dplot){
  # Function for creating the omega matrix once we have spline coefficient parameters.
  #    For use in approx AIC
  # Args:
  #  par: fitted spline coefficients
  #  sppData.s: data for survival
  #  sppdata.g: data for growth
  #  B1: spline basis evaluated at dplot
  #  nonCompLength.s: the number of columns in sppData.s that aren't showing competition
  #  nonCompLength.g: the number of columns in sppData.g that aren't showing competition
  # Returns
  #  omega: VECTOR OF VALUES FOR DIAGONAL MATRIX where the ith entry is the inverse of the variance of that
  #     observation.  See FittingDistanceKernel.pdf, currently sec 4.5
  dist.wts=mk_dist_wts(par,midRings,B1,dplot); 
  W.s=sppData.s%*%c(rep(0,nonCompLength.s),dist.wts);
  out.s=glm(survives ~ Group + logarea + W.s + factor(year), data=data.frame(sppData.s),family=binomial);
  W.g=sppData.g%*%c(rep(0,nonCompLength.g),dist.wts);
  out.g=lm(logarea.t1 ~ Group + logarea.t0 + W.g + factor(year), data=data.frame(sppData.g));
  g.vals=(rep(1/mean(out.g$residuals^2),dim(sppData.g)[1]));
  s.vals=1/(out.s$fitted.values*(1-out.s$fitted.values));
  omega=c(g.vals,s.vals);
  return(omega);
}

#############################################################
# Objective function for spline kernel: survival and growth #
#############################################################
objfun_Ksg=function(par,sppData.s,sppData.g,B1,P,lambda,nonCompLength.s,nonCompLength.g,dplot) { 
  # Function for evaluating the fit of our spline to the data by evaluating
  # the LM and GLM fits; also encorporates roughness penalty.
  # Args:
  #  par: spline coefficients
  #  sppData.s: data for survival
  #  sppData.g: data for growth
  #  B1: Our spline basis evaluated at values in dplot 
  #  P: Penalty matrix used for calculating roughness penalty
  # Returns:
  #  returns the combined -logLikelihoods and roughness penalty
  dist.wts=mk_dist_wts(par,midRings,B1,dplot); 
  W.s=sppData.s%*%c(rep(0,nonCompLength.s),dist.wts)
  out.s=speedglm(survives ~ Group + logarea + W.s + factor(year), data=data.frame(sppData.s),family=binomial(),
            sparse=FALSE,eigendec=FALSE)
  W.g=sppData.g%*%c(rep(0,nonCompLength.g),dist.wts)
  out.g=lm(logarea.t1 ~ Group + logarea.t0 + W.g + factor(year), data=data.frame(sppData.g))
  pen=t(par)%*%P%*%par; # Derivative penalty 
  val = -(logLik(out.s)+logLik(out.g))+lambda*pen;
  ofit <- tryCatch(get("ofit",.GlobalEnv),error=function(p) 0);
  assign("ofit",ofit+1,.GlobalEnv); 
  if(exists("report")) {if(report) cat(ofit, val, log10(lambda),"\n")} 
  return(val)
} 

###########################################
# Approx AIC function                     #
###########################################

approx_AIC=function(par,sppData.g,sppData.s,midRings,B1,P,nonCompLength.s,
                    nonCompLength.g,dplot,lambda){
  # Function for finding approximate AIC using linearization and effective degrees
  #   of freedom (calculated within this function).
  #
  # Args:
  #   par: fitted values for spline coefficients
  #   sppData.g: growth data for plants
  #   sppData.s: survival data for plants
  #   midRings: midRing locations
  #   B1: bspline evaluated at dplot
  #   P: penalty matrix
  #   nonCompLengths.s: number of initial columns in sppData.s that don't contain
  #     information on competitor sizes
  #   nonCompLength.g: same but for sppData.g
  #   specIdent: string for regular expression used in grep to identify the 
  #      "competition" columns of sppData.  Replace PSSP with the four-letter
  #      identifier for whatever species is being tested.
  # Returns: 
  #   aic:  approximated AIC.
  #   dof: effective degrees of freedom
  #   emptyCols: number of empty columns from X that were removed.
  ##########
  # First let's find dimensions we're working with
  rnum.g=nrow(sppData.g)
  rnum.s=nrow(sppData.s)
  numGroups.g=length(unique(sppData.g[,which(colnames(sppData.g)=="Group")]))
  numGroups.s=length(unique(sppData.s[,which(colnames(sppData.s)=="Group")]))
  # Calculate the Ws and carry out regressions
  dist.wts=mk_dist_wts(par,midRings,B1,dplot); 
  W.s=sppData.s%*%c(rep(0,nonCompLength.s),dist.wts);
  out.s=glm(survives ~ Group + logarea + W.s + factor(year), data=data.frame(sppData.s),family=binomial);
  W.g=sppData.g%*%c(rep(0,nonCompLength.g),dist.wts);
  out.g=lm(logarea.t1 ~ Group + logarea.t0 + W.g + factor(year), data=data.frame(sppData.g));
  
  
  #Create the X matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  X=matrix(0,rnum.s+rnum.g,length(par)+3+numGroups.g-1+3+numGroups.s-1)
  # Prepare to put in the terms of the first column in pdf, eq 4.10
  # Note: due to the way g is calculated, some columns of W.pregrad may be zero
  C.g=sppData.g[,-(1:nonCompLength.g)]
  C.s=sppData.s[,-(1:nonCompLength.s)]
  W.pregrad=mk_grad_W(p=par,xvals=midRings,B1=B1,dplot=dplot)  
  compcoef.g=out.g$coefficients[which(names(out.g$coefficients)=="W.g")] #competition coefficient in out.g
  compcoef.s=out.s$coefficients[which(names(out.s$coefficients)=="W.s")] #competition coefficient in out.s  
  X[1:rnum.g,1:length(par)]=compcoef.g*C.g%*%W.pregrad
  X[(1+rnum.g):(rnum.g+rnum.s),1:length(par)]=compcoef.s*C.s%*%W.pregrad
  #Add in the columns of 1, z, hat(W)  for growth row from 4.10
  X[1:rnum.g,length(par)+1]=rep(1,rnum.g)
  X[1:rnum.g,length(par)+2]=sppData.g[,which(colnames(sppData.g)=="logarea.t0")]
  X[1:rnum.g,length(par)+3]=W.g
  
  #create G.g matrix of group indicators; add in all but the first column
  #   because of intercept term
  G.g=matrix(0,rnum.g,numGroups.g)
  quadGp.g=sppData.g[,which(colnames(sppData.g)=="Group")]
  for(i.col in 1:ncol(G.g)){G.g[,i.col]=(quadGp.g==i.col)*1}
  X[1:rnum.g,(length(par)+4):(length(par)+3+ncol(G.g[,-1]))]=G.g[,-1]
  #Add in columns of 1, z, hat(W) for survival row from 4.10
  X[(1+rnum.g):(rnum.g+rnum.s),length(par)+4+ncol(G.g[,-1])]=1
  X[(1+rnum.g):(rnum.g+rnum.s),length(par)+5+ncol(G.g[,-1])]=sppData.s[,which(colnames(sppData.s)=="logarea")]
  X[(1+rnum.g):(rnum.g+rnum.s),length(par)+6+ncol(G.g[,-1])]=W.s
  
  #create G.s matrix of group indicators, add in all but the first column
  #   because of intercept term
  G.s=matrix(0,rnum.s,numGroups.s)
  quadGp.s=sppData.s[,which(colnames(sppData.s)=="Group")]
  for(i.col in 1:ncol(G.s)){G.s[,i.col]=(quadGp.s==i.col)*1}
  X[(1+rnum.g):(rnum.g+rnum.s),(length(par)+7+ncol(G.g[,-1])):ncol(X)]=G.s[,-1]  
  
  ##### Add columns for year 
  year.s=sppData.s[,which(colnames(sppData.s)=="year")]
  X.s = mk_Ymat(year.s); 
  year.g=sppData.g[,which(colnames(sppData.g)=="year")]
  X.g = mk_Ymat(year.g); 
  
  X.yr = matrix(0,nrow(X.s)+nrow(X.g),ncol(X.s)+ncol(X.g)); 
  X.yr[1:nrow(X.s),1:ncol(X.s)] = X.s; 
  X.yr[-(1:nrow(X.s)),-(1:ncol(X.s))]=X.g; 
  
  X=cbind(X,X.yr);  
  
  # Count columns that are entirely empty (because of how g is created)
  # We don't remove them, instead use pinv() which is equivalent for df
  e = which(apply(abs(X),2,sum)==0)
  emptyCols=length(e); 
  #X=X[,!e]; 
  #We have now completed X ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  omega=mk_omega(par=par,sppData.s=sppData.s,sppData.g=sppData.g,B1=B1,
                 nonCompLength.s=nonCompLength.s,nonCompLength.g=nonCompLength.g,dplot=dplot);
  omega = pmin(omega,50); ## just a kludge for now
  xtomega=scale(t(X),center=FALSE,scale=omega); #calculate t(X)%*%omega (for matrix omega)
  #create M matrix
  M=matrix(0,ncol(X),ncol(X));
  M[1:nrow(P),1:ncol(P)]=P;
  #calculate hat matrix (two steps for readability)
  inverse=ginv(xtomega%*%X+lambda*M)
  hat=X%*%inverse%*%xtomega
  #calculate approx degrees of freedom
  dof=sum(diag(hat))
  #use fit of regression and dof to find AIC
  AIC=-2*(logLik(out.s)+logLik(out.g))+2*dof;
  return(list(AIC=AIC,dof=dof,emptyCols=emptyCols))
}



