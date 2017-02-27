require(MASS); require(speedglm);  

################################################# 
# Compute the distance weighting function
# B = B-spline basis from eval.basis
# par = matrix of spline coefficients  
#################################################
mk_dist_wts <- function(p,xvals,B1,dplot) {
  # Function for relative competition at desired distances
  # such that the interaction kernel is monotonic decreasing. Uses right-hand
  # rule Reimann sum, with intervals specified by dplot. 
  #
  # Arguments
  #   p: coefficients for the spline basis
  #   xvals: distances at which we want to evaluate our kernel
  #   B1: the basis evaluated at the values in dplot 
  #   dplot: bin boundaries for Riemann sum integration 
  # Value
  #   Relative competition values at the given distances.
  g <- B1%*%p; 
  u=pmax(g,0); 
  gplus2 = u^2/{0.02+u};  #gplus2 = pmax(0,g)^2; 
  h=dplot[2]-dplot[1]; 
  cdf=cumsum(gplus2)*h; cdf=cdf-cdf[1]; 
  cdffun=approxfun(dplot,cdf); 
  return(exp(-cdffun(xvals))); 
}	

####################################################
# Use finite difference to compute gradient of W 
####################################################
# Args
#   p: coefficients for the spline basis
#   xvals: distances at which we want to evaluate our kernel
#   B1: the basis evaluated at the values in dplot 
# Returns
#   matrix dW with length(xvals) rows, ncol(B1) columns
#   j^th column is derivative of W(xvals) wrt p[j]
mk_grad_W <- function(p,xvals,B1,dplot,eps=0.01) { 
  dW <- matrix(0,length(xvals),ncol(B1)) 
  for(j in 1:ncol(B1)) {
    p1 <- p; p1[j]=p1[j]+eps;
    Wup <- mk_dist_wts(p1,xvals,B1,dplot);
    p1[j] <- p1[j]-2*eps; 
    Wdown <- mk_dist_wts(p1,xvals,B1,dplot);
    dW[,j] <- (Wup-Wdown)/(2*eps);
  }
  return(dW);
}

#######################################################################
# Spline basis setup 
#######################################################################
mk_basis = function(midRings,norder,breaks,dplotBy) {
  # Function for creating spline basis ranging from the smallest to largest  
  # mid-annulus radius. 
  #
  # Arguments:
  #  midRings: average of outer and inner radius for each annulus
  #  norder: The order of the bspline to be used - one higher than the degree 
  #    of the spline.
  #  breaks: location of the breakpoints for the spline basis
  #  dplotBy: the distance between dplot points
  #
  # Returns:
  #  Rbasis: The basis object for the b-spline.
  #  B1: The basis evaluated at a fine grid of points, dplot
  #  dplot:  Fine grid of points used for creating B1
  Rbasis=create.bspline.basis(rangeval=range(breaks), breaks=breaks,norder=norder)
  dplot=seq(min(midRings),max(midRings),by=dplotBy)
  B1=eval.basis(Rbasis,dplot); 
  return(list(Rbasis=Rbasis,B1=B1,dplot=dplot))
}

###############################
# Penalty matrix construction #
###############################
mk_penalty = function(Rbasis,porder,maxDist=max(midRings),hSize,minDist=min(midRings)){
  #Function for creating penalty matrix. t(b) %*% pen %*% b gives specific penalty, 
  #   where b is the vector of specific spline coefficients and pen is the output of
  #   this function.
  #
  #Args:
  #  Rbasis: basis object for the spline 
  #  porder: Order of penalty to use
  #  maxDist: upper limit of integration.
  #  hSize: size of segments to use for midpoint rule integration
  #  minDist: lower limit of integration.  Default is zero
  #
  #Returns:
  #  Matrix of penalty terms
  mGrid = ceiling((maxDist-minDist)/hSize); 
  hSize=(maxDist-minDist)/mGrid; 
  dplot=seq(minDist+(hSize/2),maxDist-(hSize/2),by=hSize);
  rawPen=eval.basis(Rbasis,dplot,Lfdobj=porder);
  pen=t(rawPen) %*% (rawPen * hSize); 
  return(pen);
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
  out.s=glm(survives ~ Group + logarea + W.s, data=data.frame(sppData.s),family=binomial);
  W.g=sppData.g%*%c(rep(0,nonCompLength.g),dist.wts);
  out.g=lm(logarea.t1 ~ Group + logarea.t0 + W.g, data=data.frame(sppData.g));
  g.vals=(rep(1/mean(out.g$residuals^2),dim(sppData.g)[1]));
  s.vals=1/(out.s$fitted.values*(1-out.s$fitted.values));
  omega=c(g.vals,s.vals);
  return(omega);
}

###########################################
# Approx AIC function                     #
###########################################
approx_AIC=function(par,sppData.g,sppData.s,midRings,B1,P,nonCompLength.s,
                    nonCompLength.g,dplot,lambda,gamma=1){
  #  Function for finding approximate AIC using linearization and effective degrees
  #  of freedom (calculated within this function).
  #
  # Argiments:
  #   par: fitted values for spline coefficients
  #   sppData.g: growth data for plants
  #   sppData.s: survival data for plants
  #   midRings: midRing locations
  #   B1: bspline evaluated at dplot
  #   P: penalty matrix
  #   nonCompLengths.s: number of initial columns in sppData.s that don't contain
  #     information on competitor sizes
  #   nonCompLength.g: same, for sppData.g
  #   dplot: fine grid of points for numerical integration
  #   lambda: smoothing parameter value 
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
  
  # Calculate the Ws and carry out regressions (redundant, but takes essentially no time) 
  dist.wts=mk_dist_wts(par,midRings,B1,dplot); 
  W.s=sppData.s%*%c(rep(0,nonCompLength.s),dist.wts);
  out.s=glm(survives ~ Group + logarea + W.s, data=data.frame(sppData.s),family=binomial);
  W.g=sppData.g%*%c(rep(0,nonCompLength.g),dist.wts);
  out.g=lm(logarea.t1 ~ Group + logarea.t0 + W.g, data=data.frame(sppData.g));
  
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
  
  #create G.g matrix of group indicators, add in all but the first column
  #   (because of intercept term)
  G.g=matrix(0,rnum.g,numGroups.g)
  quadGp.g=sppData.g[,which(colnames(sppData.g)=="Group")]
  for(i.col in 1:ncol(G.g)){G.g[,i.col]=(quadGp.g==i.col)*1}
  X[1:rnum.g,(length(par)+4):(length(par)+3+ncol(G.g[,-1]))]=G.g[,-1]
  #Add in columns of 1, z, hat(W) for survival row from 4.10
  X[(1+rnum.g):(rnum.g+rnum.s),length(par)+4+ncol(G.g[,-1])]=1
  X[(1+rnum.g):(rnum.g+rnum.s),length(par)+5+ncol(G.g[,-1])]=sppData.s[,which(colnames(sppData.s)=="logarea")]
  X[(1+rnum.g):(rnum.g+rnum.s),length(par)+6+ncol(G.g[,-1])]=W.s
  
  #create G.s matrix of group indicators, add in all but the first column
  #   (because of intercept term)
  G.s=matrix(0,rnum.s,numGroups.s)
  quadGp.s=sppData.s[,which(colnames(sppData.s)=="Group")]
  for(i.col in 1:ncol(G.s)){G.s[,i.col]=(quadGp.s==i.col)*1}
  X[(1+rnum.g):(rnum.g+rnum.s),(length(par)+7+ncol(G.g[,-1])):ncol(X)]=G.s[,-1]  
  
  # Count columns that are entirely zero (because of how g is created)
  # We don't remove them, instead use pinv() which is equivalent for df
  e = which(apply(abs(X),2,sum)==0)
  emptyCols=length(e); 
  #We have now completed X ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  omega=mk_omega(par=par,sppData.s=sppData.s,sppData.g=sppData.g,B1=B1,
                 nonCompLength.s=nonCompLength.s,nonCompLength.g=nonCompLength.g,dplot=dplot);
  omega = pmin(omega,50); ## just a kludge for now
  xtomega=scale(t(X),center=FALSE,scale=omega); #calculate t(X)%*%omega (for matrix omega)
  
  ############ create M matrix
  M=matrix(0,ncol(X),ncol(X));
  M[1:nrow(P),1:ncol(P)]=P;
  
  ############# calculate hat matrix (two steps for readability)
  inverse=ginv(xtomega%*%X+lambda*M)
  hat=X%*%inverse%*%xtomega
  ############# calculate approx degrees of freedom and AIC 
  dof=sum(diag(hat))
  AIC=-2*(logLik(out.s)+logLik(out.g))+2*gamma*dof;
  return(list(AIC=AIC,dof=dof,emptyCols=emptyCols))
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
  out.s=speedglm(survives ~ Group + logarea + W.s, data=data.frame(sppData.s),family=binomial(),
            sparse=FALSE,eigendec=FALSE)
  W.g=sppData.g%*%c(rep(0,nonCompLength.g),dist.wts)
  out.g=lm(logarea.t1 ~ Group + logarea.t0 + W.g, data=data.frame(sppData.g))
  pen=t(par)%*%P%*%par; # Derivative penalty 
  val = -(logLik(out.s)+logLik(out.g))+lambda*pen;
  ofit <- tryCatch(get("ofit",.GlobalEnv),error=function(p) 0);
  assign("ofit",ofit+1,.GlobalEnv); 
  if(exists("report")) {if(report) cat(ofit, val, log10(lambda),"\n")} 
  return(val)
} 

##########################################################################
# Function for fitting spline distance kernel, looping over multiple lambda 
# values and computing AIC for each lambda 
##########################################################################
fit_spline_aic<-function(sppData.s,sppData.g,nonCompLength.s,nonCompLength.g,
                         lambdaList,breaks,savName,description,
                         dplotBy=.02,hSize=.01,report=FALSE,
                         midRings=c(seq(1,19,2),seq(22.5,47.5,5),seq(55,145,10))){
  #Function for fitting competition kernel using splines. Takes data, either real or generated,
  # and fits using penalized maximum likelihood
  #
  #Args:
  #  sppData.s: annulus information for survival data.
  #  sppData.g: annulus information for growth data
  #  nonCompLength.s: int for the number of columns of sppData.s that aren't annulus data
  #  nonCompLength.g:int for the number of columns of sppData.g that aren't annulus data
  #  lambdaList: list of lambda values to test.
  #  breaks: vector of locations for knots
  #  savName: name to give file saving fit information.
  #  description: string saying where the data is from, for use in saving results
  #  dplotBy: distance between points for making dplot.  Default is .02
  #  hSize: approximate distance between points for making the penalty matrix.  Default is .01
  #  report: argument for diagnostics. Default is FALSE
  #  midRings: vector of midRing locations
  #Returns:   
  #  Y: matrix where each column is a vector of the fitted kernel for corresponding lambda
  #  Pars: matrix of fitted spline coefficients, each column corresponds to lambdaList
  #  aic: vector of approx aic for each fit, position corresponds to lambdaList
  
  #######################################################################################
  # Initial setup: define midRings, choose breaks, make spline basis and penalty matrix #
  #######################################################################################
  out=mk_basis(midRings,norder=4,breaks=breaks,dplotBy=dplotBy); # cubic splines
  dplot=out$dplot; # dense grid for plotting and numerical integration in penalty and weights. 
  Rbasis=out$Rbasis # spline basis object for creating penalty
  B1=out$B1; # Basis evaluated at dplot, used to evaluate the spline kernel 
  P=mk_penalty(Rbasis,porder=2,maxDist=max(midRings),hSize=hSize); #penalty matrix for 2nd deriv penalty
  
  #######################################################
  # For each listed lambda, fit spline from naive start #
  #######################################################
  Y=matrix(0,length(dplot),length(lambdaList)); #For storing fine grid
  Ymid = matrix(0,length(midRings),length(lambdaList)); # values at midRings 
  Pars=matrix(0,dim(B1)[2],length(lambdaList)); #for storing actual fitted parameters
  midVals=matrix(0,length(midRings),length(lambdaList)); #For storing the estimated values at the midRing points
  aic=dof=emptyCols=rep(0,length(lambdaList));
  for(i.lambda in 1:length(lambdaList)){
    lambda=lambdaList[i.lambda];
    gwd = getwd(); cat(log10(lambda),gwd,"\n");  

	dfit2 <- fit_spline(sppData.s=sppData.s,sppData.g=sppData.g,
		B1=B1,P=P,lambda=lambda,nonCompLength.s=nonCompLength.s,
        nonCompLength.g=nonCompLength.g,dplot=dplot)
    
	Y[,i.lambda]=mk_dist_wts(dfit2$par,dplot,B1,dplot); #Save the fitted kernel for plotting
    Ymid[,i.lambda]=mk_dist_wts(dfit2$par,midRings,B1,dplot); #Save the fitted kernel for plotting
    
    Pars[,i.lambda]=dfit2$par; #save spline coefficients
    out=approx_AIC(par=dfit2$par,sppData.g=sppData.g, #find approximate aic
                   sppData.s=sppData.s,midRings=midRings,B1=B1,
                   P=P,nonCompLength.s=nonCompLength.s,
                   nonCompLength.g=nonCompLength.g,dplot=dplot,lambda=lambda);
    aic[i.lambda]=out$AIC; #save aic
    dof[i.lambda]=out$dof; #save approx degrees of freedom
    emptyCols[i.lambda]=out$emptyCols; #save the number of empty columns (for diagnostics)
    cat(i.lambda,lambda, "  Done", "\n")
  }
  
  return(list(Y=Y, Ymid=Ymid,Pars=Pars,aic=aic,dof=dof,dfit2=dfit2,emptyCols=emptyCols))
}

##########################################################################
# Function to fit the spline distance kernel for one lambda,  
# given a starting guess at coefficient values 
# Optimization alternates between Nelder-Mead (to avoid or escape local
# traps) and NEWUOA to converge faster on local minimum.  
##########################################################################
fit_spline<-function(sppData.s,sppData.g,B1,P,lambda,nonCompLength.s,
                     nonCompLength.g,dplot){
  assign("report",FALSE,.GlobalEnv); 
  gwd=getwd(); cat(log10(lambda),gwd,"\n"); 
  
  
  fopt=10^50; parsave=rep(0,ncol(B1)); 
  for(j in 1:5) {  
     frough=10^50; parrough=rep(0,ncol(B1)); 
     start=0.5*qnorm(seq(1,249,by=1)/250);
     start=matrix(start,ncol(B1),249,byrow=TRUE); 
     for(k in 1:ncol(B1)) start[k,]=sample(start[k,],replace=FALSE);     
    
     for(k in 1:249) {
       fval<- objfun_Ksg(start[,k],sppData.s,sppData.g,B1,P,lambda,nonCompLength.s,nonCompLength.g,dplot) 
       if(fval<frough) {frough=fval; parrough=start[,k]}
       cat(k,fval,frough,"\n"); 
     }
    
     dfit = newuoa(par=parrough,fn=objfun_Ksg,
             sppData.s=sppData.s,sppData.g=sppData.g,B1=B1,P=P,lambda=lambda,
             nonCompLength.s=nonCompLength.s,nonCompLength.g=nonCompLength.g,
             dplot=dplot,control=list(rhobeg=0.5,iprint=25,maxfun=1000,npt=2*ncol(B1)+1)); 
     if(dfit$fval<fopt) {fopt=fval; parsave=dfit$par}
   
   }
             
  cat(log10(lambda),gwd,"\n"); 
  dfit2 = newuoa(par=parsave,fn=objfun_Ksg,
             sppData.s=sppData.s,sppData.g=sppData.g,B1=B1,P=P,lambda=lambda,
             nonCompLength.s=nonCompLength.s,nonCompLength.g=nonCompLength.g,
             dplot=dplot,control=list(rhobeg=0.5,iprint=25,maxfun=2500,npt=2*ncol(B1)+1)); 
 
  
  cat(log10(lambda),gwd,"\n");                  
  dfit3 = bobyqa(par=dfit2$par,fn=objfun_Ksg,
             sppData.s=sppData.s,sppData.g=sppData.g,B1=B1,P=P,lambda=lambda,
             nonCompLength.s=nonCompLength.s,nonCompLength.g=nonCompLength.g,
             dplot=dplot,control=list(rhobeg=0.25,iprint=25,maxfun=7500,npt=2*ncol(B1)+1)); 
             
          
  dfit4 = newuoa(par=dfit3$par,fn=objfun_Ksg,
             sppData.s=sppData.s,sppData.g=sppData.g,B1=B1,P=P,lambda=lambda,
             nonCompLength.s=nonCompLength.s,nonCompLength.g=nonCompLength.g,
             dplot=dplot,control=list(rhobeg=0.1,iprint=25,maxfun=20000,npt=2*ncol(B1)+1));           
                
  return(dfit4)
}

#####################################################
# Objective function for fitting exponential kernel #
#####################################################
objfun_exp_fit=function(par,nonCompLength.s,nonCompLength.g,midRings,sppData.s,sppData.g) {
  alpha=par
  dist.wts=exp(-alpha*midRings) 
  W.s=sppData.s%*%c(rep(0,nonCompLength.s),dist.wts)
  out.s=glm(survives ~ Group + logarea + W.s, data=data.frame(sppData.s),family=binomial)
  W.g=sppData.g%*%c(rep(0,nonCompLength.g),dist.wts)
  out.g=lm(logarea.t1 ~ Group + logarea.t0 + W.g, data=data.frame(sppData.g))
  return(-(logLik(out.s)+logLik(out.g)))
}    

#####################################################
# Objective function for fitting Gaussian kernel #
#####################################################
objfun_gau_fit=function(par,nonCompLength.s,nonCompLength.g,midRings,sppData.s,sppData.g) {
  alpha=par; 
  dist.wts=exp(-(alpha*midRings^2)) 
  W.s=sppData.s%*%c(rep(0,nonCompLength.s),dist.wts)
  out.s=glm(survives ~ Group + logarea + W.s, data=data.frame(sppData.s),family=binomial)
  W.g=sppData.g%*%c(rep(0,nonCompLength.g),dist.wts)
  out.g=lm(logarea.t1 ~ Group + logarea.t0 + W.g, data=data.frame(sppData.g))
  return(-(logLik(out.s)+logLik(out.g)))
}   

peaks<-function(series,span=3) {
  if((span <- as.integer(span)) %% 2 != 1) stop("'span' must be odd")
  z = embed(series, span)
  s=(span-1)/2; 	
  zmid=z[,(1+span)/2]; zother=z[,-(1+span)/2]; zomax=apply(zother,1,max);
  v=(zmid>zomax);
  pad = rep(FALSE, s)
  result = c(pad, v, pad)
  return(result) 
} 

