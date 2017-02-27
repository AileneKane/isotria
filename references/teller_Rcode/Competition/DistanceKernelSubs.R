
###############################################################
## Kernels that can be fitted to data or used to compute W
###############################################################
linekernel <- function(x,alpha) {
  a=alpha[1]; b=alpha[2]; v=alpha[3]; f=exp(v)/(1+exp(v)); 
  k=f*exp(-a*x) + (1-f)*pmax(0,1-b*x)
  return(pmax(0,k));
}  

expkernel <- function(x,alpha){exp(-alpha*x)}; 
gaukernel <- function(x,alpha){z=pmax(0,x-1); exp(-alpha*z^2)}; 
stepkernel = function(x,alpha) 1/(1+exp(0.5*(x-alpha)))  # sigmoid 

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
  dist.wts=gaukernel(midRings,alpha) 
  W.s=sppData.s%*%c(rep(0,nonCompLength.s),dist.wts)
  out.s=glm(survives ~ Group + logarea + W.s, data=data.frame(sppData.s),family=binomial)
  W.g=sppData.g%*%c(rep(0,nonCompLength.g),dist.wts)
  out.g=lm(logarea.t1 ~ Group + logarea.t0 + W.g, data=data.frame(sppData.g))
  return(-(logLik(out.s)+logLik(out.g)))
}   


########################################################################
## Importance of W.s: %residual deviance explained by M1 vs. M0 (reduced) 
########################################################################
PDE <- function(M0,M1) {
    d0 <- deviance(M0); d1 <- deviance(M1); 
    round((d0-d1)/d0,digits=3);
}

###############################################$$########################
## Importance of W.g: %residual variance explained by M1 vs. M0 (reduced) 
#################################################$$######################
PVE <- function(M0,M1) {
    d0 <- summary(M0)$r.squared; d1 <- summary(M1)$r.squared; 
    round( (d1-d0)/(1-d0),digits=3);
}

###############################################################
# Adjust coefficient of W.s to give a target value of PDE 
###############################################################
findBonusWs <- function(bonus,pde.target,D.s,W.s,nreps=10) {
    D.s.temp <- D.s; 
    fitS=glm(survives ~ W.s + Group + logarea,  data=D.s, family=binomial)    
    u = predict(fitS); # linear predictor
    u = u + bonus*coef(fitS)["W.s"]*(W.s-mean(W.s)); 
    p = exp(u)/(1+exp(u)); 
    
    pde.S.fake=0; 
    for(j in 1:nreps) {
    D.s.temp$survives <- as.numeric(runif(nrow(D.s))<p);
    fitS.fake=glm(survives ~ W.s + Group + logarea, data=D.s.temp, family=binomial)    
    fitS0.fake=glm(survives ~ Group + logarea, data=D.s.temp, family=binomial)   
    pde.S.fake = PDE(fitS0.fake,fitS.fake) + pde.S.fake;
    } 
    return(pde.S.fake/nreps-pde.target)
}

addBonusWs <- function(bonus,D.s,W.s) {
    D.s.temp <- D.s; 
    fitS=glm(survives ~ W.s + Group + logarea,  data=D.s, family=binomial)    
    u = predict(fitS); # linear predictor
    u = u + bonus*coef(fitS)["W.s"]*(W.s-mean(W.s)); 
    p = exp(u)/(1+exp(u)); 
    D.s.temp$survives <- as.numeric(runif(nrow(D.s))<p);
    return(D.s.temp)
}

