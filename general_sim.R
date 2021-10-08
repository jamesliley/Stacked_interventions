################################################################
################################################################
## Successive adjuvancy basic simulation of theorem 2         ##
################################################################
################################################################
## 
## James Liley, April 27
##



################################################################
## Seed, functions, switches, parameters                      ##
################################################################

# Random seed
seed=7364

# Set to TRUE to draw figures
drawfigs=TRUE

# Set to TRUE to save figures
savefigs=FALSE

# Set to TRUE to force a re-run; otherwise, existing objects will be used to draw figures.
rerun=FALSE



################################################################
## Packages, scripts and basic functions                      ##
################################################################

# Packages
library(matrixStats)
library(randomForest)
library(shape)

# Logistic and logit functions
logistic=function(x,a=1) 1/(1+exp(-a*x))
logit=function(x,a=1) (1/a)*log(-1 + 1/x)


################################################################
## Control parameters                                         ##
################################################################


## Settings for main simulation
npop=10000
nmed=10
pmed0=runif(nmed,0,0.5)

## General parameter controlling strength of intervention; higher=heavier intervention.
gamma=0.5

# Eventual number of epochs
e=20 

## Absolute size of coefficients for modifiable and non-modifiable covariates respectively
mod_scale=1
nmod_scale=0.5

# Degree of randomness in intervention. If >3*gamma, then intervention can move covariates in the 'wrong' direction.
int_r0=2

# Method for model fitting; "rf" or "glm"
fmodel="rf"


################################################################
## Simulate covariates of a random population                 ##
################################################################


##' @name newX
##' @description generates a random matrix of size npop containing
##'  medical-like covariates
##' @param seed random seed
##' @param npop size of population to simulate
##' @param pmed frequencies at which medications are prescribed (length is number of separate medications); defaults to 10 'drugs' at random frequencies
##' @return a data frame containing values
newX=function(seed=1,npop=10000,pmed=pmed0) {
  
  # random seed
  set.seed(seed)
  
  ## Covariates 1; non-modifiable. Age, deprivation, medical history
  age=sample(1:100,npop,rep=TRUE,prob=100 - (0.8*(1:100)))
  sex=sample(0:1,npop,rep=TRUE)
  deprivation=sample(1:5,npop,rep=TRUE)
  medhx=sample(1:10,npop,re=TRUE)
  
  ## Covariates 2; modifiable. Medications 1 to 10 (prescribed or not), diet, 
  ##  smoking, alcohol, exercise. Drugs/medications may either increase or 
  ##  decrease risk; with lifestyle factors, higher value = higher risk.
  if (is.null(pmed)) {
    nmed=10
    pmed=runif(nmed,0,0.5) # frequency with which drugs are prescribed
  } else nmed=length(pmed)
  meds=c(); for (i in 1:nmed) meds[[i]]=sample(0:1,npop,rep=T,prob=c(1-pmed[i],pmed[i]))
  meds=as.data.frame(meds); colnames(meds)=paste0("med",1:nmed)
  
  diet=runif(npop,0,10)
  smoking=sample(0:1,npop,rep=T,prob=c(0.7,0.3)) # 30% smoke
  alcohol=runif(npop,0,10)
  
  # Full data matrix
  data=cbind(age,sex,deprivation,medhx,meds,diet,smoking,alcohol)
  
  return(data)
}

################################################################
## Define ground-truth probability of outcome (f)             ##
################################################################

## Uses a linear model, with covariates set semi-randomly
set.seed(3726)

## Values of each covariate corresponding to 'reasonable' overall risk,
##  once age and medhx have been substituted to values corresponding to 
##  subcohort.
xsc=0.9 # Controls the difference between the baseline mean risk and the target risk. Lower=lower target
midvals=c(50,0.5,2,5,xsc*pmed0,xsc*5,xsc*0.5,xsc*5)
names(midvals)=colnames(newX(1,npop=10))

## Set true coefficients
med_coef_gen=rnorm(nmed,sd=0.1*mod_scale)
coef_gen=c(nmod_scale*c(age=0.01,sex=0.5,deprivation=0.2,medhx=0.1),
  mod_scale*c(meds=med_coef_gen,diet=0.1,smoking=0.2,alcohol=0.2))

## Constant
const_gen=-(t(midvals) %*% coef_gen)[1,1]-1


##' @name f
##' @description Ground truth function; probability of outcome given data dat
##' @param dat input matrix with same column names, types, and ranges as above
##' @return vector of probabilities of outcome
f=function(dat,coef=coef_gen,const=const_gen) {
  return(logistic((as.matrix(dat) %*% coef) + const))
}





################################################################
## Split population into subcohorts with various acceptable   ##
##  outcome probabilities                                     ##
################################################################

# indicator function: cohorts based on age and medical history
age_breaks=c(0,40,70,100)
medhx_breaks=c(0,3,7,10)
age_mean=(age_breaks[1:(length(age_breaks)-1)] + age_breaks[2:length(age_breaks)])/2
medhx_mean=(medhx_breaks[1:(length(medhx_breaks)-1)] + medhx_breaks[2:length(medhx_breaks)])/2

# Acceptable values (rho_eq) in each subcohort

##' @name rho_eq_sub
##' @description Generates a list of 'acceptable' values for risk in each age/medical-history subcohort.
##' @param coef the coefficients to use in the estimation of acceptable risk
##' @param const the constant to use in the estimation of acceptable risk
##' @return a vector of acceptable risks with names corresponding to 2d intervals.
rho_eq_sub=function(coef=coef_gen,const=const_gen) {
  rho_eq_val=matrix(0,length(age_breaks),length(medhx_breaks))
  for (i in 1:length(age_mean))
    for (j in 1:length(medhx_mean)) {
      mid0=midvals
      mid0["age"]=age_mean[i]; mid0["medhx"]=medhx_mean[j];
      rho_eq_val[i,j]=logistic(const + (t(mid0)%*%coef)[1,1])
    }
  xa=cut(-1:101,breaks=age_breaks); 
  xb=cut(-1:11,breaks=medhx_breaks);
  rnames=as.vector(outer(levels(xa),levels(xb),function(x,y) paste0(x,"_",y)))
  rvec=as.vector(rho_eq_val[1:(length(age_breaks)-1),1:(length(medhx_breaks)-1)])
  names(rvec)=rnames
  return(rvec)
}

# Values for const=const_gen,coef=coef_gen
rho_eq0=rho_eq_sub()


##' @name subcohorts
##' @description generates subcohort assignments for a matrix of data
##' @param X matrix of data
##' @return subcohort indices, as factor. These can be used to look up acceptable risk values.
subcohorts=function(X) {
  age_f=cut(X$age,breaks=age_breaks)
  medhx_f=cut(X$medhx,breaks=medhx_breaks)
  subc=as.factor(paste0(as.character(age_f),"_",as.character(medhx_f)))
  return(subc)
}




# Draw figure
if (drawfigs) {
  if (savefigs) pdf("./rho_eq_values_by_cohort.pdf",width=4,height=4)
  plot(0,xlim=range(age_breaks),ylim=range(medhx_breaks),type="n",xlab="Age",ylab="Med. Hx.",
    xaxs="i",yaxs="i",xaxt="n",yaxt="n")
  rho_eqm=matrix(rho_eq_sub(),3,3)
  axis(1,at=age_breaks); axis(2,at=medhx_breaks)
  for (i in 1:(length(age_breaks)-1))
    for (j in 1:(length(medhx_breaks)-1)) {
      polygon(age_breaks[c(i,i+1,i+1,i)],medhx_breaks[c(j,j,j+1,j+1)],
        col=gray(1-rho_eqm[i,j]))
      text((age_breaks[i]+age_breaks[i+1])/2,(medhx_breaks[j]+medhx_breaks[j+1])/2,
        round(100*rho_eqm[i,j]))
    }
  if (savefigs) dev.off()
}



################################################################
## Functions to fit risk score and make predictions           ##
################################################################

##' @name fit_risk_score
##' @description Fits a risk score to a design matrix and target, using a simple random forest
##' @param X design matrix
##' @param target target outcome
##' @param ntree number of trees; default 100
##' @param ... passed to randomForest()
##' @return an object which can be used to make predictions with predict()
fit_risk_score=function(X,target,...) {
  if (fmodel=="glm") return(glm(target~.,data=cbind(X,target),family=binomial(link="logit")))
  if (fmodel=="rf") return(randomForest(X,as.factor(target),...))
}

##' @name predict_risk
##' @description Makes a prediction given a model; wrapper for predict function
##' @param X design matrix
##' @param fit fitted risk score object 
##' @param ... passed to predict()
##' @return vector of predictions
predict_risk=function(X,fit,...) {
  if (fmodel=="glm") return(predict(fit,X,type="response"))
  if (fmodel=="rf") return(predict(fit,X,type="prob",...)[,2])
}


## Draw 
if (drawfigs) {
  if (savefigs) pdf("./example_fitted_vs_true.pdf",width=4,height=4)
  ## Generate a random population
  set.seed(seed)
  data=newX(seed,npop,pmed0)
  train=1:(npop/2); test=(npop/2 + 1):npop
  target_prob=f(data)
  target=rbinom(npop,1,prob=target_prob)
  scm=fit_risk_score(data[train,],target[train])
  xprob=predict_risk(data[test,],scm)
  plot(xprob,target_prob[test],xlab="True probability",ylab="Risk score",pch=16,cex=0.5)
  abline(0,1,col="red")
  if (savefigs) dev.off()
}


################################################################
## Define an intervention function (random-valued)            ##
################################################################


##' @name intervention
##' @description Defines a set of interventions on variables X, returning a new matrix X.
##'  The intervention is 'realistic' in the sense of modelling the action of a medical practitioner.
##' @param X input matrix of covariates
##' @param Y risk scores
##' @param rho_eq vector of acceptable risks, attainable from rho_eq_sub()
##' @param med_effect we assume the actor knows the sign of effects of medications on risk.
##' @param int_r parameter controlling degree of randomness in intervention. Inevitably still randomness unless int_r=0; decrease in randomness is not continuous.
##' @return a new matrix similar to X with appropriate changes made.
intervention=function(X,Y,rho_eq=rho_eq0,med_effect=sign(med_coef_gen),int_r=int_r0) {
  
  # General information
  N=dim(X)[1]
  
  # Establish subcohorts of individuals in X
  subc=subcohorts(X)
  
  # Acceptable risks for each individual
  r_eq=rho_eq[as.character(subc)]
  
  # Difference between reported risk and acceptable risk
  delta=Y-r_eq
  
  ### Interventions to lifestyle
  
  # Change to diet (probabilistic). If risk is high, try and improve diet; if risk is low, 
  # do not give specific advice (diet will probabilistically worsen slightly)
  dm=X$diet-(3*gamma*delta); ds=abs(int_r*delta)
  newdiet=runif(N,
    dm - ds,   # Probably change X in the 'right' direction: 
    dm + ds    # Possibly change X in the wrong direction
  )
  newdiet=pmax(pmin(newdiet,10),0)
  
  # Change to alcohol (probabilistic). If risk is high, try and improve alcohol; if risk is low, 
  # do not give specific advice (alcohol will probabilistically worsen slightly)
  am=X$alcohol-(3*gamma*delta); as = abs(int_r*delta)
  newalcohol= runif(N,
    am - as,  # Probably change X in the 'right' direction: 
    am + as   # Possibly change X in the wrong direction
  )
  newalcohol=pmax(pmin(newdiet,10),0)
  
  # Change to smoking probability (probabilistic). 
  # If the person already smokes, then decrease/slightly increase their probability of smoking after
  #  intervention according to their score.
  # Denote d=delta; v0=smoking status pre-intervention; p1(d,v0)=probability of smoking post-intervention
  # p1(0,v0) ≈ 1_(v0=1)               # if at equivocal risk, they stay as they are
  # p1(x,v0) ≈ (1-x) 1_(v0=1) if x>0  # if at higher-than equivocal risk, reduce probability of smoking
  # p1(x,v0) ≈ |x| 1_(v0=0) if x<0    # if at lower-than equivocal risk, probability of smoking passively increases
  v0=X$smoking; x=delta; gx=2*logistic(4*gamma*x)-1
  p1= (1-gx)*(v0==1 & x >= 0) + abs(gx)*(v0==0 & x<0) + 1*(v0==1 & x <= 0)
  if (int_r>0) newsmoke=rbinom(N,1,prob=p1) else newsmoke=p1
  
  # Change to medication probabilities (probabilistic). Similar to smoking.
  medx=setdiff(grep("med",colnames(X),value=T),"medhx")
  nmed=length(medx)
  newmed=X[,medx]
  for (i in 1:length(nmed)) {
    smed=med_effect[i]
    c0=X[,medx[i]]
    gx=2*logistic(4*gamma*x)-1
    if (smed>0) pnew= (1-gx)*(c0==1 & x >= 0) + abs(gx)*(c0==0 & x<0) + 1*(c0==1 & x <= 0)
    if (smed<=0) pnew= (1-abs(gx))*(c0==1 & x <= 0) + abs(gx)*(c0==0 & x>0) + 1*(c0==1 & x > 0)
    if (int_r>0) newmed[,i]=rbinom(N,1,prob=pnew) else newmed[,i]=pnew
  }
  colnames(newmed)=medx
  
  newX=cbind(X$age,X$sex,X$deprivation,X$medhx,newmed,newdiet,newsmoke,newalcohol)
  
  colnames(newX)=colnames(X)
  
  return(newX)    
}

## Draw some examples of interventions
if (drawfigs) {
  if (savefigs) pdf("./intervention1.pdf",width=4,height=4)
  plot(0,type="n",xlim=c(-1,1),ylim=c(0,10),
    xlab=expression(paste(rho," - ",rho[eq])),ylab="Diet/alcohol",yaxt="n") 
  orig=3; xr0=-1; xr1=1; yr0=0; yr1=10
  axis(2,at=c(0,2,3,4,6,8,10),label=c("0","2","Orig.","4","6","8","10"),las=2)
  g1=-(3*gamma - int_r0); g2=-(3*gamma + int_r0)
  polygon(c(-1,1,1,-orig/g2,-1,-1),c(orig - g1,orig + g1,0,0,orig - g2,orig-g1),
    col="gray",border=NA)
  abline(orig,g1,col="red");   abline(orig,g2,col="red")
  abline(v=0,lty=2); abline(h=yr0,lty=2); abline(h=yr1,lty=2)
  abline(h=orig,lwd=1)
  xex=0.6
  Arrows(xex,orig + g1*xex,xex,orig + g2*xex,code=3,arr.type="simple",arr.length=0.2)
  if (savefigs) dev.off()
  
  if (savefigs) pdf("./intervention2.pdf",width=4,height=4)
  plot(0,type="n",xlim=c(-1,1),ylim=c(0,1),
    xlab=expression(paste(rho," - ",rho[eq])),ylab="P(smoking)/P(drug)") 
  xr0=-1; xr1=1; yr0=0; yr1=1; xp=seq(0,1,length=200); xn=-xp
  gp=2*logistic(4*gamma*xp)-1;   gn=2*logistic(4*gamma*xn)-1
  abline(v=0,lty=2); abline(h=yr0,lty=2); abline(h=yr1,lty=2)
  lines(xp,1-gp,col="red"); lines(xn,rep(1,length(xn)),col="red")
  lines(xp,rep(0,length(xp)),col="black"); lines(xn,abs(gn),col="black")
  legend("bottomleft",c("True","False"),lty=c(1,1), col=c("red","black"),title="Orig. value")
  if (savefigs) dev.off()
}



################################################################
## Plotting function for scores                               ##
################################################################

##' @name plotscores
##' @description Draws a plot of risk scores across risk categories.
##' @param X matrix of covariates
##' @param legend draw a legend
plotscores=function(X,legend=FALSE) {
  target_prob=f(X)
  strat=subcohorts(X)
  sval=levels(strat)
  rho_eq0=rho_eq_sub()
  plot(0,xlim=c(0,length(sval)),ylim=c(-0.5,1.1),type="n",xaxt="n",yaxt="n",
    ylab="Risk",xlab="",bty="n")
  dsc=0.05; rl=0.7
  for (i in 1:length(sval)) {
    w=which(strat==sval[i])
    xx=i + rl*((1:length(w))/(1+length(w)))
    dx=density(target_prob[w])
    points(xx-1,target_prob[w],pch=".",col="gray")
    polygon(i + rl/2 + dsc*c(dx$y,-rev(dx$y))-1,c(dx$x,rev(dx$x)),
      col=rgb(0,0,1,alpha=0.5),border=NA)
    lines(c(i,i+rl)-1,rep(rho_eq0[sval[i]],2),col="red",lwd=2)
    points(i + rl/2-1,mean(target_prob[w]),pch=16,cex=0.6)
  }
  sval0=strsplit(levels(strat),"_")
  s1=unlist(lapply(sval0,function(x) x[1]))
  s2=unlist(lapply(sval0,function(x) x[2]))
  #s1=gsub("(","A: ",s1,fixed=T); s1=gsub(",","-",s1); s1=gsub("]","",s1)
  #s2=gsub("(","M: ",s2,fixed=T); s2=gsub(",","-",s2); s2=gsub("]","",s2)
  #sval2=paste0(s1,"; ",s2)
  axis(2,at=seq(0,1,length=6),las=2)
  mline=-0.18; spc=0.15
  mtext("Med.Hx:",side=2,line=3,at=mline,las=2,adj=0)
  mtext("Age:",side=2,line=3,at=mline-spc,las=2,adj=0)
  text((3-(1-rl))/2,mline,s2[1])
  text(3+(3-(1-rl))/2,mline,s2[2])
  text(6+(3-(1-rl))/2,mline,s2[3])
  for (i in 1:3) for (j in 1:3) text(3*i-3 + j-1 + rl/2,mline-spc,s1[3*j -2],srt=90,adj=c(0.8,0.5))
  #axis(1,at=1:length(sval),labels=sval2,las=2)
  segments(c(0,3,6),rep(mline-spc/2,3),c(3,6,9)-(1-rl),rep(mline-spc/2,3),lwd=2)
  
  
  if (legend) 
    legend("topleft",bty="n",c("Scores","Mean","Distribution",expression(rho[eq])),
      col=c("gray","black",rgb(0,0,1,alpha=0.5),"red"),pch=c(16,16,16,NA),lty=c(NA,NA,NA,1))
}





################################################################
## Run simulation                                             ##
################################################################

if (rerun | !exists(paste0("rho",e))) {
  data0=newX(seed,npop,pmed0) # Covariates, pre-intervention
  target0=rbinom(npop,1,prob=f(data0))
  rho0=fit_risk_score(data0,target0)
  print(0)
  for (i in 1:e) {
    data0=newX(seed + i,npop,pmed0) # Covariates, pre-intervention
    riskx=matrix(0,npop,i)
    trueriskx=riskx
    for (j in 1:i) {
      riskx[,j]=predict_risk(data0,get(paste0("rho",j-1)))
    }
    ndata=data0
    for (j in 1:i) {
      ndata=intervention(ndata,riskx[,j])
      trueriskx[,j]=f(ndata)
    }
    targetn=rbinom(npop,1,prob=f(ndata))
    assign(paste0("rho",i),fit_risk_score(data0,targetn))
    print(i)
  }
}

## Draw figures
if (drawfigs) {
  if (!savefigs) par(mfrow=c(1,2))
  if (savefigs) {
    pdf("./scores_pre_intervention.pdf",width=5,height=5)
    par(mar=c(1.1, 4.1, 1.1, 1.1))
  }
  plotscores(data0,legend=TRUE)
  if (savefigs) dev.off()
  if (savefigs) {
    pdf("./scores_post_intervention.pdf",width=5,height=5)
    par(mar=c(1.1, 4.1, 1.1, 1.1))
  }
  plotscores(ndata,legend=TRUE)
  if (savefigs) dev.off()
  if (!savefigs) par(mfrow=c(1,1))
}

