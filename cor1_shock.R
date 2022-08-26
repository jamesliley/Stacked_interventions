################################################################
################################################################
## Stacked interventions simulation of corollary 1            ##
################################################################
################################################################
## 
## James Liley, April 27 2022
##




################################################################
## Seed, functions, parameters                                ##
################################################################

# Random seed
set.seed(1)

# Logistic function
logistic=function(x,a=1) 1/(1+exp(-a*x))

# Parameters
n=1000 # number of samples
p=3 # dimension of x

# Test point
x=c(2,3,-1)

# Function g moves x by this much
mxs0=1

# Compute (and plot) for this many epochs
nlimit=200

# Number of 'shocks' to system
nchanges=5

# Equivocal risk value
rho_eq=0.2

# Set to TRUE to save figure to file rather than drawing in R
save_plot=FALSE

# Function governing how to compute rho from values z=P(Y_e|X_e(0)=x) (oracle is mean of z)
rho_fit=function(z,e) mean(f(z,e=e)) 


################################################################
## Define functions f and g                                   ##
################################################################

# f, parametrised by beta (beta=(0.62,0.18,0.84) if run with seed). Varies over epochs.
beta=abs(rnorm(p)) # WOLOG assume elements of beta are positive

changepoints=round(seq(1,nlimit+1,length.out=nchanges))

beta_change=matrix(0,nlimit+1,p)
for (j in 1:(nchanges-1)) {
  for (i in 1:p) beta_change[(changepoints[j]:(changepoints[j+1]-1)),i]=runif(1,-2,2)
}

beta_change=t(beta + t(beta_change))

f=function(x,e=0) as.numeric(logistic(x %*% (beta_change[e+1,])))

xabsmin=-5; xabsmax=5
## Random valued g:  g(r,x) = x - U(-mxs/2,mxs)*sign(r-r_eq)*|r-r_eq|^gamma. Identical over epochs. 
g=function(r,x,r_eq=rho_eq,mxs=mxs0,e=0,gamma=1/3) {
  dx=matrix(runif(length(x),-mxs/2,mxs),nrow=nrow(x),ncol=ncol(x))
  sx=sign(r-r_eq)*abs(r-r_eq)^gamma
  dx=t(t(dx)*sx*sign(beta_change[e+1,]))
  out=matrix(pmin(xabsmax,pmax(xabsmin,x-dx)),dim(x-dx))
  return(out)
}


################################################################
## Functions to compute X_e and (empirically) rho_e           ##
################################################################

# Compute rho and x_e at epoch e. Since g changes with e, this is laborious.
x_e = function(x,e,nt=1000,verbose=TRUE) {
  
  xe=array(0,dim=c(nt,length(x),e-1))
  rho=rep(0,e-1)
  
  # epoch 0
  xe0a=outer(rep(1,nt),x)
  rho0=rho_fit(xe0a,0)
  
  xe[,,1]=xe0a
  rho[1]=rho0
  
  for (ex in 1:(e-1)) {
    xe0=xe0a
    rho0=rho_fit(xe0a,0)
    for (ex2 in 1:ex) {
      xe1=g(rho0,xe0,e=ex) # This is X_e
      rho1=rho_fit(xe1,ex2) # this is E(Y|X(0)=x) or E(f(Xe))
      rho0=rho1; xe0=xe1
    }
    xe[,,ex]=xe1
    rho[ex]=rho1
    if (verbose & ((ex %% 10)==0)) print(ex)
  }    
  return(list(rho=rho,xe=xe))
}



################################################################
## Compute rho and x_e                                        ##
################################################################

nstart=50 # Number of random starting values
rho_seq=c(); xt=c()
for (i in 1:nstart) {
  x=rnorm(3); xt=rbind(xt,x)
  xer=x_e(x,nlimit,verbose=F)
  rho_seq0=xer$rho
  rho_seq=rbind(rho_seq,rho_seq0)
  print(i)
}



################################################################
## Bounding interval                                          ##
################################################################

# S_inf and S_sup: to a reasonable approximation, script-X is 
#  the set of Gaussians
ntrial=10000
inf_set=c()
sup_set=c()
for (i in 1:ntrial) {
  X=c(); for (j in 1:p) X=cbind(X,pmin(xabsmax,pmax(xabsmin,rnorm(1000,mean=runif(1,-5,5),sd=runif(1,0.1,2)))))
  s_inf=FALSE; s_sup=FALSE
  if (mean(f(X))> rho_eq) s_inf=TRUE
  if (mean(f(X))< rho_eq) s_sup=TRUE
  if (s_inf) {
    inf0=mean(f(g(mean(f(X)),X)))
    inf_set=c(inf_set,inf0)
    if (inf0<0.05) print(c(colMeans(X),colSds(X)))
  }
  if (s_sup) {
    sup0=mean(f(g(mean(f(X)),X)))
    sup_set=c(sup_set,sup0)
  }
}
xlim=c(min(inf_set),max(sup_set))




################################################################
## Draw figures                                               ##
################################################################

# Set up
if (save_plot) pdf(paste0("cor1demo.pdf"),width=8,height=4)

# Set up plot parameters
par(mar=c(1,4,1,1))
layout(mat=rbind(matrix(1,4,4),matrix(2,2,4)))

# Limits
ylim1=c(min(c(0.1,rho_seq)),max(c(0.6,rho_seq)))

# Value of rho almost converges
plot(0,type="n",xlab="",xaxt="n",xlim=c(0,dim(rho_seq)[2]),
     ylab=expression(paste("P(Y"[e],"|X"[e],"(0) = x)")),
     ylim=ylim1) # converges slowly to rho_eq

for (i in 1:dim(rho_seq)[1]) lines(rho_seq[i,],col="gray")
abline(h=rho_eq,col="black")
abline(h=xlim,lty=2)
abline(v=changepoints,col="red")

legend("topright",legend=c(expression(rho[eq]),"Change point","Limit bounds"),
       col=c("black","red","black"),lty=c(1,NA,2),pch=c(NA,"|",NA),bg="white")


par(mar=c(4,4,1,1))

# Draw beta values
plot(0,type="n",xlab="Epoch e",xlim=c(0,dim(rho_seq)[2]),
     ylab=expression(beta),
     ylim=range(beta_change)) # converges slowly to rho_eq
for (i in 1:p) lines(beta_change[,i],col=i+1)

legend("bottomright",c(expression(beta[1]),expression(beta[2]),expression(beta[3])),
       lty=c(1,1,1),col=2:4,bg="white")

if (save_plot) dev.off()

