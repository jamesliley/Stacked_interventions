################################################################
################################################################
## Stacked interventions simulation of theorem 5              ##
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
logit=function(x,a=1) -(1/a)*log((1/x)-1)

# Parameters
n=1000 # number of samples
p=3 # dimension of x
alpha=0.005 # q(f_e) will differ from q(f) by at most alpha (alpha-boundedness)
rho_eq=0.2

# Test point
x=c(2,3,-1)

# Function g moves x by this much
mxs0=1

# Compute (and plot) for this many epochs
nlimit=25


# Set to TRUE to save figure to file rather than drawing in R
save_plot=TRUE


# Draw a magnification circle
magnify=TRUE

# Function governing how to compute rho from values z=P(Y_e|X_e(0)=x) (oracle is mean of z)
rho_fit=function(z,e) mean(f(z,e=e)) # effectively unbiased, no variance




################################################################
## Define functions f and g                                   ##
################################################################

# f, parametrised by beta (beta=(0.62,0.18,0.84) if run with seed). Varies over epochs.
beta=abs(rnorm(p)) # WOLOG assume elements of beta are positive
epoch_change=c(0,runif(nlimit,-alpha,alpha))
f=function(x,e=0) as.numeric(logistic((x %*% beta) + epoch_change[e+1]))

## Random valued g:  g(r,x) = x - U(-mxs/2,mxs)*sign(r-r_eq)*|r-r_eq|^gamma. Identical over epochs. 
exp0=1/3
g=function(r,x,r_eq=rho_eq,mxs=mxs0,e=0,ex=exp0) {
  dx=matrix(runif(length(x),-mxs/2,mxs),nrow=nrow(x),ncol=ncol(x))
  sx=sign(r-r_eq)*abs(r-r_eq)^ex
  dx=dx*sx
  return(x-dx)
}


# Empirically estimate gamma
rseq=seq(0,1,length=1000)-rho_eq
gamma=quantile((mxs0/4)*(sign(rseq)*abs(rseq)^exp0)/rseq,0.01) # from g; effect of g is at least gamma*(rho-rho_eq)


################################################################
## Functions to compute X_e and (empirically) rho_e           ##
################################################################

# Compute rho and x_e at epoch e with (not-quite) oracle predictor
x_e = function(x,e,nt=1000) {
  
  xe=array(0,dim=c(nt,length(x),e-1))
  rho=rep(0,e-1)
  
  # epoch 0
  xe0=outer(rep(1,nt),x)
  rho0=rho_fit(xe0,0)
  
  xe[,,1]=xe0
  rho[1]=rho0
  
  for (ex in 1:(e-1)) {
    xe1=g(rho0,xe0) # This is X_e
    rho1=rho_fit(xe1,ex) # this is E(Y|X(0)=x) or E(f(Xe))
    rho0=rho1; xe0=xe1
    xe[,,ex]=xe0
    rho[ex]=rho0
  }
  
  return(list(rho=rho,xe=xe))
}



################################################################
## Compute rho and x_e                                        ##
################################################################

nstart=50 # Number of random starting values
rho_seq=c()
for (i in 1:nstart) {
  set.seed(i)
  x=rnorm(3)
  xer=x_e(x,nlimit)
  rho_seq0=xer$rho
  rho_seq=rbind(rho_seq,rho_seq0)
}


################################################################
## Intervals I and I0                                         ##
################################################################

# Function q
Q=function(x) logit(pmin(pmax(x,1e-4),0.999))
Qinv=function(x) logistic(x)

# Interval I_rho
I_rho=c(rho_eq - 2*alpha/gamma,rho_eq + 2*alpha/gamma)

# S_inf and S_sup: to a reasonable approximation, script-X is 
#  the set of Gaussians
ntrial=1000
inf_set=c()
sup_set=c()
for (i in 1:ntrial) {
  X=c(); for (j in 1:p) X=cbind(X,rnorm(1000,mean=runif(1,-5,5),sd=runif(1,0.1,10)))
  s_inf=FALSE; s_sup=FALSE
  if (Q(mean(f(X)))> -alpha + Q(rho_eq - alpha*(1+gamma)/gamma)) s_inf=TRUE
  if (Q(mean(f(X)))<  alpha + Q(rho_eq + alpha*(1+gamma)/gamma)) s_sup=TRUE
  if (s_inf) {
    inf0=Qinv(-alpha + Q(mean(f(g(mean(f(X)),X)))))
    inf_set=c(inf_set,inf0)
  }
  if (s_sup) {
    sup0=Qinv(alpha + Q(mean(f(g(mean(f(X)),X)))))
    sup_set=c(sup_set,sup0)
  }
}
xlim=c(min(inf_set),max(sup_set))


################################################################
## Draw figures                                               ##
################################################################

# Set up
if (save_plot) pdf(paste0("thm5demo.pdf"),width=4,height=4)

par(mar=c(4.1,4.1,1,1))
#par(mfrow=c(1,2))

# Limits
ylim1=c(min(c(0.1,rho_seq)),max(c(0.6,rho_seq)))

# Value of rho almost converges
plot(0,type="n",xlab="Epoch e",xlim=c(0,dim(rho_seq)[2]),
     ylab=expression(paste("P(Y"[e],"|X"[e],"(0) = x)")),
     ylim=ylim1) # converges slowly to rho_eq
abline(h=rho_eq,col="black")
abline(h=I_rho,col="red",lty=2)
abline(h=xlim,col="blue",lty=2)

for (i in 1:dim(rho_seq)[1]) lines(rho_seq[i,],col="gray")

if (magnify) {
  ## Draw circle with magnification
  tt=seq(0,2*pi,length=100); xcirc=round(3*nlimit/5); rx=round(nlimit/6); sc=3
  polygon(xcirc + sin(tt)*rx,rho_eq+cos(tt)*rx*(ylim1[2]-ylim1[1])/nlimit,col="white")
  segments(xcirc-rx,rho_eq,xcirc+rx,rho_eq,col="black")
  
  segments(xcirc-rx,rho_eq + sc*(I_rho[1]-rho_eq),xcirc+rx,rho_eq + sc*(I_rho[1]-rho_eq),col="red",lty=2)
  segments(xcirc-rx,rho_eq + sc*(I_rho[2]-rho_eq),xcirc+rx,rho_eq + sc*(I_rho[2]-rho_eq),col="red",lty=2)
  segments(xcirc-rx,rho_eq + sc*(xlim[1]-rho_eq),xcirc+rx,rho_eq + sc*(xlim[1]-rho_eq),col="blue",lty=2)
  segments(xcirc-rx,rho_eq + sc*(xlim[2]-rho_eq),xcirc+rx,rho_eq + sc*(xlim[2]-rho_eq),col="blue",lty=2)
  
  xind=(xcirc-rx):(xcirc+rx)
  for (i in 1:dim(rho_seq)[1]) lines(xind,rho_eq + sc*(rho_seq[i,xind]-rho_eq),col="gray")
  
  tn=c(0,tt[1:50]); rw=0.6; eps=0.1
  polygon(c(xcirc+rx+eps,xcirc+rx+eps,xcirc + cos(tn)*rx,xcirc-rx-eps,xcirc-rx-eps),
          c(rho_eq+rw,rho_eq,rho_eq+sin(tn)*rx*(ylim1[2]-ylim1[1])/nlimit,rho_eq,rho_eq+rw),
          col="white",border=NA,fillOddEven=fALSE)
  polygon(c(xcirc+rx+eps,xcirc+rx+eps,xcirc + cos(tn)*rx,xcirc-rx-eps,xcirc-rx-eps),
          c(rho_eq-rw,rho_eq,rho_eq-sin(tn)*rx*(ylim1[2]-ylim1[1])/nlimit,rho_eq,rho_eq-rw),
          col="white",border="white",fillOddEven=fALSE)
  polygon(xcirc + sin(tt)*rx,rho_eq+cos(tt)*rx*(ylim1[2]-ylim1[1])/nlimit,col=NA)
  
}


legend("topright",legend=c(expression(rho[eq]),expression("I"[rho]),"Limit bounds"),
       col=c("black","red","blue"),lty=c(1,2,2))



if (save_plot) dev.off()