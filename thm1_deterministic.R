################################################################
################################################################
## Stacked interventions simulation of theorem 1              ##
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

# Function g moves x by this much
mxs0=0.6

# Compute (and plot) for this many epochs
nlimit=12

# Draw circle with magnification
magnify=TRUE

# Function governing how to compute rho from values z=P(Y_e|X_e(0)=x) (oracle is mean of z)
rho_fit=function(z,e) mean(f(z,e=e)) # effectively unbiased, no variance

# Set to TRUE to save figure to file rather than drawing in R
save_plot=TRUE



################################################################
## Define functions f and g                                   ##
################################################################

## f, parametrised by beta (beta=(0.62,0.18,0.84) if run with seed). Identical over epochs.
##  Function f has limits 0,1 as any element of x-> -infinity, infinity. It is continuous in both arguments.
beta=abs(rnorm(p)) # WOLOG assume elements of beta are positive
f=function(x,e=0) as.numeric(logistic(x %*% beta))

## Deterministically valued g:  g(r,x) = x + s*mxs*sign(r-r_eq)*|r-r_eq|^gamma. Identical over epochs. 
g=function(r,x,r_eq=0.2,mxs=mxs0,e=0,gamma=1/3) {
  dx=matrix(mxs,nrow=nrow(x),ncol=ncol(x))
  sx=sign(r-r_eq)*abs(r-r_eq)^gamma
  dx=dx*sx
  return(x-dx)
}

# Plot g
if (save_plot) pdf(paste0("g_determ.pdf"),width=4,height=4)

plot(0,type="n",xlim=c(-0.5,0.5),ylim=c(-1,1),
     xlab=expression(paste(rho,"-",rho[eq])),
     ylab=expression(paste("g(",rho,",",x[i],")-x"[i])))
xseq=seq(-0.5,0.5,length=100)
yseq=0*xseq; 
for (i in 1:length(xseq)) yseq[i]=g(xseq[i]+0.2,matrix(c(1,1,1),1,3))[1]-1
abline(h=0,col="red",lty=2)
abline(v=0,col="red",lty=2)
lines(xseq,yseq)

if (save_plot) dev.off()



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
## Theoretical bounds, empirical                              ##
################################################################

ntx=100 # grid spacing
x_vals=seq(-5,5,length=ntx)

xall=cbind(rep(x_vals,times=ntx^2),rep(x_vals,times=ntx,each=ntx),rep(x_vals,each=ntx^2))
fxall=f(xall)
fgxall=f(g(fxall,xall))

ubound=max(fgxall[which(fxall<0.2)])
lbound=min(fgxall[which(fxall>0.2)])



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
## Draw figures                                               ##
################################################################

# Set up
if (save_plot) pdf(paste0("thm1demo.pdf"),width=4,height=4)
par(mar=c(4.1,4.1,0.1,0.1))

#par(mfrow=c(1,2))

# Limits
ylim1=c(min(c(0.1,rho_seq)),max(c(0.6,rho_seq)))
rho_eq=0.2

# Value of rho almost converges
plot(0,type="n",xlab="Epoch e",xlim=c(0,nlimit-1),
     ylab=expression(paste("P(Y"[e],"|X"[e],"(0) = x)")),
     ylim=ylim1) # converges slowly to rho_eq
abline(h=rho_eq,col="red")
abline(h=ubound,col="blue",lty=2)
abline(h=lbound,col="blue",lty=2)

for (i in 1:dim(rho_seq)[1]) lines(rho_seq[i,],col="gray")
legend("topright",legend=c(expression(rho[eq]),"Limit bounds"),col=c("red","blue"),lty=c(1,2))

if (magnify) {
  ## Draw circle with magnification
  tt=seq(0,2*pi,length=100); xcirc=round(3*nlimit/5); rx=round(nlimit/6); sc=3
  polygon(xcirc + sin(tt)*rx,rho_eq+cos(tt)*rx*(ylim1[2]-ylim1[1])/nlimit,col="white")
  segments(xcirc-rx,rho_eq,xcirc+rx,rho_eq,col="red")
  segments(xcirc-rx,rho_eq + sc*(ubound-rho_eq),xcirc+rx,rho_eq + sc*(ubound-rho_eq),col="blue",lty=2)
  segments(xcirc-rx,rho_eq + sc*(lbound-rho_eq),xcirc+rx,rho_eq + sc*(lbound-rho_eq),col="blue",lty=2)
  xind=(xcirc-rx):(xcirc+rx)
  for (i in 1:dim(rho_seq)[1]) lines(xind,rho_eq + sc*(rho_seq[i,xind]-rho_eq),col="gray")
  
  tn=c(0,tt[1:50])
  polygon(c(xcirc+rx,xcirc+rx,xcirc + cos(tn)*rx,xcirc-rx,xcirc-rx),
          c(rho_eq+0.2,rho_eq,rho_eq+sin(tn)*rx*(ylim1[2]-ylim1[1])/nlimit,rho_eq,rho_eq+0.2),
          col="white",border=NA,fillOddEven=fALSE)
  polygon(c(xcirc+rx,xcirc+rx,xcirc + cos(tn)*rx,xcirc-rx,xcirc-rx),
          c(rho_eq-0.2,rho_eq,rho_eq-sin(tn)*rx*(ylim1[2]-ylim1[1])/nlimit,rho_eq,rho_eq-0.2),
          col="white",border=NA,fillOddEven=fALSE)
  
}

if (save_plot) dev.off()

