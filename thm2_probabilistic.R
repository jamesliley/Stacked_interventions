################################################################
################################################################
## Successive adjuvancy basic simulation of theorem 2         ##
################################################################
################################################################
## 
## James Liley, April 27
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
nlimit=40

# Function governing how to compute rho from values z=P(Y_e|X_e(0)=x) (oracle is mean of z)
rho_fit=function(z,e) mean(f(z,e=e)) # effectively unbiased, no variance

# Set to TRUE to save figure to file rather than drawing in R
save_plot=TRUE



################################################################
## Define functions f and g                                   ##
################################################################

# f, parametrised by beta (beta=(0.62,0.18,0.84) if run with seed). Identical over epochs.
beta=abs(rnorm(p)) # WOLOG assume elements of beta are positive
f=function(x,e=0) as.numeric(logistic(x %*% beta))

## Random valued g:  g(r,x) = x - U(-mxs/2,mxs)*sign(r-r_eq)*|r-r_eq|^gamma. Identical over epochs. 
##  We may take Xi as the set of random variables with support in the sphere |X|<A. In Xi
##   g(r_eq,x) = x 
##   E_g {g(r,x)} = x - (mxs/4)*sign(r-r_eq)*|r-r_eq|^gamma
##   so if r>r_eq, g(r,x) < x componentwise, and hence E_g {f(g(r,x))} < E_g {f(x)}, and likewise if r<r_eq E_g {f(g(r,x))} > E_g {f(x)}
##   For any x with delta_x in Xi, E_g {f(g(r,x))} - E_g {f(x)}=o(r-r_eq) if gamma<1
g=function(r,x,r_eq=0.2,mxs=mxs0,e=0,gamma=0.25) {
  dx=matrix(runif(length(x),-mxs/2,mxs),nrow=nrow(x),ncol=ncol(x))
  sx=sign(r-r_eq)*abs(r-r_eq)^gamma
  dx=dx*sx
  return(x-dx)
}


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

xer=x_e(x,nlimit)
rho_seq=xer$rho



################################################################
## Theoretical bounds, Monte Carlo                            ##
################################################################

nxx=10000 # this many random values of x and rho
xvals=matrix(runif(nxx*3,-5,5),nxx,3)
rvals=runif(nxx,0,1)

fxvals=f(xvals)

nmean=1000 # mean over this many observations of f(g(.,.))
fgxvals=0*fxvals; for (i in 1:nmean) fgxvals=fgxvals + f(g(rvals,xvals))
fgxvals=fgxvals/nmean


ubound=max(fgxvals[which(fxvals<0.2 & rvals<0.2)])
lbound=min(fgxvals[which(fxvals>0.2 & rvals>0.2)])


################################################################
## Draw figures                                               ##
################################################################

# Set up
if (save_plot) pdf(paste0("thm2demo.pdf"),width=4,height=4)

#par(mfrow=c(1,2))

# Limits
ylim1=c(min(c(0.1,rho_seq)),max(c(0.6,rho_seq)))

# Value of rho almost converges
plot(0,type="n",xlab="Epoch e",xlim=c(0,length(rho_seq)),
  ylab=expression(paste("P(Y"[e],"|X"[e],"(0) = x)")),
  ylim=ylim1) # converges slowly to rho_eq
abline(h=0.2,col="red")
abline(h=ubound,col="blue",lty=2)
abline(h=lbound,col="blue",lty=2)

lines(rho_seq)
legend("topright",legend=c(expression(rho[eq]),"Limit bounds"),col=c("red","blue"),lty=c(1,2))


if (save_plot) dev.off()

