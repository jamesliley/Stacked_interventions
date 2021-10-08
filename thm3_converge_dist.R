################################################################
################################################################
## Successive adjuvancy basic simulation of theorem 3         ##
################################################################
################################################################
## 
## James Liley, July 27
##



################################################################
## Seed, functions, parameters                                ##
################################################################

# Random seed
set.seed(1)

# Logistic function
logistic=function(x,a=1) 1/(1+exp(-a*x))

# Converge in distribution or not
#  Set to FALSE to demonstrate failure of X_e(1) to converge in 
#  distribution despite P(Y_e|X_e=x) converging 
converge=FALSE

# Parameters
n=1000 # number of samples
p=3 # dimension of x
if (converge) mxs0=0.1 else mxs0=0.5 # g shifts each element of x by up to this much
if (converge) nlimit=1000 else nlimit=1000 # g shifts each element of x by up to this much

# Test point
x=c(2,3,-1)

# Function governing how to compute rho from values z=P(Y_e|X_e(0)=x) (oracle is mean of z)
rho_fit=function(z,e) mean(f(z,e=e)) # effectively unbiased, no variance

# Set to TRUE to save figure to file rather than drawing in R
save_plot=TRUE


################################################################
## Define functions f and g                                   ##
################################################################

# f, parametrised by beta (beta=(0.62,0.18,0.84) if run with seed)
beta=abs(rnorm(p)) # WOLOG assume elements of beta are positive
f=function(x,e=0) as.numeric(logistic(x %*% beta))

## Random valued g:  
## If 'converge' is FALSE, sets g such that Xe does not converge in distribution
if (converge) 
{
  ## In this case, g(r,x) = x + sU; s=(r-r_eq)*sign(r-r_eq),U~U(-mxs/2,mxs)
  ##  We have E{f(g(r,x))-f(x)}/(r-r_eq) < E{f(x + (r-r_eq)*mxs)-f(x)}/(r-r_eq) < mxs*max(grad(f)) < 2-e
  ##  Also for E{f(x)}>r_eq, we have g(E{f(x)},x)-x = (E{f(x)}-r_eq)*U, which depends on x only through f(x), so zeta=(E{f(x)}-r_eq)*U, Z=0
  ##  Now E(zeta|rho)=E(zeta|E{f(x)}) = (rho-r_eq)*mxs/4 = O(rho-r_eq)
  g=function(r,x,r_eq=0.2,mxs=mxs0,e=0) {
    #dx=matrix(runif(length(x),0,mxs),nrow=nrow(x),ncol=ncol(x))
    dx=matrix(runif(length(x),-mxs/2,mxs),nrow=nrow(x),ncol=ncol(x))
    sx=(r-r_eq)
    dx=dx*sx
    return(x-dx)
  }
} else {
  ##  In this case, however, zeta=(E{f(x)}-r_eq)^(1/4) *U
  ##  so E(zeta|rho)=E(zeta|E{f(x)}) = (rho-r_eq)^(1/4) *mxs/4 != O(rho-r_eq)
  g=function(r,x,r_eq=0.2,mxs=mxs0,e=0) {
    dx=matrix(runif(length(x),-mxs/2,mxs),nrow=nrow(x),ncol=ncol(x))
    sx=(r-r_eq)
    dx=dx*sign(sx)*abs(sx)^(1/4)
    return(x-dx)
  }
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

# Compute distribution of X_e at various values of e
xx=x_e(x,nlimit)

ntrial=20 # this many plots
for (ii in 1:ntrial) {
  i=min(nlimit-1,ceiling((ii^(2.5))*(nlimit-1)/(ntrial^(2.5))))
  xe=xx$xe[,,i]
  dx=density(xe[,1])
  assign(paste0("dx",ii),dx)
}

rho_seq=xx$rho; 


################################################################
## Draw figures                                               ##
################################################################

# Set up
if (save_plot) pdf(paste0("thm3demo_conv",converge,".pdf"),width=4,height=4)

# With this g, xe converges in distribution if converge=TRUE and does not otherwise
if (converge) {
  xlim0=c(0,2.5)
  ylim0=c(0,14)
  ylim1=c(0.15,1)
} else {
  xlim0=c(-3,2)
  ylim0=c(0,1.5)
  ylim1=c(0.15,0.65)
}
plot(0,xlim=xlim0,ylim=ylim0,type="n",xlab=expression("First element of X"[e]),ylab="Density")
cx=colorRampPalette(c("red","blue"))(ntrial)
abline(v=x[1])
for (ii in 1:ntrial) {
  dx=get(paste0("dx",ii))
  lines(dx,col=cx[ii])
}
legend("topleft",title="Epoch",legend=c("1","...",nlimit),col=cx[c(1,10,20)],lty=1)

if (save_plot) dev.off()

if (save_plot) pdf(paste0("thm3demo_conv",converge,"_rho.pdf"),width=4,height=4)

#... and rho converges if converge=TRUE, and almost converges otherwise.
plot(0,type="n",xlab="Epoch",xlim=c(0,length(rho_seq)),
  ylab=expression(paste("P(Y"[e],"|X"[e],"(0) = x)")),
  ylim=ylim1) # converges slowly to rho_eq
abline(h=0.2,col="red")
lines(rho_seq)
legend("topright",legend=expression(rho[eq]),col="red",lty=1)

if (save_plot) dev.off()

