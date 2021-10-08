################################################################
################################################################
## Successive adjuvancy basic simulation of theorem 4         ##
################################################################
################################################################
## 
## James Liley, April 27
##



################################################################
## Seed, functions, parameters                                ##
################################################################

# Random seed
seed=1
set.seed(seed)


# Logit function
logistic=function(x,a=1) 1/(1+exp(-a*x))
logit=function(x,a=1) (-1/a)*log((1/x) - 1) 

# Parameters
n=100 # Use this many samples to fit rho_e
p=3 # dimension of x

# Test point
x=c(2,3,-1)

# Function g moves x by this much
mxs0=1

# Compute (and plot) for this many epochs
e_max=50

# Set to TRUE to save figure to file rather than drawing in R
save_plot=FALSE

# rho-eq value
rho_eq=0.2

# Also fit n_sd other equivalent rho_e functions in order to estimate expected value later
n_sd=5


################################################################
## Define functions f and g                                   ##
################################################################

# f, parametrised by beta 
set.seed(seed)
beta=abs(rnorm(p)) # WOLOG assume elements of beta are positive
f=function(x,e=0) as.numeric(logistic((x %*% beta)))

## Random valued g: g(r,x) = x - U(-mxs/2,mxs)*sign(r-r_eq)*|r-r_eq|^g_ex. Identical over epochs. 
g=function(r,x,r_eq=rho_eq,mxs=mxs0,e=0,g_ex=0.5) {
  dx=matrix(runif(length(x),-mxs/2,mxs),nrow=nrow(x),ncol=ncol(x))
  sx=sign(r-r_eq)*abs(r-r_eq)^g_ex
  dx=dx*sx
  return(x-dx)
}


#' build a new funcion with a smaller environment
#' @param f input function
#' @param variables names we are allowing to be captured in the closere
#' @return new function with closure restricted to varaibles
#' @example 
#' x=5; f=function() print(x)
#'  f(); x=6; f()
#'  
#' g=restrictEnvironment(f,"x")
#'  g(); x=7; g()
#' @export
restrictEnvironment <- function(f,varList) {
  oldEnv <- environment(f)
  newEnv <- new.env(parent=parent.env(oldEnv))
  for(v in varList) {
    assign(v,get(v,envir=oldEnv),envir=newEnv)
  }
  environment(f) <- newEnv
  f
}


################################################################
## Delta, gamma, and c0                                       ##
################################################################

# Use this many trials to estimate delta
n_est=50

np=100; nt=100; nx=np*nt
xs0=matrix(rnorm(np*length(x)),np,length(x))
xs=c()
for (i in 1:nt) xs=rbind(xs,xs0)

ipred=matrix(0,nx,n_est)
xlhs=matrix(0,nx,n_est)
for (i in 1:n_est) {
  dat=matrix(rnorm(n*length(x)),n,length(x))
  Y=rbinom(n,1,prob=f(dat))
  mod=glm(Y~.,data=data.frame(cbind(dat,Y)),family=binomial(link=logit))
  ipred[,i]=predict(mod,as.data.frame(xs),type="response")
  xlhs[,i]=f(g(ipred[,i],xs))
}
imean=rowMeans(ipred)
lhs0=f(g(imean,xs))
rhs0=rowMeans(xlhs)

# Now mean over g
lhs=rep(0,np); for (i in 1:nt) lhs=lhs+lhs0[((i-1)*np+1):(i*np)]
rhs=rep(0,np); for (i in 1:nt) rhs=rhs+rhs0[((i-1)*np+1):(i*np)]
lhs=lhs/nt; rhs=rhs/nt

delta=quantile(abs(lhs-rhs),0.99) # typically small

gamma=mxs0/4 # from g

c0=0.8 # unbiased






################################################################
## Fit functions rho_e                                        ##
################################################################

# Epoch 0
D0=matrix(rnorm(length(x)*n),n,length(x))
X0=D0
Y0=rbinom(n,1,prob=f(X0))

model0=glm(Y0~.,data=data.frame(cbind(D0,Y0)),family=binomial(link=logit))
rho0=restrictEnvironment(function(x) predict(model0,as.data.frame(x),type="response"),"model0")

for (ii in 1:n_sd) {
  D0=matrix(rnorm(length(x)*n),n,length(x))
  X0=D0
  Y0=rbinom(n,1,prob=f(X0))
  
  model0=glm(Y0~.,data=data.frame(cbind(D0,Y0)),family=binomial(link=logit))
  assign(paste0("rho0_",ii),restrictEnvironment(function(x) predict(model0,as.data.frame(x),type="response"),"model0"))
}

for (i in 1:e_max) {
  
  D_i=matrix(rnorm(length(x)*n),n,length(x)) # Data for epoch i
  
  X=D_i
  for (j in 0:(i-1)) {
    rho_j=get(paste0("rho",j))
    X=g(rho_j(D_i),X)
  }
  
  Y=rbinom(n,1,prob=f(X))
  mod=glm(Y~.,data=data.frame(cbind(D_i,Y)),family=binomial(link=logit))
  rho_i=function(x) predict(mod,as.data.frame(x),type="response")
  assign(paste0("rho",i),restrictEnvironment(rho_i,"mod"))
  
  for (ii in 1:n_sd) {
    D_i=matrix(rnorm(length(x)*n),n,length(x)) # Data for epoch i
    
    X=D_i
    for (j in 0:(i-1)) {
      rho_j=get(paste0("rho",j))
      X=g(rho_j(D_i),X)
    }
    
    Y=rbinom(n,1,prob=f(X))
    mod=glm(Y~.,data=data.frame(cbind(D_i,Y)),family=binomial(link=logit))
    rho_i=function(x) predict(mod,as.data.frame(x),type="response")
    assign(paste0("rho",i,"_",ii),restrictEnvironment(rho_i,"mod"))
  }
  
  
  print(i) # progress indicator
}


################################################################
## Compute values rho_e^o(x)                                  ##
################################################################


x0=t(matrix(x,length(x),n_oracle))
rho_e_o=rep(0,e_max+1)
E_rho_e1=rep(0,e_max)

rho_e_o[1]=mean(f(x0))
E_rho_e1[1]=mean(f(x0))

X=x0
for (i in 1:e_max) {
  rho=get(paste0("rho",i))
  X1=g(rho(x0),X)
  rho_e_o[i+1]=mean(f(X1))
  
  # Expected value
  exr=rep(0,n_sd)
  for (ii in 1:n_sd) {
    exr[ii] = mean(f(g(get(paste0("rho",i,"_",ii))(x0),X)))
  }
  E_rho_e1[i+1]=mean(exr)
  
  X=X1
}

# Draw plot
if (save_plot) pdf("./thm4demo.pdf",width=4,height=4)
plot(0,type="n",xlim=c(0,e_max),ylim=c(0,1),xlab="e",ylab=expression(rho[e]^o))
lines(0:e_max,rho_e_o)
Arrows(0:(e_max-1),rho_e_o[1:e_max],1:e_max,E_rho_e1[2:(e_max+1)],
  col="gray",arr.length=0.1)
abline(h=rho_eq,col="red",lty=2)
abline(h=rho_eq+ delta/(gamma*c0),col="red")
abline(h=rho_eq- delta/(gamma*c0),col="red")
legend("topright",c(expression(rho[e]^o),expression(paste("Exp. ",rho[e]^o)),
  expression(rho[eq]),expression(paste("I"[1],"/","I"[2]))),lty=c(1,1,2,1),
  col=c("black","gray","red","red"),bty="n")
if (save_plot) dev.off()