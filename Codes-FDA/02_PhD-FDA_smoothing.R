################################################################################
# Short Course: "Statistical methods of data science -                         #
#                An introduction to Functional Data Analysis"                  #
# Università degli Studi di Bergamo                                            #
# Alessandra Menafoglio, MOX-Department of Mathematics, Politecnico di Milano  # 
################################################################################

##-------------------------------------#
## From rough data to smooth functions #
##-------------------------------------#

setwd("~/Documents/Politecnico/Didattica/2024 PhD Bergamo/Codes-FDA")


# Partly based on Ramsay, Hooker e Graves, "Functional Data Analysis with R and Matlab", Springer, 2009
# using the R package fda; part of the codes are courtesy of Prof. Laura Sangalli


##########################################################################################
##########################################################################################

# Upload noisy data
noisycurve=read.table("noisycurvebis.txt",header=T)
head(noisycurve)

Xobs0=noisycurve$X0
abscissa=noisycurve$Abscissa
NT=length(abscissa)

par (mfrow=c(1,1))
plot(abscissa,Xobs0,xlab="t",ylab="observed data")

# Upload true data
truecurve=read.table("truecurve.txt",header=T)
head(truecurve)

points(abscissa,truecurve$X0vera,type="l")

# compute the difference quotient
rappincX1=(Xobs0[3:NT]-Xobs0[1:(NT-2)])/(abscissa[3:NT]-abscissa[1:(NT-2)])
rappincX2=((Xobs0[3:NT]-Xobs0[2:(NT-1)])/(abscissa[3:NT]-abscissa[2:(NT-1)])-(Xobs0[2:(NT-1)]-Xobs0[1:(NT-2)])/(abscissa[2:(NT-1)]-abscissa[1:(NT-2)]))*2/(abscissa[3:(NT)]-abscissa[1:(NT-2)])

par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")

######
# REGRESSION SPLINES
######

# Load package fda
library(fda)

# Set parameters
m=5           # spline order 
degree=m-1    # spline degree 

nbasis=9

# Create the basis
basis = create.bspline.basis(rangeval=c(0,1), nbasis=nbasis, norder=m)
# If breaks are not provided, equispaced knots are created
basis
plot(basis)

# Evaluate the basis on the grid of abscissa
basismat   = eval.basis(abscissa, basis)
basismat

# Fit via LS
lsfit(basismat, Xobs0, intercept=FALSE)$coef


par(mfrow=c(1,1))
Xsp0=basismat %*% lsfit(basismat, Xobs0, intercept=FALSE)$coef
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsp0 ,type="l",col="blue",lwd=2)
abline(v=basis$params)

# to obtain the first derivative (argument Lfdobj=1)
basismat1   = eval.basis(abscissa, basis,Lfdobj=1)
Xsp1=basismat1 %*% lsfit(basismat, Xobs0, intercept=FALSE)$coef

# to obtain the second derivative (argument Lfdobj=2)
basismat2   = eval.basis(abscissa, basis,Lfdobj=2)
Xsp2=basismat2 %*% lsfit(basismat, Xobs0, intercept=FALSE)$coef

par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsp0 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xsp1 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(abscissa,Xsp2 ,type="l",col="blue",lwd=2)
plot(basis)


############################
# alternative code 
Xsp = smooth.basis(argvals=abscissa, y=Xobs0, basis)
Xsp0bis  = eval.fd(abscissa, Xsp$fd) #  the curve smoothing the data
Xsp1bis  = eval.fd(abscissa, Xsp$fd, Lfd=1)
Xsp2bis = eval.fd(abscissa, Xsp$fd, Lfd=2)
df  = Xsp$df   #  the degrees of freedom in the smoothing curve  
df
############################

par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsp0bis ,type="l",col="green",lwd=2)
points(abscissa,Xsp0 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xsp1bis ,type="l",col="green",lwd=2)
points(abscissa,Xsp1 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(abscissa,Xsp2bis ,type="l",col="green",lwd=2)
points(abscissa,Xsp2 ,type="l",col="blue",lwd=2)

# alternative way to compute derivatives:
# derivatives of splines are still splines with the same knot vector, the appropriate order
# and coefficients directly computable from the coefficients of the original spline

# approximate pointwise confidence intervals
# As in linear models we can estimate the variance of x(t) as
# sigma^2*diag[phi*(phi'phi)^{-1}(phi)']
S=basismat%*%solve(t(basismat)%*%basismat)%*%t(basismat)
sigmahat=sqrt(sum((Xsp0-Xobs0)^2)/(NT-nbasis))
lb=Xsp0-qnorm(0.975)*sigmahat*sqrt(diag(S))
ub=Xsp0+qnorm(0.975)*sigmahat*sqrt(diag(S))
plot(abscissa,Xsp0,type="l",col="blue",lwd=2,ylab="")
points(abscissa,lb,type="l",col="blue",lty="dashed")
points(abscissa,ub,type="l",col="blue",lty="dashed")
points(abscissa,truecurve$X0vera,type="l")

#


# Oversmoothing

nbasis=7

basisbis = create.bspline.basis(c(0,1), nbasis, m)
basismatbis   = eval.basis(abscissa, basisbis)

Xsp0bis=basismatbis %*% lsfit(basismatbis, Xobs0, intercept=FALSE)$coef

basismat1bis   = eval.basis(abscissa, basisbis,Lfdobj=1)
Xsp1bis=basismat1bis %*% lsfit(basismatbis, Xobs0, intercept=FALSE)$coef

basismat2bis   = eval.basis(abscissa, basisbis,Lfdobj=2)
Xsp2bis=basismat2bis %*% lsfit(basismatbis, Xobs0, intercept=FALSE)$coef



par(mfrow=c(1,1))
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsp0bis ,type="l",col="green",lwd=2)
abline(v=basisbis$params)



par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsp0bis ,type="l",col="green",lwd=2)
points(abscissa,Xsp0 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xsp1bis ,type="l",col="green",lwd=2)
points(abscissa,Xsp1 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(abscissa,Xsp2bis ,type="l",col="green",lwd=2)
points(abscissa,Xsp2 ,type="l",col="blue",lwd=2)
plot(basisbis)



# undersmoothing

nbasis=30

basister = create.bspline.basis(c(0,1), nbasis, m)
basismatter   = eval.basis(abscissa, basister)

Xsp0ter=basismatter %*% lsfit(basismatter, Xobs0, intercept=FALSE)$coef


par(mfrow=c(1,1))
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsp0ter ,type="l",col="red",lwd=2)
abline(v=basister$params)


basismat1ter   = eval.basis(abscissa, basister,Lfdobj=1)
Xsp1ter=basismat1ter %*% lsfit(basismatter, Xobs0, intercept=FALSE)$coef

basismat2ter   = eval.basis(abscissa, basister,Lfdobj=2)
Xsp2ter=basismat2ter %*% lsfit(basismatter, Xobs0, intercept=FALSE)$coef



par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsp0ter ,type="l",col="red",lwd=2)
points(abscissa,Xsp0bis ,type="l",col="green",lwd=2)
points(abscissa,Xsp0 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xsp1ter ,type="l",col="red",lwd=2)
points(abscissa,Xsp1bis ,type="l",col="green",lwd=2)
points(abscissa,Xsp1 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(abscissa,Xsp2ter ,type="l",col="red",lwd=2)
points(abscissa,Xsp2bis ,type="l",col="green",lwd=2)
points(abscissa,Xsp2 ,type="l",col="blue",lwd=2)

# generalized cross-validation
nbasis=6:30
gcv=numeric(length(nbasis))
for (i in 1:length(nbasis))
{basis = create.bspline.basis(c(0,1), nbasis[i], m)
gcv[i] = smooth.basis(abscissa, Xobs0, basis)$gcv
}
plot(nbasis,gcv)
nbasis[which.min(gcv)]


########
# Bias-Variance tradeoff

sigma=0.003
nbasis=9:15
integrationinterval=11:90
bias=rep(NA,len=length(nbasis))
var=rep(NA,len=length(nbasis))
for (j in 1:length(nbasis))
{
basis = create.bspline.basis(c(0,1), nbasis[j], m)
basismat   = eval.basis(abscissa, basis)
S=basismat%*%solve(t(basismat)%*%basismat)%*%t(basismat)
bias[j]=sum((truecurve$X0vera-S%*%truecurve$X0vera)[integrationinterval])
var[j]=(sigma^2)*sum(diag(S[integrationinterval,integrationinterval]))
}
mse=var+bias^2

par(mfrow=c(1,1))
plot(nbasis,bias^2,ylim=c(0,max(mse)),type="l",ylab="",main="Bias-Variance tradeoff")
points(nbasis,var,col="red",type="l")
points(nbasis,mse,col="green",type="l")
legend('topright', c("Bias","Var","MSE"), col=c("black","red","green"), 
       lty=1, cex=.5)

########


############################################################################
############################################################################

######
# SMOOTHING SPLINES
######

breaks=abscissa[((0:50)*2)+1]

basis = create.bspline.basis(breaks, norder=m)
functionalPar = fdPar(fdobj=basis, Lfdobj=3, lambda=1e-8)  
# functional parameter, having arguments: basis, order of the derivative
# to be penalized, smoothing parameter

Xss=smooth.basis(abscissa, Xobs0, functionalPar)

Xss0 = eval.fd(abscissa, Xss$fd, Lfd=0)
Xss1 = eval.fd(abscissa, Xss$fd, Lfd=1)
Xss2 = eval.fd(abscissa, Xss$fd, Lfd=2)

df  = Xss$df   #  the degrees of freedom in the smoothing curve
df
gcv = Xss$gcv  #  the value of the gcv statistic
gcv

par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xss0 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xss1 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(abscissa,Xss2 ,type="l",col="blue",lwd=2)
plot(basis)

functionalParbis = fdPar(fdobj=basis, Lfdobj=3, lambda=1e-5)  

Xssbis=smooth.basis(abscissa, Xobs0, functionalParbis)

Xss0bis = eval.fd(abscissa, Xssbis$fd, Lfd=0)
Xss1bis = eval.fd(abscissa, Xssbis$fd, Lfd=1)
Xss2bis = eval.fd(abscissa, Xssbis$fd, Lfd=2)

dfbis  = Xssbis$df   #  the degrees of freedom in the smoothing curve
dfbis
gcvbis = Xssbis$gcv  #  the value of the gcv statistic
gcvbis


functionalParter = fdPar(fdobj=basis, Lfdobj=3, lambda=1e-12)  

Xsster=smooth.basis(abscissa, Xobs0, functionalParter)

Xss0ter = eval.fd(abscissa, Xsster$fd, Lfd=0)
Xss1ter = eval.fd(abscissa, Xsster$fd, Lfd=1)
Xss2ter = eval.fd(abscissa, Xsster$fd, Lfd=2)

dfter  = Xsster$df   #  the degrees of freedom in the smoothing curve
dfter
gcvter = Xsster$gcv  #  the value of the gcv statistic
gcvter


par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xss0ter ,type="l",col="red",lwd=2)
points(abscissa,Xss0bis ,type="l",col="green",lwd=2)
points(abscissa,Xss0 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xss1ter ,type="l",col="red",lwd=2)
points(abscissa,Xss1bis ,type="l",col="green",lwd=2)
points(abscissa,Xss1 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(abscissa,Xss2ter ,type="l",col="red",lwd=2)
points(abscissa,Xss2bis ,type="l",col="green",lwd=2)
points(abscissa,Xss2 ,type="l",col="blue",lwd=2)

# Recommendation: when choosing the smoothing parameter, look at the
# derivatives vs central differences

# generalized cross-validation
lambda=c(1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12)
gcv=numeric(length(lambda))
for (i in 1:length(lambda))
{functionalPar = fdPar(fdobj=basis, Lfdobj=3, lambda=lambda[i])  
gcv[i] = smooth.basis(abscissa, Xobs0, functionalPar)$gcv
}
plot(log10(lambda),gcv)
lambda[which.min(gcv)]


######
# Local Polynomial Regression
######

library(KernSmooth)

m=5           # order of the polynomial
degree=m-1    # degree of the polynomial


bw=0.05

Xsm0 <- locpoly(abscissa,Xobs0,degree=degree,
                bandwidth=bw,gridsize=length(abscissa), 
                range.x=range(abscissa))
Xsm0 <- Xsm0$y


par(mfrow=c(1,1))
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsm0 ,type="l",col="blue")


Xsm1 <- locpoly(abscissa,Xobs0,drv=1,degree=degree,bandwidth=bw,
                gridsize=length(abscissa), range.x=range(abscissa))
Xsm1 <- Xsm1$y

Xsm2 <- locpoly(abscissa,Xobs0,drv=2,degree=degree,bandwidth=bw,
                gridsize=length(abscissa), range.x=range(abscissa))
Xsm2 <- Xsm2$y

par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsm0 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xsm1 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(abscissa,Xsm2 ,type="l",col="blue",lwd=2)




bw=0.15

Xsm0bis <- locpoly(abscissa,Xobs0,drv=0,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Xsm0bis <- Xsm0bis$y

Xsm1bis <- locpoly(abscissa,Xobs0,drv=1,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Xsm1bis <- Xsm1bis$y

Xsm2bis <- locpoly(abscissa,Xobs0,drv=2,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Xsm2bis <- Xsm2bis$y


par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsm0bis ,type="l",col="green",lwd=2)
points(abscissa,Xsm0 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xsm1bis ,type="l",col="green",lwd=2)
points(abscissa,Xsm1 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(abscissa,Xsm2bis ,type="l",col="green",lwd=2)
points(abscissa,Xsm2 ,type="l",col="blue",lwd=2)

# a too large bandwidth leads to oversmoothing


 
bw=0.015 # a too small bandwidth leads to undersmoothing
Xsm0ter <- locpoly(abscissa,Xobs0,drv=0,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Xsm0ter <- Xsm0ter$y

Xsm1ter <- locpoly(abscissa,Xobs0,drv=1,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Xsm1ter <- Xsm1ter$y

Xsm2ter <- locpoly(abscissa,Xobs0,drv=2,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Xsm2ter <- Xsm2ter$y

par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsm0ter ,type="l",col="red",lwd=2)
points(abscissa,Xsm0bis ,type="l",col="green",lwd=2)
points(abscissa,Xsm0 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xsm1ter ,type="l",col="red",lwd=2)
points(abscissa,Xsm1bis ,type="l",col="green",lwd=2)
points(abscissa,Xsm1 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(abscissa,Xsm2ter ,type="l",col="red",lwd=2)
points(abscissa,Xsm2bis ,type="l",col="green",lwd=2)
points(abscissa,Xsm2 ,type="l",col="blue",lwd=2)

# recommendation: when choosing the bandwidth, look at the derivatives 
# vs central differences

#############################################################################################################

# Constrained functions


######
# SMOOTHING for positive curves
######

# y_j = exp(w(t_j)) + e_j

# f(t) = exp(w(t)) 

# The function w(t) is unconstrained
# The function f(t) is monotone increasing

# w(t) is modeled via a basis expansion 
# numerical methods are used to compute the coefficients of the basis expansion


# smooth.pos




######
# SMOOTHING for monotone curves
######


# Example: Berkeley growth data
 

growth



# If we neglet considering that the curves must be monotone...

age     = growth$age
heightbasis12 = create.bspline.basis(rangeval = c(1,18), nbasis = 12, norder = 6)
basismat   = eval.basis(evalarg = growth$age, basisobj = heightbasis12)
heightmat  = growth$hgtf
heightcoef = lsfit(x = basismat, y = heightmat, intercept=FALSE)$coef

height=basismat %*% lsfit(basismat, heightmat, intercept=FALSE)$coef

basismat1   = eval.basis(evalarg = growth$age, basisobj = heightbasis12,
                         Lfdobj=1)
heightvelocity=basismat1 %*% lsfit(x = basismat, y = heightmat, 
                                   intercept=FALSE)$coef

basismat2   = eval.basis(evalarg = growth$age, basisobj = heightbasis12,Lfdobj=2)
heightacceleration=basismat2 %*% lsfit(x=basismat, y= heightmat, intercept=FALSE)$coef


par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
matplot(age,height,type="l" )
matplot(age,heightvelocity,type="l" )
abline(h=0)
matplot(age,heightacceleration,type="l")


par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
matplot(age,height,type="l" )
matplot(age[-c(1,2,3,31)],heightvelocity[-c(1,2,3,31),],type="l" )
abline(h=0)
matplot(age[-c(1,2,3,31)],heightacceleration[-c(1,2,3,31),],type="l")


# a negative velocity does not make any sense (girls height does not 
# decrease over this age interval)


# A model for monotone curves

# y_j = b_0 + b_1 * f(t_j) + e_j

# f(t) = int_(t_0)^t exp(w(u)) du

# The function w(t) is unconstrained
# The function f(t) is monotone increasing
# b_1>0 for monotone increasing functions
# b_1<0 for monotone decreasing functions
# b_0 is the value of the function at t_0

# w(t) is modeled via a basis expansion 
# numerical methods are used to compute the coefficients of the basis expansion, as well as b_0, b_1


nage    = length(age)
ageRng  = range(age)
nfine   = 101
agefine = seq(ageRng[1], ageRng[2], length=nfine)

# Let's consider only the first 5 girls 

hgtf   = growth$hgtf[,1:5]
ncasef = dim(hgtf)[2]

# We set up an order 6 spline basis with knots at ages of observations

norder = 6
nbasis = nage - 2 + norder 
wbasis = create.bspline.basis(rangeval = ageRng, nbasis = nbasis, 
                              norder = norder, breaks = age)

# We construct the functional parameter with penalization of the third 
# derivative

Lfdobj    = 3          
lambda    = 10^(-0.5)  
cvecf     = matrix(0, nbasis, ncasef) # this is used as initial value 
                                      # for the numerical techniques
Wfd0      = fd(coef = cvecf, basisobj = wbasis)
growfdPar = fdPar(fdobj = Wfd0, Lfdobj = Lfdobj, lambda = lambda)

# We carry out a monotone smoothing 

growthMon = smooth.monotone(argvals = age, y = hgtf, WfdParobj = growfdPar)


Wfd        = growthMon$Wfd
betaf      = growthMon$beta
hgtfhatfd  = growthMon$yhatfd

velocfdUN     = deriv.fd(expr = hgtfhatfd, Lfdobj = 1)
velocmeanfdUN = mean(velocfdUN)

accelfdUN     = deriv.fd(expr = hgtfhatfd, Lfdobj = 2)
accelmeanfdUN = mean(accelfdUN)


par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
plot(hgtfhatfd, xlim=c(1,18), lty=1, lwd=2,
     cex=2, xlab="Age", ylab="Growth (cm)")
plot(velocfdUN, xlim=c(1,18),  lty=1, lwd=2,
     cex=2, xlab="Age", ylab="Velocity (cm/yr)")
plot(accelfdUN, xlim=c(1,18), ylim=c(-4,3), lty=1, lwd=2,
     cex=2, xlab="Age", ylab="Acceleration (cm/yr/yr)")

plot(wbasis)

##################################################################################################################



######
# Extension to multidimensional curves
######


noisycurve3D=read.table("noisycurve3D.txt",header=T)
Xobs0=noisycurve3D$X0
Yobs0=noisycurve3D$Y0
Zobs0=noisycurve3D$Z0
obs0=rbind(Xobs0,Yobs0,Zobs0)
abscissa=noisycurve3D$Abscissa
NT=length(abscissa)


# let us consider the easier case where we have the vector of
# times/locations of the observations

truecurve3D=read.table("truecurve3D.txt",header=T)
Xtrue0=truecurve3D$X0
Ytrue0=truecurve3D$Y0
Ztrue0=truecurve3D$Z0
true0=rbind(Xtrue0,Ytrue0,Ztrue0)
abscissa=noisycurve3D$Abscissa

library(misc3d)
library(rgl)

open3d()
lines3d(t(true0[1,]),t(true0[2,]),t(true0[3,]),xlab="",ylab="",zlab="",size=3,axes=F)
points3d(t(obs0[1,]),t(obs0[2,]),t(obs0[3,]),xlab="",ylab="",zlab="",size=2,axes=F,pch=19,cex=2)
box3d()


# compute the difference quotient
rappincX1=(Xobs0[3:NT]-Xobs0[1:(NT-2)])/(abscissa[3:NT]-abscissa[1:(NT-2)])
rappincY1=(Yobs0[3:NT]-Yobs0[1:(NT-2)])/(abscissa[3:NT]-abscissa[1:(NT-2)])
rappincZ1=(Zobs0[3:NT]-Zobs0[1:(NT-2)])/(abscissa[3:NT]-abscissa[1:(NT-2)])

rappincX2=((Xobs0[3:NT]-Xobs0[2:(NT-1)])/(abscissa[3:NT]-abscissa[2:(NT-1)])-(Xobs0[2:(NT-1)]-Xobs0[1:(NT-2)])/(abscissa[2:(NT-1)]-abscissa[1:(NT-2)]))*2/(abscissa[3:(NT)]-abscissa[1:(NT-2)])
rappincY2=((Yobs0[3:NT]-Yobs0[2:(NT-1)])/(abscissa[3:NT]-abscissa[2:(NT-1)])-(Yobs0[2:(NT-1)]-Yobs0[1:(NT-2)])/(abscissa[2:(NT-1)]-abscissa[1:(NT-2)]))*2/(abscissa[3:(NT)]-abscissa[1:(NT-2)])
rappincZ2=((Zobs0[3:NT]-Zobs0[2:(NT-1)])/(abscissa[3:NT]-abscissa[2:(NT-1)])-(Zobs0[2:(NT-1)]-Zobs0[1:(NT-2)])/(abscissa[2:(NT-1)]-abscissa[1:(NT-2)]))*2/(abscissa[3:(NT)]-abscissa[1:(NT-2)])




par (mfrow=c(3,3),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)

plot(abscissa,obs0[1,],xlab=expression(tilde(s)),ylab="observed data x",cex=0.1,asp=1)
plot(abscissa,obs0[2,],xlab=expression(tilde(s)),ylab="observed data y",cex=0.1,asp=1)
plot(abscissa,obs0[3,],xlab=expression(tilde(s)),ylab="observed data z",cex=0.1,asp=1)
plot(abscissa[2:(NT-1)],rappincX1,xlab=expression(tilde(s)),ylab="first differences x",type="l",asp=1)
plot(abscissa[2:(NT-1)],rappincY1,xlab=expression(tilde(s)),ylab="first differences y",type="l",asp=1)
plot(abscissa[2:(NT-1)],rappincZ1,xlab=expression(tilde(s)),ylab="first differences z",type="l",asp=1)
plot(abscissa[2:(NT-1)],rappincX2,xlab=expression(tilde(s)),ylab="second differences x",type="l")
plot(abscissa[2:(NT-1)],rappincY2,xlab=expression(tilde(s)),ylab="second differences y",type="l")
plot(abscissa[2:(NT-1)],rappincZ2,xlab=expression(tilde(s)),ylab="second differences z",type="l")




bw=0.05

Xsm0 <- locpoly(abscissa,Xobs0,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Xsm0 <- Xsm0$y

Xsm1 <- locpoly(abscissa,Xobs0,drv=1,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Xsm1 <- Xsm1$y

Xsm2 <- locpoly(abscissa,Xobs0,drv=2,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Xsm2 <- Xsm2$y


Ysm0 <- locpoly(abscissa,Yobs0,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Ysm0 <- Ysm0$y

Ysm1 <- locpoly(abscissa,Yobs0,drv=1,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Ysm1 <- Ysm1$y

Ysm2 <- locpoly(abscissa,Yobs0,drv=2,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Ysm2 <- Ysm2$y


Zsm0 <- locpoly(abscissa,Zobs0,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Zsm0 <- Zsm0$y

Zsm1 <- locpoly(abscissa,Zobs0,drv=1,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Zsm1 <- Zsm1$y

Zsm2 <- locpoly(abscissa,Zobs0,drv=2,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Zsm2 <- Zsm2$y



par (mfrow=c(3,3),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)

plot(abscissa,obs0[1,],xlab="s",ylab="x",cex=0.1,asp=1,xlim=c(0,1))
points(abscissa,Xsm0,type="l",col="blue",lwd=2)
plot(abscissa,obs0[2,],xlab="s",ylab="y",cex=0.1,asp=1,xlim=c(0,1))
points(abscissa,Ysm0,type="l",col="blue",lwd=2)
plot(abscissa,obs0[3,],xlab="s",ylab="z",cex=0.1,asp=1,xlim=c(0,1))
points(abscissa,Zsm0,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX1,xlab="s",ylab="x'",type="l",ylim=c(-0.5,0.5),xlim=c(0,1))
points(abscissa,Xsm1,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincY1,xlab="s",ylab="y'",type="l",ylim=c(-0.5,0.5),xlim=c(0,1))
points(abscissa,Ysm1,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincZ1,xlab="s",ylab="z'",type="l",ylim=c(-0.5,0.5),xlim=c(0,1))
points(abscissa,Zsm1,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX2,xlab="s",ylab="x''",type="l")
points(abscissa,Xsm2,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincY2,xlab="s",ylab="y''",type="l")
points(abscissa,Ysm2,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincZ2,xlab="s",ylab="z''",type="l")
points(abscissa,Zsm2,type="l",col="blue",lwd=2)


open3d()
lines3d(t(true0[1,]),t(true0[2,]),t(true0[3,]),xlab="",ylab="",zlab="",size=3,axes=F)
points3d(t(obs0[1,]),t(obs0[2,]),t(obs0[3,]),size=2,pch=19,cex=2)
lines3d(t(Xsm0),t(Ysm0),t(Zsm0),size=3,col="blue")
box3d()


# or we may use another among the techniques seen above
