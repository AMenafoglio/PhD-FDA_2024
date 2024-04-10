################################################################################
# Short Course: "Statistical methods of data science -                         #
#                An introduction to Functional Data Analysis"                  #
# Università degli Studi di Bergamo                                            #
# Alessandra Menafoglio, MOX-Department of Mathematics, Politecnico di Milano  # 
################################################################################

##-----------------------------------------##
## Functional principal component Analysis ##
##-----------------------------------------##

setwd("~/Documents/Politecnico/Didattica/2024 PhD Bergamo/Codes-FDA")


# Partly based on Ramsay, Hooker e Graves, "Functional Data Analysis with R and Matlab", Springer, 2009
# using the R packge fda; part of the codes are courtesy of Mara Bernardi and Alessia Pini

library(fda)

# First dataset: canadian weather
# daily temperatures recorded in 35 weather station of Canada 
# (data are averages over 35 years - 1960 to 1994)
data_W <- CanadianWeather$dailyAv[,,1]
head(data_W)
matplot(data_W,type='l',main='Canadian temperature',xlab='Day',ylab='Temperature')

# First of all we smooth the data. We choose a Fourier basis
# (periodic). We need to set the dimension of the basis

# Choice 1: we set a high dimensional basis (interpolating)
# Pros: no loss of information
# Cons: possible overfitting 
time <- 1:365
basis.1 <- create.fourier.basis(rangeval=c(0,365),nbasis=365)
data_W.fd.1 <- Data2fd(y = data_W,argvals = time,basisobj = basis.1)
plot.fd(data_W.fd.1)

# Choice 2: reduced dimensionality
# Pros: the data are much smoother and the measurement error is filtered
# Cons: I could have lost important information
basis.2 <- create.fourier.basis(rangeval=c(0,365),nbasis=21)
data_W.fd.2 <- Data2fd(y = data_W,argvals = time,basisobj = basis.2)
plot.fd(data_W.fd.2)

# Choice 3: compromise between 1 and 2
basis.3 <- create.fourier.basis(rangeval=c(0,365),nbasis=109)
data_W.fd.3 <- Data2fd(y = data_W,argvals = time,basisobj = basis.3)
plot.fd(data_W.fd.3)

# estimate of the mean and of the covariance kernel
library(fields)
X11(width = 14, height = 7)
layout(cbind(1,2))
plot.fd(data_W.fd.1)
lines(mean(data_W.fd.1),lwd=2)
eval.1 <- eval.fd(time,data_W.fd.1)
image.plot(time,time,(cov(t(eval.1))[1:365,]))

X11(width = 14, height = 7)
layout(cbind(1,2))
plot.fd(data_W.fd.2)
lines(mean(data_W.fd.2),lwd=2)
eval.2 <- eval.fd(time,data_W.fd.2)
image.plot(time,time,(cov(t(eval.2))[1:365,]))

X11(width = 14, height = 7)
layout(cbind(1,2))
plot.fd(data_W.fd.3)
lines(mean(data_W.fd.3),lwd=2)
eval.3 <- eval.fd(time,data_W.fd.3)
image.plot(time,time,(cov(t(eval.3))[1:365,]))


#####################################################

# Second dataset: lip
# 51 measurements of the position of the lower lip every 7 
# milliseconds for 20 repitions of the syllable 'bob'.
data_L <- lip

layout(1)
time <- seq(0,350,by=7)
matplot(time,data_L,type='l',main='Lip data',ylab='Position',
        xlab='Time (millisec.)')

# smoothing. In this case, the data are already quite smooth. The smoothing
# does not pose problems

# Since the Fourier basis creates periodical data, it modifies the data
# near the boundaries of the domain, to render them periodical
basis <- create.fourier.basis(rangeval=c(0,350),nbasis=51)
data_L.fd <- Data2fd(data_L,time,basis)
plot.fd(data_L.fd)

# It would be better to use a b-spline basis
basis <- create.bspline.basis(rangeval=c(0,350),nbasis=21)
data_L.fd <- Data2fd(y = data_L,argvals = time,basisobj = basis)
plot.fd(data_L.fd)

# estimate of the mean and of the covariance function

X11(width = 14, height = 7)
layout(cbind(1,2))
plot.fd(data_L.fd,xaxs='i')
lines(mean(data_L.fd),lwd=2)
eval <- eval.fd(time,data_L.fd)
image.plot(time, time, (cov(t(eval))[1:51,]))

############################################################################################################

dev.off()

##### FPCA #####

# First dataset: canadian weather

# interpolated data (Choice 1)
layout(1)
plot.fd(data_W.fd.1,ylab='temperature')

pca_W.1 <- pca.fd(data_W.fd.1,nharm=5,centerfns=TRUE)

# scree plot
# pca.fd computes alle the 365 eigenvalues, but only the first 
# N-1=34 are non-null
X11(width = 7, height = 10)
layout(rbind(1,2))
plot(pca_W.1$values[1:34],xlab='j',ylab='Eigenvalues')
plot(cumsum(pca_W.1$values)[1:34]/sum(pca_W.1$values),xlab='j',ylab='CPV',ylim=c(0.8,1))

# first two FPCs
X11(width = 14, height =5)
layout(cbind(1,2))
plot(pca_W.1$harmonics[1,],col=1,ylab='FPC1')
plot(pca_W.1$harmonics[2,],col=2,ylab='FPC2')

# plot of the FPCs as perturbation of the mean
media<-mean(data_W.fd.1)

plot(media,lwd=2,ylim=c(-25,20),ylab='temperature',main='FPC1')
lines(media+pca_W.1$harmonics[1,]*sqrt(pca_W.1$values[1]), col=2)
lines(media-pca_W.1$harmonics[1,]*sqrt(pca_W.1$values[1]), col=3)
# Variation in amplitude (more in winter than in summer)

plot(media,lwd=2,ylim=c(-20,20),ylab='temperature',main='FPC2')
lines(media+pca_W.1$harmonics[2,]*sqrt(pca_W.1$values[2]), col=2)
lines(media-pca_W.1$harmonics[2,]*sqrt(pca_W.1$values[2]), col=3)
# temperate climate or not

# Command of the library fda that automatically does these plots
par(mfrow=c(1,2))
plot.pca.fd(pca_W.1, nx=100, pointplot=TRUE, harm=c(1,2), expand=0, cycle=FALSE)

graphics.off()

###
# data smooth (Choice 2)
layout(1)
plot.fd(data_W.fd.2)

pca_W.2<-pca.fd(data_W.fd.2,nharm=5,centerfns=TRUE)

# scree plot
X11(width = 7, height =10)
layout(rbind(1,2))
plot(pca_W.2$values,xlab='j',ylab='Eigenvalues')
plot(cumsum(pca_W.2$values)/sum(pca_W.2$values),xlab='j',ylab='CPV',ylim=c(0.8,1))

# first two FPCs
X11(width = 14, height =5)
layout(cbind(1,2))
plot(pca_W.2$harmonics[1,],col=1,ylab='FPC1')
plot(pca_W.2$harmonics[2,],col=2,ylab='FPC2')

# plot of the FPCs as perturbation of the mean
media<-mean(data_W.fd.2)

plot(media,lwd=2,ylim=c(-25,20),ylab='temperature',main='PC1')
lines(media+pca_W.2$harmonics[1,]*sqrt(pca_W.2$values[1]), col=2)
lines(media-pca_W.2$harmonics[1,]*sqrt(pca_W.2$values[1]), col=3)

plot(media,lwd=2,ylim=c(-20,20),ylab='temperature',main='PC2')
lines(media+pca_W.2$harmonics[2,]*sqrt(pca_W.2$values[2]), col=2)
lines(media-pca_W.2$harmonics[2,]*sqrt(pca_W.2$values[2]), col=3)

# similar interpretations as before

# scatter plot of the scores
x11(width = 14, height = 7)
par(mfrow=c(1,2))
plot(pca_W.1$scores[,1],pca_W.1$scores[,2],xlab="Scores FPC1",ylab="Scores FPC2",lwd=2)
points(pca_W.1$scores[35,1],pca_W.1$scores[35,2],col=2, lwd=4)

plot(pca_W.1$scores[,1],pca_W.1$scores[,2],type="n",xlab="Scores FPC1",
     ylab="Scores FPC2",xlim=c(-400,250))
text(pca_W.1$scores[,1],pca_W.1$scores[,2],dimnames(data_W)[[2]], cex=1)

# outlier: Resolute (35)

matplot(eval.1,type='l')
lines(eval.1[,35],lwd=4, col=2)

coord<-CanadianWeather$coordinates
plot(coord[,2:1],col=0)
text(coord[,2:1],rownames(coord))

# Exercise: perform FPCA with Choice 3 of smoothing

graphics.off()

###
# Second dataset: lip

layout(1)
plot.fd(data_L.fd)

pca_L<-pca.fd(data_L.fd,nharm=5,centerfns=TRUE)

# scree plot
X11(width = 7, height = 10)
layout(rbind(1,2))
plot(pca_L$values,xlab='j',ylab='Eigenvalues')
plot(cumsum(pca_L$values)/sum(pca_L$values),xlab='j',ylab='CPV',ylim=c(0.8,1))

# First three FPCs
X11(width = 21, height = 7)
layout(cbind(1,2,3))
plot(pca_L$harmonics[1,],col=1,ylab='FPC1')
plot(pca_L$harmonics[2,],col=2,ylab='FPC2')
plot(pca_L$harmonics[3,],col=3,ylab='FPC3')

# plot of the principal components as perturbation of the mean
media<-mean(data_L.fd)

plot(media,lwd=2,ylim=c(-10,12),main='FPC1')
lines(media+pca_L$harmonic[1,]*sqrt(pca_L$values[1]), col=2)
lines(media-pca_L$harmonic[1,]*sqrt(pca_L$values[1]), col=3)
# variation in amplitude in the centre

plot(media,lwd=2,ylim=c(-10,12),main='PC2')
lines(media+pca_L$harmonic[2,]*sqrt(pca_L$values[2]), col=2)
lines(media-pca_L$harmonic[2,]*sqrt(pca_L$values[2]), col=3)
# horizontal translation

plot(media,lwd=2,ylim=c(-10,12),main='PC3')
lines(media+pca_L$harmonic[3,]*sqrt(pca_L$values[3]), col=2)
lines(media-pca_L$harmonic[3,]*sqrt(pca_L$values[3]), col=3)
# variation in amplitude at the boundaries of the domain

# Command of the library fda that automatically does these plots
par(mfrow=c(1,3))
plot.pca.fd(pca_L, nx=100, pointplot=TRUE, harm=c(1,2,3), expand=0, cycle=FALSE)

# Scores
X11(width = 21, height = 7)
layout(cbind(1,2,3))
plot(pca_L$scores[,1],pca_L$scores[,2],xlab="Scores PC1",ylab="Scores PC2",lwd=2)
points(pca_L$scores[12,1],pca_L$scores[12,2],col=2, lwd=4)
points(pca_L$scores[9,1],pca_L$scores[9,2],col=3, lwd=4)
plot(pca_L$scores[,1],pca_L$scores[,3],xlab="Scores PC1",ylab="Scores PC3",lwd=2)
points(pca_L$scores[12,1],pca_L$scores[12,3],col=2, lwd=4)
points(pca_L$scores[9,1],pca_L$scores[9,3],col=3, lwd=4)
plot(pca_L$scores[,2],pca_L$scores[,3],xlab="Scores PC2",ylab="Scores PC3",lwd=2)
points(pca_L$scores[12,2],pca_L$scores[12,3],col=2, lwd=4)
points(pca_L$scores[9,2],pca_L$scores[9,3],col=3, lwd=4)

X11(width = 10, height = 7)
matplot(eval,type='l')
lines(eval[,12],lwd=4, col=2)
lines(eval[,9],lwd=4, col=3)

graphics.off()
