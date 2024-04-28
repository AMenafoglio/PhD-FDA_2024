################################################################################
# Short Course: "Statistical methods of data science -                         #
#                An introduction to Functional Data Analysis"                  #
# Università degli Studi di Bergamo                                            #
# Alessandra Menafoglio, MOX-Department of Mathematics, Politecnico di Milano  # 
################################################################################

##--------------------------##
## Functional Linear Models ##
##--------------------------##

setwd("~/Documents/Politecnico/Didattica/2024 PhD Bergamo/Codes-FDA")

library(fda)
data("CanadianWeather")

matplot(CanadianWeather$dailyAv[,,1],type='l') # temperature
matplot(CanadianWeather$dailyAv[,,2],type='l') # precipitation
matplot(CanadianWeather$dailyAv[,,3],type='l') # log precipitation

# fda ------------------------------------------------------------------
#  data are in Canadian Weather object
#  print the names of the data
print(names(CanadianWeather))


###  functional response with vector explanatory variables  
###
help(fRegress)
daybasis65 <- create.fourier.basis(rangeval=c(0, 365), nbasis=65,
                                   axes=list('axesIntervals'))
Temp.fd <- with(CanadianWeather, smooth.basisPar(day.5,
                                                 dailyAv[,,'Temperature.C'], 
                                                 daybasis65)$fd)
TempRgn.f <- fRegress(Temp.fd ~ region, CanadianWeather)
plot(TempRgn.f$betaestlist$const)
plot(TempRgn.f$betaestlist$region.Atlantic)
plot(TempRgn.f$betaestlist$region.Continental)
plot(TempRgn.f$betaestlist$region.Pacific)

# estimated mean of the four regions:
Arctic.mean = eval.fd(evalarg = day.5,fdobj=TempRgn.f$betaestlist$const$fd)
Atlantic.mean = Arctic.mean + eval.fd(evalarg = day.5,fdobj=TempRgn.f$betaestlist$region.Atlantic$fd)
Continental.mean = Arctic.mean + eval.fd(evalarg = day.5,fdobj=TempRgn.f$betaestlist$region.Continental$fd)
Pacific.mean = Arctic.mean + eval.fd(evalarg = day.5,fdobj=TempRgn.f$betaestlist$region.Pacific$fd)

region.means = rbind(t(Arctic.mean),t(Atlantic.mean),t(Continental.mean),t(Pacific.mean))
Temp.eval = eval.fd(day.5,Temp.fd)
regions = CanadianWeather$region

matplot((Temp.eval),col=factor(regions),type='l',lty=1,lwd=0.5)
matlines(t(region.means),type='l',lwd=2,lty=1)

###
###
###  functional response with functional explanatory variable  
###  concurrent model
###
LogPrec.fd <- with(CanadianWeather, smooth.basisPar(day.5,
                                                 dailyAv[,,'log10precip'], 
                                                 daybasis65)$fd)


# without penalization
Temp.Prec.regress.conc <- fRegress(LogPrec.fd ~ Temp.fd)
plot(Temp.Prec.regress.conc$betaestlist$const)
plot(Temp.Prec.regress.conc$betaestlist$Temp.fd)


# with penalization
betabasis = daybasis65
beta0Par = fdPar(betabasis,lambda=1e5)
beta1Par = fdPar(betabasis,lambda=1e7)
betaList = list(beta0Par, beta1Par)

Temp.Prec.regress.conc <- fRegress(LogPrec.fd ~ Temp.fd,betalist = betaList)
plot(Temp.Prec.regress.conc$betaestlist[[1]])
plot(Temp.Prec.regress.conc$betaestlist[[2]])


### total model
betabasis = daybasis65
tempBeta1fd = bifd(matrix(0,65,65), betabasis, betabasis)

beta0Par = fdPar(betabasis,lambda=1e-5)
beta1stPar = bifdPar(bifdobj=tempBeta1fd, 
                     lambdas=1e3,lambdat=1e3) #1e10

betaList = list(beta0Par, beta1stPar)
linmodSmooth = linmod(LogPrec.fd, Temp.fd, betaList)

afd <- linmodSmooth$beta0estfd   #  The intercept function
bfd <- linmodSmooth$beta1estbifd    #  The bivariate regression function

#  plot the intercept function

plot(afd, xlab="Day", ylab="Intercept function")

#  plot the regression function as a surface

bfdmat <- eval.bifd(weeks, weeks, bfd)
persp(weeks, weeks, bfdmat, xlab="Day(t)", ylab="Day(s)")


library(fields)
image.plot(weeks, weeks, bfdmat)

library(plotly)
fig <- plot_ly(z= bfdmat) %>% add_surface()
fig


