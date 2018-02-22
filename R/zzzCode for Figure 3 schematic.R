# loading data and preliminaries

library(plyr)
library(dplyr)
library(magrittr)
options(scipen=3,stringsAsFactors=FALSE)
setwd("C:\\Dropbox\\sex decomp Tim\\Data\\")

# HRS data
load("gyall.Rdata")
PrevTTD$TTD <- as.numeric(PrevTTD$TTD)

# HMD data from three countries (USA, JPN, and CAN)
DxA <- read.table("C:\\hmd_statistics\\deaths\\Deaths_1x1\\USA.Deaths_1x1.txt",
                  head=T,skip=2)
NxA <- read.table("C:\\hmd_statistics\\exposures\\Exposures_1x1\\USA.Exposures_1x1.txt",
                  head=T,skip=2)
DxB <- read.table("C:\\hmd_statistics\\deaths\\Deaths_1x1\\JPN.Deaths_1x1.txt",
                  head=T,skip=2)
NxB <- read.table("C:\\hmd_statistics\\exposures\\Exposures_1x1\\JPN.Exposures_1x1.txt",
                  head=T,skip=2)
DxC <- read.table("C:\\hmd_statistics\\deaths\\Deaths_1x1\\CAN.Deaths_1x1.txt",
                  head=T,skip=2)
NxC <- read.table("C:\\hmd_statistics\\exposures\\Exposures_1x1\\CAN.Exposures_1x1.txt",
                  head=T,skip=2)


# Building the output for the figure panels
# Panels are numbered left to right, top to bottom,
# Panel 1 = blanck

#-------------------------------------------------------------------
# PANEL 2: Smoothed death densities corresponding to e60 = 20,25,30
#-------------------------------------------------------------------

# Building smooth death densities assuming Gompertz mortality, based roughly
# on empirical data for remaining life expectancies at age 60 of 20, 25, and 30

# some helper functions to get us from a Gompertz fit to a death density

mu.gompertz <- function(alpha, beta, x) {
  return(alpha * exp(beta*x))
}

loglike.poisson.gompertz <- function(theta, dx, nx, x) {
  alpha <- theta[1]
  beta <- theta[2]
  answer <- -(sum(dx * log(mu.gompertz(x=x, alpha=alpha, beta=beta)) -
                    mu.gompertz(x=x, alpha=alpha, beta=beta)*nx))
  return(answer)
}

mx2qx <- function(mx){
  mx / (1 + .5 * mx)
}

qx2lx <- function(qx){
  cumprod(1-c(0,qx))
}

lx2dx <- function(lx){
  -diff(c(lx,0))
}

mx2dx <- function(mx){
  mx %>% mx2qx %>% qx2lx %>% lx2dx
}


# Death density for e60=20, based on USA males 2002
Dx <- filter(DxA, Year==2002)$Male[51:91]
Nx <- filter(NxA, Year==2002)$Male[51:91]
gomp.est <- optim(par=c(0.00001, 0.14), fn=loglike.poisson.gompertz,
                  dx=Dx, nx=Nx, x=50:90)
Mx.ex20 = mu.gompertz(gomp.est$par[1],gomp.est$par[2],60:109) # 19.9 years
Dx.ex20 = mx2dx(Mx.ex20)

# Death density for e60=25, based on CAN females 2005
Dx <- filter(DxC, Year==2005)$Female[51:91]
Nx <- filter(NxC, Year==2005)$Female[51:91]
gomp.est <- optim(par=c(0.00001, 0.14), fn=loglike.poisson.gompertz,
                  dx=Dx, nx=Nx, x=50:90)
Mx.ex25 = mu.gompertz(gomp.est$par[1],gomp.est$par[2],60:109) # 25.0 years
Dx.ex25 = mx2dx(Mx.ex25)

# Death density for e60=30, based roughly on JPN females 2014
# with minor adjustments to a and b, to increase e60 from 28.7 (highest empirically observed) to 30 years
Dx <- filter(DxB, Year==2014)$Female[51:91]
Nx <- filter(NxB, Year==2014)$Female[51:91]
gomp.est <- optim(par=c(0.00001, 0.14), fn=loglike.poisson.gompertz,
                  dx=Dx, nx=Nx, x=50:90)
Mx.ex29 = mu.gompertz(gomp.est$par[1],gomp.est$par[2],60:109) # 28.7 years
Mx.ex30 = mu.gompertz(gomp.est$par[1]*0.78,gomp.est$par[2]*1.01,60:109) # 30.0 years
Dx.ex30 = mx2dx(Mx.ex30)

# plot
Age = 60:110
plot(Age,Dx.ex20,t="l",lwd=2,xlab="Age",ylab="Death density",
     ylim=c(0,0.045),col="darkorange")
lines(Age,Dx.ex25,lwd=2,col="darkred")
lines(Age,Dx.ex30,lwd=2,col="darkblue")
grid()


#-------------------------------------------------------------------
# PANEL 3,5,7: Proportion disabled by time-to-death 
#-------------------------------------------------------------------

# These three plots estimate the proportion disabled by time to death
# for 3 different disabling processes, all of which are experienced
# by half of the population at the time of death, but which differ
# in the timing of onset prior to death and in the steepness of the
# curve with the approach to death. The ???rst type of disability is virtually
# nonexistent 5 years prior to death, but then increases very rapidly as death
# approaches. The middle variant of disability is rare 15 years before death,
# but increases to about 20 percent of the population 5 years before death
# and rises sharply thereafter. The bottom ???gure depicts a disabling process
# that although still strictly determined by time-to-death, is common and
# accumulates very slowly starting from about 50 years before death.


# To base these curves on reality, we created the middle variant to 
# approach difficulties in performing at least one of 5 IADLs. 


# averaging over cohorts
meanprev <- ddply(PrevTTD,.(Morbidity,Sex,TTD), summarize, Morb=mean(Prev,na.rm=T))

#--!!!!check definition of iadl5 and change section in text from bathing 
# to iadl5!!!!

x=filter(meanprev,Morbidity=="iadl5_",Sex=="f")$TTD
y=filter(meanprev,Morbidity=="iadl5_",Sex=="f")$Morb
plot(x,y,xlab="TTD",ylab="Morbidity prevalence",ylim=c(0,0.6),col="red")

# fiting an exponential curve through it with intercept
yfit = 0.5*exp(-0.2*x)
lines(x,yfit)

# pretty close, so we'll use this line as the mid variant, but push things out by 50 years
# plot 1, steep variant
x <- 0:50
ysteep = 0.5*exp(-0.7*x)
plot(x,ysteep,xlab="TTD",ylab="Morbidity prevalence",ylim=c(0,0.5),t="l",lwd=2)

# plot 2, mid-variant, modeled on IADL5
ymed = 0.5*exp(-0.2*x)
plot(x,ymed,xlab="TTD",ylab="Morbidity prevalence",ylim=c(0,0.5),t="l",lwd=2)

# plot 3, gentle variant, prevalence is still above 10% 15 years before death
ygentle = 0.5*exp(-0.05*x)
plot(x,ygentle,xlab="TTD",ylab="Morbidity prevalence",ylim=c(0,0.5),t="l",lwd=2)



#-------------------------------------------------------------------
# PANEL 4,6,8: Interaction of death density and proportion disabled
# by time-to-death, i.e. g*(a,y) schedules that apply to g(y) and d(x)
#--------------------------------------------------------------------





# combining Dx and Age in data frames for female data at different e60 levels
DF20 <- data.frame(Dx=Dx.ex20,Age=60:110)
DF25 <- data.frame(Dx=Dx.ex25,Age=60:110)
DF30 <- data.frame(Dx=Dx.ex30,Age=60:110)

# vectors of remaining lifespans  
ages = 60:110 
rls20 <- rls25 <- rls30 <-list()
for (i in 1:length(ages)){
  rls20[[i]] <- (DF20$Dx[DF20$Age>=ages[i]])/sum(DF20$Dx[DF20$Age>=ages[i]])  
  rls25[[i]] <- (DF25$Dx[DF25$Age>=ages[i]])/sum(DF25$Dx[DF25$Age>=ages[i]])  
  rls30[[i]] <- (DF30$Dx[DF30$Age>=ages[i]])/sum(DF30$Dx[DF30$Age>=ages[i]])  
}

# glueing NAs to the end of each vector to make them the same length
rls.wNA20 <- lapply(rls20, function(x) {c(x, rep(NA, length(ages) - length(x)))})
rls.wNA25 <- lapply(rls25, function(x) {c(x, rep(NA, length(ages) - length(x)))})
rls.wNA30 <- lapply(rls30, function(x) {c(x, rep(NA, length(ages) - length(x)))})

# putting remaining lifespans into matrix format
rlsmat20 <- do.call(rbind, rls.wNA20)
rlsmat25 <- do.call(rbind, rls.wNA25)
rlsmat30 <- do.call(rbind, rls.wNA30)

# calculating the disability prevalence by chronological age 
# for steep, medium and gentle disability curves for each remaining life
# expectancy level



Age = 60:111
TTD = 0:50
DP20.stp <- colSums(t(rlsmat20[,1:51])*ysteep,na.rm=T)
DP20.med <- colSums(t(rlsmat20[,1:51])*ymed,na.rm=T)
DP20.gen <- colSums(t(rlsmat20[,1:51])*ygentle,na.rm=T)

DP25.stp <- colSums(t(rlsmat25[,1:51])*ysteep,na.rm=T)
DP25.med <- colSums(t(rlsmat25[,1:51])*ymed,na.rm=T)
DP25.gen <- colSums(t(rlsmat25[,1:51])*ygentle,na.rm=T)

DP30.stp <- colSums(t(rlsmat30[,1:51])*ysteep,na.rm=T)
DP30.med <- colSums(t(rlsmat30[,1:51])*ymed,na.rm=T)
DP30.gen <- colSums(t(rlsmat30[,1:51])*ygentle,na.rm=T)


# example plot
Age <- 60:110
plot(Age,DP20.stp,t="l",lwd=2,xlab="Age",ylab="Disability prevalence",
     ylim=c(0,0.5),
     col="darkorange",cex.axis=0.9)
lines(Age,DP25.stp,lwd=2,col="darkred")
lines(Age,DP30.stp,lwd=2,col="darkblue")
grid()
text(70,0.4,bquote(paste(e[60]," = 20")),col="darkorange",cex=1.5)
text(70,0.3,bquote(paste(e[60]," = 25")),col="darkred",cex=1.5)
text(70,0.2,bquote(paste(e[60]," = 30")),col="darkblue",cex=1.5)


# Now creating the whole schematic

pdf("C:\\Dropbox\\sex decomp Tim\\postWP\\schematic3.pdf",width=8,height=8)

Age <- 60:110
layout(matrix(c(1:8), 4, 2, byrow=TRUE), 
       widths=c(lcm(6),lcm(5.7)), 
       heights=c(lcm(4),lcm(3.7),lcm(3.7),lcm(3.7)))
#layout.show(8)
par(mar=c(1,3,1,1))
plot.new()

#panel 2
par(mar=c(1,1,1,1))
plot(Age,Dx.ex20,t="l",lwd=2,xlab="",ylab="",ylim=c(0,0.045),col="darkorange",cex.axis=0.9)
lines(Age,Dx.ex25,lwd=2,col="darkred")
lines(Age,Dx.ex30,lwd=2,col="darkblue")
grid()
mtext("Death density",side=2,line=2.5,outer=F,cex=0.8)

#panel 3
par(mar=c(1,3,1,1))
plot(TTD[1:21],ysteep[1:21],t="l",lwd=2,xlab="",ylab="",
     ylim=c(0,0.5),col="cornflowerblue",cex.axis=0.9)
grid()
mtext("Proportion disabled",side=2,line=2.5,outer=F,cex=0.8)

#panel 4
par(mar=c(1,1,1,1))
plot(Age,DP20.stp,t="l",lwd=2,xlab="",ylab="",ylim=c(0,0.5),col="darkorange",cex.axis=0.9)
lines(Age,DP25.stp,lwd=2,col="darkred")
lines(Age,DP30.stp,lwd=2,col="darkblue")
grid()
text(70,0.4,bquote(paste(e[60]," = 20")),col="darkorange",cex=1.5)
text(70,0.3,bquote(paste(e[60]," = 25")),col="darkred",cex=1.5)
text(70,0.2,bquote(paste(e[60]," = 30")),col="darkblue",cex=1.5)

#panel 5
par(mar=c(1,3,1,1))
plot(TTD[1:21],ymed[1:21],t="l",lwd=2,xlab="",ylab="",ylim=c(0,0.5),col="cornflowerblue",cex.axis=0.9)
grid()
mtext("Proportion disabled",side=2,line=2.5,outer=F,cex=0.8)

#panel 6
par(mar=c(1,1,1,1))
plot(Age,DP20.med,t="l",lwd=2,xlab="",ylab="",ylim=c(0,0.5),col="darkorange",cex.axis=0.9)
lines(Age,DP25.med,lwd=2,col="darkred")
lines(Age,DP30.med,lwd=2,col="darkblue")
grid()

#panel 7
par(mar=c(1,3,1,1))
plot(TTD[1:21],ygentle[1:21],t="l",lwd=2.5,xlab="",ylab="",ylim=c(0,0.5),col="cornflowerblue",cex.axis=0.9)
grid()
mtext("Proportion disabled",side=2,line=2.5,outer=F,cex=0.8)
mtext("Time to Death",side=1,line=3,outer=F,cex=0.8)

#panel 8
par(mar=c(1,1,1,1))
plot(Age,DP20.gen,t="l",lwd=2,xlab="",ylab="",ylim=c(0,0.5),col="darkorange",cex.axis=0.9)
lines(Age,DP25.gen,lwd=2,col="darkred")
lines(Age,DP30.gen,lwd=2,col="darkblue")
grid()
mtext("Age",side=1,line=3,outer=F,cex=0.8)


dev.off()

