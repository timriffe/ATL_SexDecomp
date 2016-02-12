 
# Author: tim
###############################################################################
setwd("/home/tim/git/ATL_SexDecomp/ATL_SexDecomp")
library(HMDHFDplus)
library(reshape2)
source("R/Functions.R")
#xyz <-"CHE"
#XYZ <- getHMDcountries()
#ab <- do.call(rbind,lapply(XYZ,function(xyz, us, pw){
#					cat(xyz,"\n")
#			Mx   <- readHMDweb(xyz,"Mx_1x1",username=us,password=pw)
# # would be better to weight by log exposure, but then there are
# # tons of ways this could be better...
# # fit to ages 40-95.
#			Mxm <- acast(Mx,Age~Year,value.var = "Male")
#			Mxf <- acast(Mx,Age~Year,value.var = "Female")
#			Mxm[Mxm == 0] <- NA
#			Mxf[Mxf == 0] <- NA
#			Mxm <- Mxm[,colSums(is.na(Mxm[41:96,])) == 0]
#			Mxf <- Mxf[,colSums(is.na(Mxf[41:96,])) == 0]
#			
#			abM <- t(apply(Mxm, 2, function(.Mx){
#								b <- 40:95
#								lm(log(.Mx[b + 1]) ~ b)$coef
#							}))
#			abF <- t(apply(Mxf, 2, function(.Mx){
#								b <- 40:95
#								lm(log(.Mx[b + 1]) ~ b)$coef
#							}))
#			
#			
#			abM[, 1] 		<- exp(abM[, 1] )
#			abF[, 1] 		<- exp(abF[, 1] )
#			colnames(abM) 	<- c("a", "b")
#			colnames(abF) 	<- c("a", "b")
#			abM 			<- as.data.frame(abM)
#			abF 			<- as.data.frame(abF)
#			rbind(abM,abF)
#		}, us = us, pw = pw))
#
#par(mfrow = c(1, 2))
#plot(1933:2013, abM[1, ], ylab = "a", type = 'l', col = "blue", ylim = c(0, .0005))
#lines(1933:2013, abF[1, ],col="red")
#
#plot(1933:2013,  abM[2, ], ylab = "b", type = 'l', col = "blue", ylim = c(0.06, .1))
#lines(1933:2013, abF[2, ], col = "red")
#range(ab$b)
#bnew <- seq(min(ab$b),max(ab$b),length=100)
#plot(ab[,1],ab[,2],xlab="a",ylab="b",pch=19,
#		col="#0000FF10")
#lines(predict(loess(a~b,data=abM,span=1),data.frame(b=bnew)),bnew,col="blue")
#points(abF[,1],abF[,2],xlab="a",ylab="b",pch=19,
#		col="#FF000050")
#lines(predict(loess(a~b,data=abF,span=1),data.frame(b=bnew)),bnew,col="red")
#
#plot(ab[,1],ab[,2],xlab="a",ylab="b",pch=19,
#		col="#0000FF10")
#abloess <- loess(a~b,data=ab,span=.4)
#getafromb <- function(b,mod){
#	predict(mod,data.frame(b=b))
#}
#lines(predict(loess(a~b,data=ab,span=.4),data.frame(b=bnew)),bnew,col="red")


JPN <- readHMDweb("JPN","mltper_1x1",username=us,password=pw)

##############################################
x  <- 70:110
y  <- 0:40

## male a,b
#bh <- .05
#ah <- .002
##getafromb(.07,abloess)
#bl <- .065
#al <- .0005
#
#plot(gompmx(ah,bh,x),type='l',col="red",log='y')
#lines(gompmx(al,bl,x),col="blue")
#plot(x,gompmx(am,bm,x),type='l',log='y')
#lines(x,gompmx(af,bf,x),col="red")




#mxl <- gompmx(al,bl,x)
#mxh <- gompmx(ah,bh,x)
mx2 <- JPN$mx[JPN$Year == 2010 & JPN$Age %in% x]
mx1 <- JPN$mx[JPN$Year == 1980 & JPN$Age %in% x]
lx2 <- mx2lx(mx2)
lx1 <- mx2lx(mx1)
dx2 <- mx2dx(mx2)
dx1 <- mx2dx(mx1)
Lx2 <- mx2Lx(mx2)
Lx1 <- mx2Lx(mx1)
gy  <- c(seq(.5,0,length=10),rep(0,31))
# slow ramp up
gy  <- cumprod(rep(.9,41)) * seq(.9,0,length=41)^8
plot(gy)

Morb <- May(gy,dx2)

# remaining life expectancy at age 70
(ex2 <- mx2ex(mx2))
(ex1  <- mx2ex(mx1))
(eH2 <- geteHx(mx2,Morb)  )
(eH1 <- geteHx(mx1,Morb))

(eH1 / ex1 )
(eH2 / ex2) # higher proportion

# unhealthy, similar but also increased
(ex1 - eH1)
(ex2 - eH2)

################################################
# buy what if we had used the sullivan method based on 
# 1980 age pattern of morbidity fixed, rather than ttd pattern fixed?

ga1         <- mxgay2gaLx(mx1, Morb)
eU2sullivan <- sum(ga1 * Lx2)
eHsullivan  <- ex2 - eU2sullivan

# sullivan increase in unhealthy expectancy
100 * eU2sullivan / (ex1 - eH1) # sullivan predicts 52% increase
100 * (ex2 - eH2) / (ex1 - eH1) # real 6% increase

##############################################################
# OLDER CODE, written on train between Prague and Rostock, now superceded by Functions.R
#library(magrittr)
#library(DecompHoriuchi)
## plot(y,gy,type='l')
#
#gompmx <- function(a,b,x){
#	a * exp(b*x)
#}
#
## some helper functions:
#mx2qx <- function(mx){
#	mx / (1 + .5 * mx)
#}
#
#
#qx2lx <- function(qx){
#	cumprod(1-c(0,qx))
#}
#
#lx2dx <- function(lx){
#	-diff(c(lx,0))
#}
#
#lx2Lx <- function(lx){
#	(lx + c(lx[-1],0)) / 2
#	
#}
#
#lx2ex <- function(lx){
#	lx <- lx / lx[1]
#	Lx <- lx2Lx(lx)
#	sum(Lx)
#}
#
#mx2ex <- function(mx){
#	mx %>% mx2qx %>% qx2lx %>% lx2ex
#}
#mx2dx <- function(mx){
#	mx %>% mx2qx %>% qx2lx %>% lx2dx
#}
#mx2Lx <- function(mx){
#	mx %>% mx2qx %>% qx2lx %>% lx2Lx
#}
##plot(x,mx2lx(gompmx(am,bm,x))[-1], type = 'l')
##lines(x, mx2lx(gompmx(af,bf,x))[-1],col="red")
##
#
##mx2ex(gompmx(am,bm,x))
##mx2ex(gompmx(af,bf,x))
#
#############################
#
#
##plot(x,mx2dx(gompmx(am,bm,x))[-1], type = 'l')
##lines(x, mx2dx(gompmx(af,bf,x))[-1],col="red")
#
#day <- function(dx){
#	N <- length(dx)
#	dmat <- matrix(0,N,N)
#	dxy <- dx[row(dmat)+col(dmat) - 1]
#	dim(dxy) <- dim(dmat)
#	dxy[is.na(dxy)] <- 0
#	dxy
#}
#
##dayf <- day(mx2dx(gompmx(af,bf,x)))
##rowSums(dayf) - colSums(dayf)
#
## turn vector of dx into triangle that conforms with morbidity.
## ncol(Morb) must equal length(dx)
#dxweight <- function(dx, Morb){
#	stopifnot(length(dx) == ncol(Morb))
#	dxM      <- dx[row(Morb) + col(Morb) - 1]
#	dim(dxM) <- dim(Morb)
#	dxM
#}
#
## weight morbidity triangle by dx (density of lifelines)
#getMorbWeighted <- function(dx, Morb){
#	dxweight(dx, Morb) * Morb
#}
#
## unhealthy expectancy
## best to do everything straight from mx, because it's cooler to perturb it rather
## than the other columns #mx <- mxf
## this gives
#getMorbage <- function(mx, Morb){
#	lx  <- qx2lx(mx2qx(mx))
#	dx  <- lx2dx(lx)
#	Mwx <- colSums(getMorbWeighted(dx, Morb), na.rm = TRUE)
#	Lx  <- lx2Lx(lx)
#	Lx  <- Lx / lx[1]
#	if (length(Lx) > length(Mwx)){
#		Mwx <- c(Mwx, Mwx[length(Mwx)])
#	}
#	Mwx / Lx
#}
#mxgy2eMwx <- function(mx,gy){
#	Morb <- matrix(gy, length(gy), length(gy))
#	Mwx  <- getMorbage(mx, Morb)
#	Mwx
#}
#geteUx <- function(mx, Morb){
#	Lx <- mx2Lx(mx)
#	Mwx <- getMorbage(mx, Morb)
#	sum(Mwx * Lx)
#}
#
## prob: the morbidity today is a function of the mortality in the future?
## then we can't use period mortality to weight current morbidity?
#
## healthy expectancy
## likewise we build from mx
#geteHx <- function(mx, Morb){
#	lx  <- qx2lx(mx2qx(mx))
#	mx2ex(mx) - geteUx(mx, Morb)
#}
#
#mxgy2eH <- function(rates){
#	mx   <- rates[1:41]
#	gy   <- rates[-c(1:41)]
#	Morb <- matrix(gy, length(gy), length(gy))
#	geteHx(mx,Morb) 
#}
#
#nogy <- DecompContinuousOrig(mxgy2eH, rates1 = c(mxm,gy), rates2 = c(mxf,gy),N=20)
#mxcontrib <- nogy[1:41]
#gycontrib <- nogy[-c(1:41)]
#####################################
#
#
#
#
#
##
#plot(mxcontrib)
#plot(mx2dx(mxm))
#lines(mx2dx(mxf))
#lines(mxcontrib)
#plot(getMorbage(mxf, Morb))
#lines(getMorbage(mxm, Morb))
#gym <- c(seq(.4,0,length=6),rep(0,36))
#gyf <- c(seq(.5,0,length=11),rep(0,31))
#
##plot(gym, ylim=c(0,1))
##lines(gyf)
#gydiff <- DecompContinuousOrig(mxgy2eH, rates1 = c(mxm,gym), rates2 = c(mxf,gyf),N=20)
#mxcontrib <- gydiff[1:41]
#gycontrib <- gydiff[-c(1:41)]
#ehM<- mxgy2eH(c(mxm,gym));ehF<- mxgy2eH(c(mxf,gyf))
#ehM;ehF
#ehF - ehM
#sum(mxcontrib);sum(gycontrib)
#sum(mxcontrib)+sum(gycontrib)
#
#
#mxMorb2eH <- function(rates){
#	mx    <- rates[1:41]
#	Morb   <- rates[-c(1:41)]
#	dim(Morb) <- rep(sqrt(length(Morb)),2)
#	geteHx(mx,Morb) 
#}
#
#closeoutmat <- function(Morbin){
#	Morbin[is.na(Morbin)] <- 0	
#	dimswewant <- matrix(0,42,42,dimnames = list(0:41,70:111))
#	dimswewant[rownames(Morbin), colnames(Morbin)] <- Morbin
#	dimswewant
#}
#
## 1915
#Results <- local(get(load("/home/tim/git/ATL_SexDecomp/ATL_SexDecomp/Data/ResultsIADL_ADL.Rdata")))
#M1915 <- closeoutmat(Results[["m"]][["1915"]][["adl5_"]])
#F1915 <- closeoutmat(Results[["f"]][["1915"]][["adl5_"]])
#
#Sex1915diff <- DecompContinuousOrig(
#		mxMorb2eH, 
#		rates1 = c(mxm,c(M1915)), 
#		rates2 = c(mxf,c(F1915)),
#		N=20)
#mxcontrib <- Sex1915diff[1:41]
#Morbcontrib <- Sex1915diff[-c(1:41)]
#dim(Morbcontrib) <- c(42,42)
#
#sum(mxcontrib)
#sum(Morbcontrib)
#
#
#
#
#par(mfrow=c(1,2))
#plot(70:111,colSums(Morbcontrib),ylim=c(-.1,.02))
#plot(0:41, rowSums(Morbcontrib),ylim=c(-.1,.02))
#source("R/SurfMap.R")
#graphics.off()
#filled.contour(70:111,0:41,t(Morbcontrib),asp=1)
#
#sum(mxcontrib)
#sum(Morbcontrib)
#
#plot(mxcontrib)
#
#mH <- mxMorb2eH(c(mxm,c(M1915)))
#fH <- mxMorb2eH(c(mxf,c(F1915)))
#
#Sex1915diff2 <- DecompContinuousOrig(
#		mxMorb2eH, 
#		rates1 = c(mxm*.5,c(M1915)), 
#		rates2 = c(mxf*.5,c(F1915)),
#		N=20)
#mxcontrib2 <- Sex1915diff2[1:41]
#Morbcontrib2 <- Sex1915diff2[-c(1:41)]
#dim(Morbcontrib2) <- c(42,42)
#
#graphics.off()
#filled.contour(70:111,0:41,t(Morbcontrib),asp=1,zlim=c(.001,-.01))
#dev.new()
#filled.contour(70:111,0:41,t(Morbcontrib2),asp=1,zlim=c(.001,-.01))
#
#sum(Morbcontrib2);sum(Morbcontrib)
#sum(mxcontrib2);sum(mxcontrib)
#
#sum(c(Morbcontrib,mxcontrib))
#sum(c(Morbcontrib2,mxcontrib2))
#
#
#gym;mxm
#mx2Lx <- function(mx){
#	mx %>% mx2qx %>% qx2lx %>% lx2Lx
#}
#plot(70:111,mxgy2eMwx(mxm, gym)/mx2Lx(mxm), type='l')
#lines(70:111,mxgy2eMwx(mxm * .75, gym)/mx2Lx(mxm*.75))
#lines(70:111,mxgy2eMwx(mxm * .5, gym)/mx2Lx(mxm*.5))
