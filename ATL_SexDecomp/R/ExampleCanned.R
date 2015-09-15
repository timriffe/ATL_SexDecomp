 
# Author: tim
###############################################################################

x <- 70:110
y <- 0:40
gy <- c(seq(1,0,length=10),rep(0,32))

# plot(y,gy,type='l')

gompmx <- function(a,b,x){
	a * exp(b*x)
}

# male a,b
am <- .0000015
bm <- .125

af <- .0000005
bf <- .135

plot(x,gompmx(am,bm,x),type='l',log='y')
lines(x,gompmx(af,bf,x),col="red")
# some helper functions:
mx2qx <- function(mx){
	mx / (1 + .5 * mx)
}


qx2lx <- function(qx){
	cumprod(1-c(0,qx))
}

lx2dx <- function(lx){
	-diff(c(lx,0))
}

lx2Lx <- function(lx){
	(lx + c(lx[-1],0)) / 2
	
}

lx2ex <- function(lx){
	lx <- lx / lx[1]
	Lx <- lx2Lx(lx)
	sum(Lx)
}
library(magrittr)
mx2ex <- function(mx){
	mx %>% mx2qx %>% qx2lx %>% lx2ex
}

mxf <- gompmx(af,bf,x)
mxm <- gompmx(am,bm,x)
lxf <- mx2lx(mxf)
lxm <- mx2lx(mxm)
dxf <- mx2dx(mxf)
dxm <- mx2dx(mxm)

#plot(x,mx2lx(gompmx(am,bm,x))[-1], type = 'l')
#lines(x, mx2lx(gompmx(af,bf,x))[-1],col="red")
#

#mx2ex(gompmx(am,bm,x))
#mx2ex(gompmx(af,bf,x))

############################


#plot(x,mx2dx(gompmx(am,bm,x))[-1], type = 'l')
#lines(x, mx2dx(gompmx(af,bf,x))[-1],col="red")

day <- function(dx){
	N <- length(dx)
	dmat <- matrix(0,N,N)
	dxy <- dx[row(dmat)+col(dmat) - 1]
	dim(dxy) <- dim(dmat)
	dxy[is.na(dxy)] <- 0
	dxy
}

#dayf <- day(mx2dx(gompmx(af,bf,x)))
#rowSums(dayf) - colSums(dayf)

# turn vector of dx into triangle that conforms with morbidity.
# ncol(Morb) must equal length(dx)
dxweight <- function(dx, Morb){
	stopifnot(length(dx) == ncol(Morb))
	dxM      <- dx[row(Morb) + col(Morb) - 1]
	dim(dxM) <- dim(Morb)
	dxM
}

# weight morbidity triangle by dx (density of lifelines)
getMorbWeighted <- function(dx, Morb){
	dxweight(dx, Morb) * Morb
}

# unhealthy expectancy
# best to do everything straight from mx, because it's cooler to perturb it rather
# than the other columns #mx <- mxf
geteUx <- function(mx, Morb){
	lx  <- qx2lx(mx2qx(mx))
	dx  <- lx2dx(lx)
	Mwx <- colSums(getMorbWeighted(dx, Morb), na.rm = TRUE)
	Lx  <- lx2Lx(lx)
	Lx  <- Lx / lx[1]
	if (length(Lx) > length(Mwx)){
		Mwx <- c(Mwx, Mwx[length(Mwx)])
	}
	sum(Mwx * Lx)
}


# healthy expectancy
# likewise we build from mx
geteHx <- function(mx, Morb){
	lx  <- qx2lx(mx2qx(mx))
	mx2ex(mx) - geteUx(mx, Morb)
}

Morb <- matrix(gy,42,42)

geteHx(mxf,Morb)

library(magrittr)

mxf %>% mx2qx %>% qx2lx %>% lx2ex -> exf
mxm %>% mx2qx %>% qx2lx %>% lx2ex -> exm

eHf <- geteHx(mxf,Morb)  
eHm <- geteHx(mxm,Morb)

eHf / exf 
eHm / exm 

exf - eHf
exm - eHm

library(DecompHoriuchi)
mxf
mxgy2eH <- function(rates){
	mx   <- rates[1:41]
	gy   <- rates[-c(1:41)]
	Morb <- matrix(gy, length(gy), length(gy))
	geteHx(mx,Morb) 
}

nogy <- DecompContinuousOrig(mxgy2eH, rates1 = c(mxm,gy), rates2 = c(mxf,gy),N=20)
mxcontrib <- nogy[1:41]
gycontrib <- nogy[-c(1:41)]

plot(mxcontrib)

gym <- c(seq(.4,0,length=6),rep(0,36))
gyf <- c(seq(.5,0,length=11),rep(0,31))
length(gy)
gy
plot(gym, ylim=c(0,1))
lines(gyf)
gydiff <- DecompContinuousOrig(mxgy2eH, rates1 = c(mxm,gym), rates2 = c(mxf,gyf),N=20)
mxcontrib <- gydiff[1:41]
gycontrib <- gydiff[-c(1:41)]
plot(gycontrib)

mxMorb2eH <- function(rates){
	mx    <- rates[1:41]
	Morb   <- rates[-c(1:41)]
	dim(Morb) <- rep(sqrt(length(Morb)),2)
	geteHx(mx,Morb) 
}

closeoutmat <- function(Morbin){
	Morbin[is.na(Morbin)] <- 0	
	dimswewant <- matrix(0,42,42,dimnames = list(0:41,70:111))
	dimswewant[rownames(Morbin), colnames(Morbin)] <- Morbin
	dimswewant
}

# 1915
M1915 <- closeoutmat(Results[["m"]][["1915"]][["adl5_"]])
F1915 <- closeoutmat(Results[["f"]][["1915"]][["adl5_"]])

Sex1915diff <- DecompContinuousOrig(
		mxMorb2eH, 
		rates1 = c(mxm,c(M1915)), 
		rates2 = c(mxf,c(F1915)),
		N=20)
mxcontrib <- Sex1915diff[1:41]
Morbcontrib <- Sex1915diff[-c(1:41)]
dim(Morbcontrib) <- c(42,42)

par(mfrow=c(1,2))
plot(70:111,colSums(Morbcontrib),ylim=c(-.1,.02))
plot(0:41, rowSums(Morbcontrib),ylim=c(-.1,.02))
source("R/SurfMap.R")
graphics.off()
filled.contour(70:111,0:41,t(Morbcontrib),asp=1)

sum(mxcontrib)
sum(Morbcontrib)

plot(mxcontrib)

mH <- mxMorb2eH(c(mxm,c(M1915)))
fH <- mxMorb2eH(c(mxf,c(F1915)))

Sex1915diff2 <- DecompContinuousOrig(
		mxMorb2eH, 
		rates1 = c(mxm*.5,c(M1915)), 
		rates2 = c(mxf*.5,c(F1915)),
		N=20)
mxcontrib2 <- Sex1915diff2[1:41]
Morbcontrib2 <- Sex1915diff2[-c(1:41)]
dim(Morbcontrib2) <- c(42,42)

graphics.off()
filled.contour(70:111,0:41,t(Morbcontrib),asp=1,zlim=c(.001,-.01))
dev.new()
filled.contour(70:111,0:41,t(Morbcontrib2),asp=1,zlim=c(.001,-.01))

sum(Morbcontrib2);sum(Morbcontrib)
sum(mxcontrib2);sum(mxcontrib)

sum(c(Morbcontrib,mxcontrib))
sum(c(Morbcontrib2,mxcontrib2))