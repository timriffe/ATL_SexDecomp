
# Author: tim
###############################################################################
library(magrittr)
library(DecompHoriuchi)


# Gompertz mx given a,b, and age vector
gompmx <- function(a,b,x){
	a * exp(b*x)
}

# other gompertz exact and approximate things,
# not used. Instead we go to lifetable mode
gompdx <- function(a,b,x){
	a*exp(b*x)*exp(-a/b *(exp(b*x)-1))
}
gomplx <- function(a,b,x){
	exp(-a/b*(exp(b*x)-1))
}
gompex <- function(a,b,x1,x2=Inf){
	integrate(gomplx, lower = x1, upper = x2, a = a, b = b)$value
}

gompLx <- function(a,b,x){
	N <- length(x)
	Lx <- rep(NA, N)
	ag <- c(x,Inf)
	for (i in 1:N){
		Lx[i] <- gompex(a,b,ag[i],ag[i+1])
	}
	Lx
}

# this works unless mx is high (infants, very old)
# some helper functions.

# no need for a0 because we assume starting age later in life, e.g. 60
mx2qx <- function(mx){
	mx / (1 + .5 * mx)
}

# cumprod of px
qx2lx <- function(qx){
	cumprod(1-c(0,qx))[1:length(qx)]
}

#
lx2dx <- function(lx){
	-diff(c(lx,0))
}

lx2Lx <- function(lx,closeout=0){
	(lx + c(lx[-1],closeout)) / 2
	
}

lx2ex <- function(lx){
	lx <- lx / lx[1]
	Lx <- lx2Lx(lx)
	sum(Lx)
}

mx2ex <- function(mx){
	mx %>% mx2lx %>% lx2ex
}
mx2dx <- function(mx){
	mx %>% mx2lx %>% lx2dx
}
mx2Lx <- function(mx){
	mx %>% mx2lx %>% lx2Lx
}
mx2lx <- function(mx){
	mx %>% mx2qx %>% qx2lx
}
#plot(x,mx2lx(gompmx(am,bm,x))[-1], type = 'l')
#lines(x, mx2lx(gompmx(af,bf,x))[-1],col="red")
#

#mx2ex(gompmx(am,bm,x))
#mx2ex(gompmx(af,bf,x))

############################
# morbidity by ttd and age, given a vector either chrono or thano
# dx is for dimensions. This function just makes a conformable matrix
# assuming you don't have one to start with. dx <- dxl;morb <- gy
May <- function(morb, dx, chrono = FALSE){
	N 				        <- length(dx)
	
	if (is.null(dim(morb))){
		if (chrono){
			stopifnot(length(morb) == N)
			Morb 			<- matrix(morb, N, N, byrow = TRUE)
		} else {
			if (length(morb) < N){
				morb        <- c(morb,rep(0, N - length(morb)))
			}
			Morb 			<- matrix(morb, N, N)
		}
	} else {
		# in this case we uglyly assign 0s to values outside the range.
		# to do something fancier, make Morb yourself!
		Morb                <- matrix(0, N, N)
		Morb[1:nrow(morb), 1:ncol(morb)] <- morb
	}
	Morb
}
#plot(x,mx2dx(gompmx(am,bm,x))[-1], type = 'l')
#lines(x, mx2dx(gompmx(af,bf,x))[-1],col="red")

day <- function(dx){
	# rescale to be sure
	dx 				<- dx / sum(dx)
	N 				<- length(dx)
	dmat 			<- matrix(0, N, N)
	day 			<- dx[row(dmat) + col(dmat) - 1]
	dim(day) 		<- dim(dmat)
	day[is.na(day)] <- 0
	# all(colSums(day) == rowSums(day))
	day
}

#dayf <- day(mx2dx(gompmx(af,bf,x)))
#rowSums(dayf) - colSums(dayf)

# Morb must either be a thano x chrono matrix
# of morbidity prevalence (cross classified)
# or it can be a thanatological vector. If you want 
# to model a strictly chronological process then
# you'd be the 2d matrix where values are equal within
# columns.
#dxweight <- function(dx, Morb){
#	
#	if (is.null(dim(Morb))){
#		# then we assume Morb is a thanatological vector
#		Morb <- matrix(Morb,ncol=length(dx),nrow=length(Morb))
#	}
#	stopifnot(length(dx) == ncol(Morb))
#	dxM      <- dx[row(Morb) + col(Morb) - 1]
#	dim(dxM) <- dim(Morb)
#	dxM
#}

# weight morbidity triangle by dx (density of lifelines)
getMorbWeighted <- function(dx, Morb){
    Day <- day(dx)
	
	if(!all(dim(Day) == dim(Morb))){
		Morb <- May(Morb,dx)
	}
	Day * Morb
}
#plot(colSums(getMorbWeighted(dxf,Morb), na.rm = TRUE))
#lines(mx2lx(mxf))
#lines(mx2Lx(mxf))
#plot(colSums(getMorbWeighted(dxf,Morb), na.rm = TRUE)/mx2lx(mxf))
# unhealthy expectancy
# best to do everything straight from mx, because it's cooler to perturb it rather
# than the other columns #mx <- mxf
# this gives morb prevalence that lines up with lx, not Lx.
# need to average in the same way to translate to Lx and get
# expectancy. This one is OK for plotting ga, but it conforms
# with lx, not Lx, and shouldn't be used for expectancies.
getMorbagelx <- function(mx, Morb){
	
	dx  <- mx2dx(mx)
	Mwx <- colSums(getMorbWeighted(dx, Morb), na.rm = TRUE)
	lx  <- mx2lx(mx)
	
	if (length(lx) > length(Mwx)){ # this is now guaranteed 
		Mwx <- c(Mwx, Mwx[length(Mwx)])
	}
	Mwx / lx
}
# This one is adequate for calculating expectancies (multiply into Lx)
mxgay2gaLx <- function(mx, Morb){
	N 			<- length(mx)
	mwx   		<- getMorbagelx(mx, Morb)
	Mwx 		<- lx2Lx(mwx, closeout = mwx[N])
	Mwx 
}

mxgy2gaLx <- function(mx, gy){
	Morb 		<- May(morb = gy, dx = mx, chrono = FALSE)
	Mwx  		<- mxgay2gaLx(mx, Morb)
	Mwx
}
geteM <- function(mx, Morb){
	Lx  <- mx2Lx(mx)
	Mwx <- mxgay2gaLx(mx, Morb)
	sum(Mwx * Lx)
}

# prob: the morbidity today is a function of the mortality in the future?
# then we can't use period mortality to weight current morbidity?

# healthy expectancy
# likewise we build from mx
geteH <- function(mx, Morb){
	mx2ex(mx) - geteM(mx, Morb)
}












