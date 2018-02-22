
# Author: tim
###############################################################################

# at one point we used Nusselder decomposition to make our point
# Then Shkolnikov and Andreev published this (poorly typeset) MPIDR Working Paper
# the_decomposition_of_the_difference_between_two_healthy_life_expectancies_which_formula_is_right_5924.htm

# they also give a spreadsheet comparing Nusselder w Andreev (2002), which we now use.
# I'd like to compare WHP vs Andreev

# for TR computers...otherwise you need to set wd yourself
if (system("hostname",intern=TRUE) %in% c("triffe-N80Vm","tim-ThinkPad-L440")){
	# if I'm on the laptop
	setwd("/home/tim/git/HLETTD")
} 
getwd()
# install.packages("lubridate")

Results               <- local(get(load("Data/resultsP.Rdata")))
library(reshape2)
head(melt(Results))

# get prevalence by age, hypothetical.
Results[[1]]$Male$Surf[,,"1915"]
ga <- colMeans(Results[[1]]$Male$Surf[,,"1915"],na.rm=TRUE)
gy <- rowMeans(Results[[1]]$Male$Surf[,,"1915"],na.rm=TRUE)
library(HMDHFDplus)
mlt <- readHMDweb("USA","mltper_1x1",us,pw)
mx  <- mlt$mx[mlt$Year == 2015]
source("R/Functions.R")
a <- 0:110

pi1 <- mxgy2gaLx(mx, gy)
pi2 <- mxgy2gaLx(mx / 2, gy)

#
pi1 <- 1
pi2 <- 1


qx1 <- mx2qx(mx)
qx2 <- mx2qx(mx / 2)
lx1 <- qx2lx(qx1)
lx2 <- qx2lx(qx2)
Lx1 <- lx2Lx(lx1)
Lx2 <- lx2Lx(lx2)
ex1 <- rev(cumsum(rev(Lx1))) / lx1
ex2 <- rev(cumsum(rev(Lx2))) / lx2
Px1 <- Lx1 / lx1
Px2 <- Lx2 / lx2
hLx1 <- pi1 * Lx1
hLx2 <- pi2 * Lx2
Thx1 <- rev(cumsum(rev(hLx1)))
Thx2 <- rev(cumsum(rev(hLx2)))
hx1  <- Thx1 / lx1
hx2  <- Thx2 / lx2
hx1[1];ex1[1];hx2[1];ex2[1]

lxavg <- (lx1 + lx2) / 2
Pxavg <- (Px2 + Px1) / 2
piavg <- (pi1 + pi2) / 2
# Andreev method:
Mort <- lxavg * piavg * (Px2 - Px1)  +
		# mean lx * Lx diff * mean prev
		1 / 2 * (lx2 * c(hx1[-1],0.2362796) + lx1 *  c(hx2[-1],0.2362796)) * (qx1 - qx2)
             # 1/2 (surv 2 * HLE1 + surv1 * HLE2) * (mort rate diff)
#Morb <- 1 / 4 * (lx1 + lx2) * (Px1 + Px2) * (pi2 - pi1)
Morb <- lxavg * Pxavg * (pi2 - pi1)

Mort2 <- 1 / 2 * lxavg * (Px2 - Px1) * (pi1 + pi2) +
		# mean lx * Lx diff * mean prev
		1 / 2 * (lx2 * hx1 + lx1 *  hx2) * (qx1 - qx2)
# 1/2 (surv 2 * HLE1 + surv1 * HLE2) * (mort rate diff)
plot(a, Mort, ylim = c(-.05,.05))
lines(a, Morb)

plot(Px2)

library(DecompHoriuchi)

mxpi2DLY <- function(mxpivec){
	N   <- length(mxpivec)
	mx  <- mxpivec[1:(N/2)]
	pie <- mxpivec[(N/2 + 1):N]
	lx  <- mx2lx(mx)
	Lx  <- lx2Lx(lx)
	Thx <- rev(cumsum(rev(Lx*pie)))
	HLE <- Thx / lx
	HLE[1]
}

dec <- DecompTest(mxpi2DLY, rates1 = c(mx, pi1), rates2 = c(mx/2, pi2),500)

plot(a, dec[1:111], ylim = c(-.05,.05),type = 'l',col = "blue")
lines(a, dec[112:222],col = "blue")
lines(a, Mort,col = "red")
lines(a, Morb,col="red")

lines(a, Mort2,col = "red",lty=2)

mxpi2DLY(c(mx/2, pi2)) - mxpi2DLY(c(mx, pi1))

sum(DecompTest(mxpi2DLY,rates1,rates2,500))


DecompTest <- function(func,rates1,rates2,N,...){
	# number of interval jumps   
	#y1 			<- func(rates1,...)
	#y2 			<- func(rates2,...)
	#N1 <- N + 1
	d 			<- rates2 - rates1
	n 			<- length(rates1)
	delta 		<- d / N
	#x <- rates1 + d * matrix(rep(.5:(N - .5) / N, n),byrow=TRUE,ncol=N)
	x <- rates1 + d * matrix(rep(((1:N)-.5) / N, n),byrow=TRUE,ncol=N)
	cc <- matrix(0,nrow=n,ncol=N)
	for (j in 1:N){
		for (i in 1:n){
			z 		<- rep(0,n)
			z[i] 	<- delta[i]/2
			cc[i,j] <- func((x[,j]+z),...)-func((x[,j]-z),...)
		}
	}	
	#cat("Error of decomposition (in %)\ne =",100*(sum(cc)/(y2-y1)-1),"\n")
	return(rowSums(cc))
}


#ex1[1] - ex2[1]
#plot(0:109,
#diff(cumsum(Lx2) - cumsum(Lx1) + (lx2 - lx1) * ex1),
#type='l')
#
#lines(0:110,lx2*(ex2-ex1) - c(lx2[-1],0) * (c(ex2[-1],0) - c(ex1[-1],0)))
