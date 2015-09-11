# Author: riffe

# Alyson, skip to line 35

###############################################################################
if (system("hostname",intern=TRUE) %in% c("triffe-N80Vm", "tim-ThinkPad-L440")){
	# if I'm on the laptop
	setwd("/home/tim/git/APCT/APCT")
} else {
	if (system("hostname",intern=TRUE) == "PC-403478"){
		# on MPIDR PC
		setwd("U://git//APCT//APCT")
	} else {
		# in that case I'm on Berkeley system, and other people in the dept can run this too
		setwd(paste0("/data/commons/",system("whoami",intern=TRUE),"/git/APCT/APCT"))
	}
}
this.wd <- getwd()
# need to grab loess results from the ThanoEmpirical repository
TEwd    <- gsub(pattern="APCT",replacement="ThanoEmpirical",this.wd )
Results <- local(get(load(file.path(
								TEwd,
								"Data",
								"LoessQuinquenal_imp.Rdata"))))
names(Results)

# make 2 arrays, a male and a female.
# think some, then send to Alyson around when she gets back from vacay. Maybe Anna O too.
ADL3 <- Results[["adl3__0.5"]]


####################################################################################
# some code to do a first approximation of decomposing sex-differences in healthy life expectancy
# first load this object ADL3, the results of predicting some values from a loess
# model, fit to some HRS data.
#save(ADL3,file="Data/ADL3.Rdata")
ADL3 <- local(get(load("Data/ADL3.Rdata"))) # modify as necessary
# can't have negatives
imp0 <- function(x){
	x[x < 0]<- 0
	x
}
Males   <- imp0(ADL3$Male$Surf)
Females <- imp0(ADL3$Female$Surf)
# these are arrays. the third level is the 1915-1919 cohort.
# that's the one we'll work with.

# quick peek to figure out how to proceed.
# image(t(Females[,,3]))

# note that we're missing a piece of the upper triangle (bottom of matrix)
# too few data points there to get anything out of it.
# let's just impute it with something reasonable-seeming
# so we can move forward.

#image(t(Males[,,3]))
#a <- 70:100
#plot(a,Males["12",,3],type='l',col="green",ylim=c(0,.06))
#lines(a,Males["11",,3],col="royalblue")
#lines(a,Males["10",,3],col="blue")   # (looks similar for females)

# let's just apply the same proportional change that we see btwn
# the last two rows to each extrapolated row. If it drops
# below zero, we impute zero. Easy, cheap.
extendOmatic <- function(M){
	# some columns may have all NAs, be sure to keep it that way
	# assume data are in this very particular shape. OK, because
	# this is a once-off hack (we hope)
	start <- M[nrow(M), ]
	ratio <- start /  M[nrow(M) - 1, ]
	
	Col <- ncol(M)
	Row <- nrow(M)
	
	# the part we need to fill in with a triangle of extrapolated data
	NewBlock <- matrix(nrow = Col - Row, ncol = Col)
	
	for (i in 1:nrow(NewBlock)){
		NewBlock[i, ] <- start * ratio ^ i # cumprod
	}
	
	# append matrix
	Mout <- rbind(M, NewBlock)
	
	# need to respect NA diagonal.
	Nums <- row(Mout) + col(Mout)
	
	# need a marker on the right side in the top row
	Ind  <- max((1:ncol(Mout))[!is.na(Mout[1,])])
	Mout[Nums > Nums[1,Ind]] <- NA
	
	# clean extrap, just in case
	Mout[Mout < 0] <- 0
	Mout[Mout > 1] <- 1
	
	# include rownames
	rownames(Mout) <- 1:nrow(Mout) - 1
	Mout
}
# cool, now we can do stuff!
M1915e <- extendOmatic(Males[,,3])     
F1915e <- extendOmatic(Females[,,3])

# make sure NA placement matches, you never know...
all(colSums(!is.na(M1915e)) == colSums(!is.na(F1915e)))
all(rowSums(!is.na(M1915e)) == rowSums(!is.na(F1915e)))

# slice off unneeded columns (i.e. ages not observed for this cohort)
M1915e <- M1915e[, colSums(!is.na(M1915e)) > 0] 
F1915e <- F1915e[, colSums(!is.na(F1915e)) > 0] 
M1915e <- M1915e[rowSums(!is.na(M1915e)) > 0, ] 
F1915e <- F1915e[rowSums(!is.na(F1915e)) > 0, ] 

# chrono ages + thano ages, should be same lengths
(ca <- as.integer(colnames(M1915e)))
(ta <- as.integer(rownames(M1915e)))

#image(t(M1915e-F1915e))
#
#par(mfrow=c(1,2))
#image(t(M1915e), zlim = c(0,.6))
#image(t(F1915e), zlim = c(0,.6))

# females have higher values than males in ALL cells, wow.
# I guess that's consistent with lit?
sum(M1915e > F1915e, na.rm = TRUE)
sum(F1915e > M1915e, na.rm = TRUE)

#######################################################################
# next: get mortality data from HMD
# We want cohort Mx data for the 1915-1919 cohort,
# ages 74 - 95

# no need to run this: I copied results in-line here:
# library(HMDHFDplus)
#cMx <- readHMD("Data/USA_cMx_1x5.txt") # could also read from web with readHMDweb()
#cMx <- cMx[cMx$Year == 1915 & cMx$Age >= 74 &  cMx$Age <= 95 & !is.na( cMx$Age),]
#cMxm <- cMx$Male
#cMxf <- cMx$Female
#
#dput(cMxm)
cMxm <- c(0.049986, 0.052877, 0.057103, 0.061325, 0.066315, 0.07132, 
		0.081139, 0.088158, 0.096355, 0.104718, 0.113543, 0.122644, 0.133971, 
		0.144406, 0.15709, 0.169257, 0.183678, 0.199131, 0.217271, 0.239125, 
		0.263696, 0.285049)
#dput(cMxf)
cMxf <- c(0.028989, 0.031342, 0.034459, 0.037594, 0.041088, 0.045432, 
		0.05254, 0.057924, 0.064801, 0.071759, 0.079933, 0.088454, 0.09763, 
		0.107266, 0.118591, 0.129736, 0.142486, 0.157165, 0.173349, 0.191498, 
		0.212159, 0.237106)

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

qxm <- mx2qx(cMxm)
qxf <- mx2qx(cMxf)
lxm <- qx2lx(qxm)
lxf <- qx2lx(qxf)
# these could of course be rescaled such that the radix is the 1992 population of people from
# the 1915-1919 cohorts. But I think it's not necessary here.
#plot(74:96,lxm,type='l')
#lines(74:96,lxf,col = "blue")

# dx is the real sugar for working with these data.
plot(lx2dx(lxm),type='l')
lines(lx2dx(lxf),col="red")

############################################################################
# note: the people alive in x are the people that die in x+

# we learn that there is another direction in which we ought to extend
# the data. Since we only have up to age 95, there are people alive in l(x)
# that we don't have in this window of d(x)...

# That last open age group of d(x)
# is therefore present in each preceding age of l(x), and if we want to see the
# chrono consequences of these morbidity patterns then we need to be able to 
# compose an l(x) out of it's various d(x) components, taking into account
# the thanatological trajectories of each d(x). So we need values for those 
# missing d(x), for the ages > 95. I'm inclined to just repeat the experience of those
# that died at age 95, shifting it right. Then the matrix will get bigger.
# another problem will be that we need to extend cMx until the cohort
# is mostly extinct. In that case, how about doing the HMD extendo-O-matic
# thing and borrowing from the preceding cohorts.

# in short, we need to 'close out' both the mortality data and the morbidity data. 
# In the end, this will have little leverage on results. If mortality were
# a lot lower, it would have leverage. Realistically, but the time mortality
# ever does get that low, we'll actually have the data and not need to extrapolate!!

############################################################
# OK, step 1, for fudging a working example:
# extend our morbidity surfaces to include people out to omega.
# let's say omega is 110. Because HMD.

openm       <- rev(diag(M1915e[nrow(M1915e):1,]))
openf       <- rev(diag(F1915e[nrow(F1915e):1,]))
# again, for the sake of closing out. More thought could go into this.
anew        <- 96:110
RightSide   <- matrix(nrow=nrow(M1915e),ncol=length(anew),dimnames=list(rownames(M1915e),anew))
M1915ee     <- cbind(M1915e, RightSide)
F1915ee     <- cbind(F1915e, RightSide)
N           <- ncol(M1915e)
for (i in 1:ncol(RightSide)){
	for (z in 1:length(openm)){
		M1915ee[z,N+1+i-z] <- openm[z]
		F1915ee[z,N+1+i-z] <- openf[z]
	}
}

M1915ee <- rbind(M1915ee,M1915ee[1:15,] * 0)
F1915ee <- rbind(F1915ee,F1915ee[1:15,] * 0)

M1915ee[upper.tri(M1915ee)[nrow(M1915ee):1,]] <- NA
F1915ee[upper.tri(F1915ee)[nrow(F1915ee):1,]] <- NA
rownames(F1915ee) <- rownames(M1915ee) <- 1:nrow(M1915ee) - 1
#####################################
# PS: this kind of extrap could maybe just use Lee Carter?
# or just spline extrap? Extrap over logit, then expit back?
# or force it out to an asymptote? For this particular case
# there is no substantive precedent I can think of, given 
# the dimensions. Note that the above double-loop extrapolation
# thingy only affects 5% of a given chrono age for males (since
# .05 were in 96+ of d96). 11% for females. So the leverage of
# decision making on this piece is at least bounded by those
# percentages. We could put in all 1s or all 0s, that is, and 
# see the unreasonable bounds to the consequences.
#####################################
# now extend Mx to 110+
#####################################
plot(log(cMxm))
plot(log(cMxf))
# cohort Mx is really dang linear, good!
# I'm going to extrapolate with a line...
# Jim wouldn't like that. Whatever. In this case
# we just want to close out. In the end, we just 
# want to draw conclusions about ages 74-95...
# these decisions boil back to the tail of d(x) for us,
# which contributes to the pattern in each preceding age.

coefm <- lm(log(cMxm)~c(74:95))$coef
coeff <- lm(log(cMxf)~c(74:95))$coef

# predict ages 96-110...
cMxme <- exp(coefm[1] + coefm[2] * 96:110)
cMxfe <- exp(coeff[1] + coeff[2] * 96:110)

# looks good.
#plot(74:95, cMxm,log='y', type = 'l', col = "blue",xlim=c(74,110),ylim=c(.02,1))
#lines(74:95, cMxf,col="red")
#points(96:110,cMxme,pch=19,col="blue" )
#points(96:110,cMxfe,pch=19,col="red" )

# so we'll redo to get the final d(x)
cMxmE <- c(cMxm,cMxme)
cMxfE <- c(cMxf,cMxfe)

lxm <- qx2lx(mx2qx(cMxmE))
lxf <- qx2lx(mx2qx(cMxfE))

dxm <- lx2dx(lxm)
dxf <- lx2dx(lxf)

# groups the last two dx elements so everything both adds and is conformable
dxm <- c(dxm[1:(length(dxm)-2)],sum(rev(dxm)[1:2]))
dxf <- c(dxf[1:(length(dxf)-2)],sum(rev(dxf)[1:2]))
# cool, we're closed out now.
#plot(dxm)
#plot(dxf)

##################################
# now we can do demography!
##################################

# Just assume we have two valid ATL patterns, 
# one for males, the other for females.
# yellow = hotter than red in this ramp....

#par(mfrow=c(1,2))
#image(74:110, 0:37, t(M1915ee), zlim = range(c(M1915ee, F1915ee), na.rm = TRUE), main = "males")
#image(74:110, 0:37, t(F1915ee), zlim = range(c(M1915ee, F1915ee), na.rm = TRUE), main = "females")

dxM      <- dxm[row(M1915ee) + col(M1915ee) - 1]
dim(dxM) <- dim(M1915ee)
dxF      <- dxf[row(F1915ee) + col(F1915ee) - 1]
dim(dxF) <- dim(F1915ee)

# verify that the column sums match lx:
all(colSums(dxM,na.rm=TRUE) == lxm[-length(lxm)])
all(colSums(dxF,na.rm=TRUE) == lxf[-length(lxf)])

# so, we can multiply in neat ways:
Mtc <- dxM * M1915ee
Ftc <- dxF * F1915ee

# the column sums of this morbidity prevalence matrix give the chronological
# age pattern of prevalence within the cohort (Assuming this radix)

Mc <- colSums(Mtc, na.rm=TRUE) / lxm[-length(lxm)]
Fc <- colSums(Ftc, na.rm=TRUE) / lxf[-length(lxf)]

#plot(74:110, Mc, type = 'l', ylim = c(0, .6), main = "marginal chronological age pattern (artifactual in this case)")
#lines(74:110, Fc, col = "red")

# Part of the difference between males and females is due to the underlying
# morbidity pattern, which translates to the difference in the distribution
# of lifelines. 

# Further (what I tried to demonstrate in the presentation),
# if one holds the morbidity pattern constant, but changes
# the underlying mortality, we will draw different conclusions about the expected
# change in prevalence than if we simply assume the marginal chronological age pattern
# and then perturb mortality. 
# leverage gained by this added insight into morbidity variation could help
# better predict convergence/divergence between males and females, for instance.


# what are the healthy and unhealthy life expectancies at age 74? (odd age, but whatever)

lx2Lx <- function(lx){
	(lx + c(lx[-1],0)) / 2
	
}

lx2ex <- function(lx){
	lx <- lx / lx[1]
	Lx <- lx2Lx(lx)
	sum(Lx)
}


lx2ex(lxm) # 2.5 years difference
lx2ex(lxf)

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
# than the other columns
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
	lx2ex(lx) - geteUx(mx, Morb)
}




getPropH <- function(mx, Morb){
	geteHx(mx,Morb) / lx2ex(qx2lx(mx2qx(mx)))
}
# males higher proportion healthy
geteUx(cMxmE,M1915ee) # males, unhealthy
geteUx(cMxfE,F1915ee) # females, unhealthy
geteHx(cMxmE,M1915ee) # males, healthy
geteHx(cMxfE,F1915ee) # females, healthy

# proportion healthy
getPropH(cMxmE,M1915ee)  # males
getPropH(cMxfE,F1915ee)  # females

#############################################
# playing more:
#############################################
# How would this change if males had female mortality?
getPropH(cMxfE,M1915ee) # it would decrease a little

# male mort with female morb
getPropH(cMxmE,F1915ee) # big decrease

# probably better to look at years of unhealthy / healthy instead of prop
# healthy years, male morb w female mort:
geteHx(cMxfE,M1915ee)
geteUx(cMxfE,M1915ee)
# female morb with male mort is the worst!
geteHx(cMxmE,F1915ee)
geteUx(cMxmE,F1915ee)

# increase male longevity with no change in morb pattern
# here implies greater prop of life healthy (but about the same)
getPropH(cMxmE, M1915ee)
getPropH(cMxmE * .8, M1915ee)

# same result for females
getPropH(cMxfE, F1915ee)
getPropH(cMxfE * .8, F1915ee)

# most improvement goes to healthy years
# females:
geteHx(cMxmE, M1915ee)
geteHx(cMxmE * .8, M1915ee)
geteUx(cMxmE, M1915ee)
geteUx(cMxmE * .8, M1915ee)

# same:
geteHx(cMxfE, F1915ee)
geteHx(cMxfE * .8, F1915ee)
geteUx(cMxfE, F1915ee)
geteUx(cMxfE * .8, F1915ee)
##########################################################
# this is all fun, but how to decompose?
##########################################################
# *use Horiuchi-Pletcher-Wilmoth pseudo continuous decomposition
# but what quantity do we decompose?
# I guess the difference in healthy years of life?

#library(devtools)
#install_github("timriffe/DecompHoriuchi/DecompHoriuchi")
library(DecompHoriuchi) # the key function is DecompHoriuchiOrig()

# need to give a function that takes a single argument as a vector
vec2obj <- function(vec,mxind=1:37,matdim=c(37,37)){
	mx        <- vec[mxind]
	morb      <- vec[-mxind]
	dim(morb) <- matdim
	list(mx=mx,morb=morb)
}
getHxvec <- function(vec,mxind=1:37,matdim=c(37,37)){
	obj  <- vec2obj(vec,mxind,matdim)
	geteHx(obj$mx, obj$morb)
}

# this does the decomp (computationally intensive...). A direct analytic decomp method
# could be devised, but this was handy.
output <- DecompContinuousOrig(func = getHxvec, 
		rates1 = c(cMxmE, c(M1915ee)),
		rates2 = c(cMxfE, c(F1915ee)),
		N = 20,
		mxind = 1:37,
		matdim = c(37, 37))
diffs <- vec2obj(output)

##############################################
# what results do we have?
plot(74:110, diffs$mx, type = 'l', ylim = c(-.05,.3), main = "mortality accounts for most sex difference\nin healthy life expectancy here")
lines(74:110, colSums(diffs$morb), col = "red")
abline(h=0)
rect(95,-1,120,1,col = "#00000020", border=NA)
# errors add to the sex gap in healthy e(74), with small resid
# one could reduce error to zero by changing the decomp function to not go over midpoints, but rather over interval
# endpoints, but it's a trivial error, so who cares.

# what about remaining lifetime contributions?
plot(0:36,rowSums(diffs$morb),type='l')
rect(12,-1,12,1, border=NA, col ="#00000020")
# appear stronger and more concentrated than marginal chrono contributions

# and disagregated morbidity contributions to difference:
filled.contour(74:110,0:36,t(diffs$morb))

# notice how despite sex differences in our assumed patterns of morbidity beyond age 100
# these have no effect on the sex differences in healty life expectancy because mortality
# gives them little weight. If mortality were to improve, then these cells would obtain weight!
# it's a cool interaction.

# note: these data are best thought of as test data, for playing. I'm awaiting some statistical
# refinement for how to read these patterns out of the data in the first place. Also, these 
# morbidity extrapolations didn't have much weight here, but their influence could also be 
# displayed in the graph with some extra computation.

# note: things would work differently for morbidity conditions with other patterns of variation.
# a decomposition like this sees through this diversity in variation, whereas as other standard
# methods do not filter it out properly. So, there's a lot to learn still.