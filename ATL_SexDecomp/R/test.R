
# Author: tim
###############################################################################
setwd("/home/tim/git/ATL_SexDecomp/ATL_SexDecomp")
library(HMDHFDplus)
LT <- readHMDweb("USA","mltper_1x1",username=us,password=pw)
mx <- LT$mx[LT$Year == 2010]


source("R/Functions.R")
source("R/SullivanDudel.R")


# now how to get prev as function of gy using matrices. In the first
# instance, assuming stationary survival.

sullivany <- function(mx, gy) {
	prevalence <- mxgay2gaLx(mx, gy)
	survival   <- 1- mx2qx(mx)
	s <- length(survival)
	
	# Transition matrix
	# survival in subdiag
	U <- diag(survival)
	U <- cbind(U,rep(0,s))
	U <- rbind(rep(0,s+1),U)
	
	# Probability of dying
	M <- 1-colSums(U)
	
	# Transition matrix
	# add mort to bottom (qx), and 0s on right, except final 1
	P <- rbind(U,M)
	P <- cbind(P,c(rep(0,s+1),1))
	
	# prevalence
	# in subdiagonal, then closed out with 0s
	R1 <- diag(prevalence)
	R1 <- cbind(R1,rep(0,s))
	R1 <- rbind(rep(0,s+1),R1) 
	R1 <- cbind(R1,rep(0,s+1))
	R1 <- rbind(R1,rep(0,s+2))
	
	# each column gives the conditional survival.
	N  <- solve(diag(1,s+1)-U)
	# 1s in subdiagonal
	Z  <- cbind(diag(1,s+1),rep(0,s+1))
	
	### Expectation
#	dim(N);dim(Z);dim(P);dim(Z %*% t(P*R1));dim(t(N) %*% Z )
	rho1 <- t(N) %*% Z %*% t(P*R1) %*% rep(1,s+2)
	
	# Output
	return(rho1)
}

survival     <- 1 - mx2qx(mx)
gy           <- cumprod(rep(.85,N))^4
prevalence   <- mxgay2gaLx(mx, gy)

# slightly different results. It's because
# of the markov approximations.
sullivan(survival, prevalence)[1]
sum(prevalence * LT$Lx[LT$Year == 2010]/1e5)
# identical for stationary population
sullivany(mx, gy) -
sullivan(1-mx2qx(mx), prevalence)

# not identical if mortality changes
sullivany(mx*.8, gy) -
		sullivan(1-mx2qx(mx*.8), prevalence)


ref <- matrix(dx,111,111,dimnames=list(0:110,0:110))
ref <- ref * lower.tri(ref,TRUE)
#LexisUtils::LexisMap(ref,log=FALSE)
plot(dx)

# point made!
# however, this is only partially valid as an estimation technique,
# because by admitting that health is due to time until death, we
# are forced to look to the future. Since mortality is changing, we
# can assume that 

