
# Author: tim
###############################################################################
setwd("/home/tim/git/HLETTD")
library(HMDHFDplus)
LT <- readHMDweb("USA","mltper_1x1",username=us,password=pw)
mx <- LT$mx[LT$Year == 2010]

library(Matrix)
source("R/Functions.R")
source("R/SullivanDudel.R")

sdiag <- function(x,shift=0){
	out <- matrix()
}
# now how to get prev as function of gy using matrices. In the first
# instance, assuming stationary survival.
#install.packages("matrixcalc")
library(matrixcalc)

# go back to matrix from vec
vecinv <- function(x){
	m <- sqrt(length(x))
	stopifnot(floor(m)==m)
	dim(x) <- c(m,m)
	x
}

makeSquare <- function(ind = 1,X){
	Trans <- X * 0
	m     <- unique(dim(X))
	
	Trans[row(X) == col(X) + (ind - 1)] <- 1
	Trans[col(X) == row(X) + m - ind + 1] <- 1
	Trans
}

vec.lower.tri.up <- function(x){
	X <- matrix(0,sqrt(length(x)),sqrt(length(x)))
	inds <- 1:ncol(X)
	Perm <- as.matrix(Matrix::bdiag(lapply(inds,makeSquare,X)))
	vecinv(c(t(Perm) %*% x))
}

lower.tri.up <- function(X){
	vec.lower.tri.up(c(X))
}

getprev <- function(survival, gy){
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
	
	# each column gives the conditional survival.
	N  <- solve(diag(1,s+1)-U)
	
	# get conditional remaining lifetime distribution
	D  <- -diff(rbind(N,0))
	D  <- D * lower.tri(D,TRUE)
	
	# permute to useable form:
	DD <- lower.tri.up(D)
	
	# colSums...
	prevalence      <- t(DD) %*% c(gy,0)
	c(prevalence)
}

sullivany <- function(survival, gy) {
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
	
	# each column gives the conditional survival.
	N  <- solve(diag(1,s+1)-U)
	
	# get conditional remaining lifetime distribution
	D  <- -diff(rbind(N,0))
	D  <- D * lower.tri(D,TRUE)
	
	# permute to useable form:
	DD <- lower.tri.up(D)

	prevalence      <- c(c(gy,0) %*% DD)
	# prevalence
	# in subdiagonal, then closed out with 0s
	R1 <- diag(prevalence)
	R1 <- cbind(R1,rep(0,s))
	R1 <- rbind(rep(0,s+1),R1) 
	R1 <- cbind(R1,rep(0,s+1))
	R1 <- rbind(R1,rep(0,s+2))
	 
	 # 1s in subdiagonal
	Z  <- cbind(diag(1,s+1),rep(0,s+1))
	
	### Expectation
#	dim(N);dim(Z);dim(P);dim(Z %*% t(P*R1));dim(t(N) %*% Z )
	rho1 <- t(N) %*% Z %*% t(P*R1) %*% rep(1,s+2)
	
	# Output
	return(rho1)
}





n            <- length(mx)
survival     <- 1 - mx2qx(mx)
gy           <- cumprod(rep(.85,n))^4
prevalence   <- getprev(survival, gy)

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

pdf("Figures/gyTim.pdf")
plot(0:110,gy,type = 'l',xlab="Time to death",ylab="prevalence")
dev.off()