# Author: Christian Dudel
###############################################################################


sullivan <- function(survival,prevalence) {
	
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
