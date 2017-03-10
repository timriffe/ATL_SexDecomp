# TODO: Add comment
# 
# Author: tim
###############################################################################


# this requires reading in Data from the stem project ThanoEmpirical.
# The dat object LoessListCohrts5 is a complex list structure from said project,
# which will change in the future when a better-informed smoothing procedure
# is determined and implemented. For now, we have simple loess smoothes of the 
# given (scaled) variable over chronological age, thanatological age, and birth cohort.
# each birth cohort has a matrix with thanatological age in rows and chronological
# age in columns. Cell elements indicate prevalence. This method doesn't account
# for the panel structure of the data, treating each observation as independent.
# It also does not take account of universe boundaries, and so can have edge artifacts.
# It also does not take account of the original coding in any intelligent way (link functions)
# These are the three things that will be accounted for when a new method is 
# developed. Presumably the results will change, but probably not qualitatively.
# I think that the range of g(y) patterns in here is likely indicative of the range
# we would see in rigorously treated data.
Dat <- local(get(load("/home/tim/git/ThanoEmpirical/ThanoEmpirical/Data/LoessListCohrts5.Rdata")))

# get gy, need data object with columns for 
# Variable, Cohort, TTD, Prevalance

X <- Dat[[1]]
# Females
# thanatological
PrevTTD <- do.call(rbind,lapply(Dat, function(X){
			coh         <- X[["Female"]]$C[,1]
			names(X[["Female"]]$Surf) <- coh
			names(X[["Male"]]$Surf)   <- coh
			prevf        <- do.call(rbind, lapply(coh,function(coh., X.){
						x <- X.[["Female"]][["Surf"]][[as.character(coh.)]]
						
						# this is a simple mean of the valid observations within each
						# thanatological age. We don't weight by lifespan distributions.
						x <- rowMeans(x, na.rm = TRUE)
						x[x < 0 & !is.nan(x)] <- 0
						x[is.nan(x)] <- NA
						data.frame(Cohort = coh., TTD = names(x), Prev = x)
					}, X. = X))
	        prevf$Morbidity <- X[["Male"]]$varname
			prevf$Sex <- "f"

			prevm        <- do.call(rbind,lapply(coh,function(coh., X.){
								x <- X.[["Male"]][["Surf"]][[as.character(coh.)]]
								x <- rowMeans(x,na.rm=TRUE)
								x[x < 0 & !is.nan(x)] <- 0
								x[is.nan(x)] <- NA
								data.frame(Cohort = coh., TTD = names(x), Prev = x)
							}, X. = X))
			prevm$Morbidity <- X[["Male"]]$varname
			prevm$Sex <- "m"
			prev <- rbind(prevf,prevm)
			prev[,c("Morbidity","Cohort","Sex","TTD","Prev")]
		}))

rownames(PrevTTD) <- NULL

save(PrevTTD, file = "/home/tim/git/ATL_SexDecomp/ATL_SexDecomp/Data/gyALL.Rdata")


