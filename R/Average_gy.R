
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

# TR: June 18, all of these artifacts now fixed!

#Dat <- local(get(load("/home/tim/git/ThanoEmpirical/ThanoEmpirical/Data/LoessListCohrts5.Rdata")))
Dat <- local(get(load("/home/tim/Dropbox/RiifevanRaalteBijlsma/Data/resultsP bootsize999.rdata")))
names(Dat[[1]])
dimnames(Dat[[1]]$Female$Surf)
# get gy, need data object with columns for 
# Variable, Cohort, TTD, Prevalance
setwd("/home/tim/git/HLETTD")
hm <- local(get(load("Data/gyALL.Rdata")))
head(hm)
X <- Dat[[1]][["Female"]]$Surf
library(reshape2)
length(unique(hm$Morbidity))
PrevTTD <- do.call(rbind,lapply(Dat, function(X){
					
					require(reshape2)
					prevf <- melt(
							apply(X$Female$Surf,3,rowMeans,na.rm=TRUE), 
							varnames = c("TTD","Cohort"), value.name = "Prev")
					prevf$Morbidity <- X$var
					prevf$Sex <- "f"
					
					prevm <- melt(
							apply(X$Male$Surf,3,rowMeans,na.rm=TRUE), 
							varnames = c("TTD","Cohort"), value.name = "Prev")
					prevm$Morbidity <- X$var
					prevm$Sex <- "m"

					prev <- rbind(prevf,prevm)
					
					prev[,c("Morbidity","Cohort","Sex","TTD","Prev")]
				}))

PrevTTD[PrevTTD$Prev < 0 & !is.na(PrevTTD$Prev),]
length(unique(PrevTTD$Morbidity))
rownames(PrevTTD) <- NULL
head(PrevTTD)
save(PrevTTD, file = "/home/tim/git/HLETTD/Data/gyALL_P.Rdata")


do.old <- FALSE
if(do.old){
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

save(PrevTTD, file = "/home/tim/git/HLETTD/Data/gyALL.Rdata")
}

