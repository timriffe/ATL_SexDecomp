# here we check to see if there is problematic sparseness that could
# lead to bias in areas where we estiamte.

if (system("hostname",intern=TRUE) %in% c("triffe-N80Vm","tim-ThinkPad-L440")){
	# if I'm on the laptop
	setwd("/home/tim/git/ATL_SexDecomp/ATL_SexDecomp")
} else {
	# in that case I'm on Berkeley system, and other people in the dept can run this too
	setwd(paste0("/data/commons/",system("whoami",intern=TRUE),"/git/ATL_SexDecomp/ATL_SexDecomp"))
}

# contains scripts to do this diagnostic.
source("R/apct.boot.R")
# load in data
Dat        <- local(get(load("Data/Data_long.Rdata")))
# no ages under 65
#oops age is in months!
Dat        <- Dat[(Dat$age / 12) >= 65, ]
# birth year must be known
Dat        <- Dat[!is.na(Dat$b_yr), ]
# group to quinquennial
Dat$Coh5   <- Dat$b_yr - Dat$b_yr %% 5 
# Coh5keep are cohorts kept for fitting
Coh5keep   <- seq(1900, 1930, by = 5) 
# select only cohorts used for fitting
Dat        <- Dat[Dat$Coh5 %in% Coh5keep, ]
## integer age representations
Dat$ta_int <- floor(Dat$ta)
Dat$ca_int <- floor(Dat$ca)
Dat$la_int <- floor(Dat$ta + Dat$ca)
# let's remove people with ta = -1
Dat        <- Dat[Dat$ta >= 0,]
# even tho most neg are very close to zero
# end data prep preamble
Dat$adl3_  <- ifelse(Dat$adl3_ > 0, 1, 0)

# -----------------------------------------

(varnames <- local(get(load("Data/varnames.Rdata"))))
varscompare <- c("back","adl3_","psych","lung")
#magfactors  1,2,5,10 
#YearsFromEdge  2,5 


library(parallel)
boot.default <- lapply(varscompare, function(varname, Dat){
			get.goods(apct.boot(Dat, varname = varname, nboot = 1000,
		              YearsFromEdge = 0, MagFactor = 1))
		}, Dat = Dat)
boot.2.2 <- lapply(varscompare, function(varname, Dat){
			get.goods(apct.boot(Dat, varname = varname, nboot = 1000,
							YearsFromEdge = 2, MagFactor = 2))
		}, Dat = Dat)
boot.5.2 <- lapply(varscompare, function(varname, Dat){
			get.goods(apct.boot(Dat, varname = varname, nboot = 1000,
							YearsFromEdge = 5, MagFactor = 2))
		}, Dat = Dat)
boot.2.5 <- lapply(varscompare, function(varname, Dat){
			get.goods(apct.boot(Dat, varname = varname, nboot = 1000,
							YearsFromEdge = 2, MagFactor = 5))
		}, Dat = Dat)
boot.5.5 <- lapply(varscompare, function(varname, Dat){
			get.goods(apct.boot(Dat, varname = varname, nboot = 1000,
							YearsFromEdge = 5, MagFactor = 5))
		}, Dat = Dat)
boot.2.10 <- lapply(varscompare, function(varname, Dat){
			get.goods(apct.boot(Dat, varname = varname, nboot = 1000,
							YearsFromEdge = 2, MagFactor = 10))
		}, Dat = Dat)
boot.5.10 <- lapply(varscompare, function(varname, Dat){
			get.goods(apct.boot(Dat, varname = varname, nboot = 1000,
							YearsFromEdge = 5, MagFactor = 10))
		}, Dat = Dat)

# so this will be a set of named lists. matrices be called out using
# get.mat(goods, column = "median", cohort = 1915)
# where possible columns are "median", "mean", "2.5%", "97.5%", "Diff"
# I've not simplified plotting, but that can happen.


