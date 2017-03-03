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

do.this <- FALSE
if (do.this){

boot.default <- lapply(varscompare, function(varname, Dat){
			get.booty(apct.boot(Dat, varname = varname, nboot = 1000,
		              YearsFromEdge = 0, MagFactor = 1))
		}, Dat = Dat)
boot.2.2 <- lapply(varscompare, function(varname, Dat){
			get.booty(apct.boot(Dat, varname = varname, nboot = 1000,
							YearsFromEdge = 2, MagFactor = 2))
		}, Dat = Dat)
boot.5.2 <- lapply(varscompare, function(varname, Dat){
			get.booty(apct.boot(Dat, varname = varname, nboot = 1000,
							YearsFromEdge = 5, MagFactor = 2))
		}, Dat = Dat)
boot.2.5 <- lapply(varscompare, function(varname, Dat){
			get.booty(apct.boot(Dat, varname = varname, nboot = 1000,
							YearsFromEdge = 2, MagFactor = 5))
		}, Dat = Dat)
boot.5.5 <- lapply(varscompare, function(varname, Dat){
			get.booty(apct.boot(Dat, varname = varname, nboot = 1000,
							YearsFromEdge = 5, MagFactor = 5))
		}, Dat = Dat)
boot.2.10 <- lapply(varscompare, function(varname, Dat){
			get.booty(apct.boot(Dat, varname = varname, nboot = 1000,
							YearsFromEdge = 2, MagFactor = 10))
		}, Dat = Dat)
boot.5.10 <- lapply(varscompare, function(varname, Dat){
			get.booty(apct.boot(Dat, varname = varname, nboot = 1000,
							YearsFromEdge = 5, MagFactor = 10))
		}, Dat = Dat)

# so this will be a set of named lists. matrices be called out using
# get.mat(goods, column = "median", cohort = 1915)
# where possible columns are "median", "mean", "2.5%", "97.5%", "Diff"
# I've not simplified plotting, but that can happen.
getwd()
save(boot.default, file = "Data/BootDiagnostics/boot.default.Rdata")
save(boot.2.2, file = "Data/BootDiagnostics/boot.2.2.Rdata")
save(boot.5.2, file = "Data/BootDiagnostics/boot.5.2.Rdata")
save(boot.2.5, file = "Data/BootDiagnostics/boot.2.5.Rdata")
save(boot.5.5, file = "Data/BootDiagnostics/boot.5.5.Rdata")
save(boot.2.10, file = "Data/BootDiagnostics/boot.2.10.Rdata")
save(boot.5.10, file = "Data/BootDiagnostics/boot.5.10.Rdata")
}
names(boot.default) <- varscompare
names(boot.2.2) <- varscompare
names(boot.5.2) <- varscompare
names(boot.2.5) <- varscompare
names(boot.5.5) <- varscompare
names(boot.2.10) <- varscompare
names(boot.5.10) <- varscompare

fit         	<- glm(back ~ 
				ns(b_yr, knots = seq(1902.5,1925.5,by=5)) + 
				ns(ta, knots = c(.5,1,2,4,7.5,10)) +  
				ns(ca, knots = seq(72.5,97.5,by=5)), 
		weights = p_wt2,
		data = Dat,
		family = quasibinomial)
out <- boot.default[["back"]][,1:4]
compare <- cbind(out, compare= predict(fit, out,type='response'))




head(boot.default[[1]])
compare    <- get.mat(compare, "compare", 1915)
def    <- get.mat(boot.default[["back"]], "median", 1915)
mag2.5 <- get.mat(boot.2.5[["back"]], "median", 1915)
mag2.10 <- get.mat(boot.2.10[["back"]], "median", 1915)
mag5.10 <- get.mat(boot.5.10[["back"]], "median", 1915)
source("R/SurfMap.R")
library(RColorBrewer)
graphics.off()
SurfMap(def)
title("default")
dev.new()
SurfMap(mag5.10 )
title("mag")
Diff <- (mag5.10 - def)

Diff <- def - compare
range(Diff, na.rm=TRUE)
mag <- max(abs(range(Diff, na.rm=TRUE)))

breaks <- seq(-mag,mag,length=11)
cols   <- colorRampPalette(brewer.pal(9,"RdBu"),space="Lab")
image(t(Diff), breaks = breaks, col = rev(cols(length(breaks)-1)))
contour(t(Diff),breaks=breaks,add=TRUE)

image(t(Diff / def))
contour(t(Diff / def),add=TRUE)






