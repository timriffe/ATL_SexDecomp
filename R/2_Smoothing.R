
# Author: tim
###############################################################################
# Instructions: 
# tend to steps 1-3 as appropriate, explained in situ.
# the script will produce a single list data object
# containing 3d arrays for 78 variables. (age,ttd,single-year cohorts), sexes separate.
# the object created at the end is called boot.results.list, saving out to ResultsP.Rdata.
# ------------------

# ------------------------
# 1) set parameters
nboot        <- 999    # ? how many should we do? 999?
do.this      <- TRUE    # change this to TRUE
make.figs    <- TRUE     # shall we make the summary historgrams?

# 2) set working directory: you'll need to modify the working dir. Possibly by generalizing the below
# code or else commenting it out altogether and setting manually. Up to you

# for TR computers...otherwise you need to set wd yourself

# ------------------------
if (system("hostname",intern=TRUE) %in% c("triffe-N80Vm","tim-ThinkPad-L440")){
	# if I'm on the laptop
	setwd("/home/tim/git/HLETTD")
	Cores <- 1 # laptop overheats...
	if (system("hostname",intern=TRUE) %in% c("tim-ThinkPad-L440")){
		Cores <- 4 # not used
	}
} else {
	if (system("hostname",intern=TRUE) == "PC-403478"){
		# on MPIDR PC
		setwd("U://git//HLETTD")
		Cores <- detectCores() # not used
	} 
}

Cores <- 32 # (you can also set this manually, override the above)

# 3) last thing to check:
# load in data, either change this path (commenting this one out so I keep it)
# or else make sure Data_long.Rdata is in a folder called Data inside the working directory...
#Dat        <- local(get(load("Data/Data_long.Rdata")))
Dat        <- local(get(load("Data/Data_longP.Rdata")))
# the rest should work fine without further ado. Two figures will also be created 
# in the working directory. These can be examined or thrown out.

# ------------------------
# get packages loaded
# ------------------------
# gives mclapply()-like functionality in Windows. 
# I usually would use parallel package, but that won't cut it on Windows
if (.Platform$OS.type == "windows"){
	# does Hydra have devools? I hope so.
	if (!"parallelsugar" %in% rownames(installed.packages())){
		library(devtools)
		install_github('nathanvan/parallelsugar')
	}
	library(parallelsugar)
} else {
	# this is for all else (Tim)
	library(parallel)
}

library(splines)
library(data.table)
library(reshape2)
library(lattice)
# ------------------------
# the following functions are defined here to 
# avoid having to source and set another path.

# ------------------------
source("R/apct.boot.R")

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

# TR: I think this list is already solidified, but just in case.
varnames <- local(get(load("Data/varnames_fit.Rdata"))) 

# -----------------------
# this is the slow part!
# -----------------------
# TR: this is toggled at the head of the script.
# just to make sure resources not too tied up
if (do.this){
	boot.results.list <- mclapply(varnames, function(varname, Dat, nboot){
				fem 	<- apct.boot.wrapper(
							Dat[Dat$sex == "f", ], 
							varname = varname, 
							sex = "f", 
							nboot = nboot)
				mal 	<- apct.boot.wrapper(
							Dat[Dat$sex == "m", ], 
							varname = varname, 
							sex = "m", 
							nboot = nboot)
				
			
                out 	<- list(Female=fem, Male = mal)
				out$var <- varname
				out
			}, Dat = Dat,  
			   nboot = nboot,    # nboot is set at the head of the script.
			   mc.cores = 32)   # ncores is set just above here
}

# Maarten: change path if necessary
save(boot.results.list, file = "Data/resultsP.rdata")

# stop here! rest of this is old.


do.old <- FALSE
if (do.old){
ResultsLong         <- do.call(rbind, boot.results.list)

ResultsLong$Cohort  <- ResultsLong$b_yr - ResultsLong$b_yr %% 5

# take simple means within ca, ta, Cohort, var, Sex
ResultsLong   		<- data.table(ResultsLong)

# choosing median here, could just as well be mean.
ResultsLong5  		<- ResultsLong[, list(pi = mean(median, na.rm=TRUE)), by = list(var, Sex, Cohort, ca, ta )]
ResultsLong5  		<- as.data.frame(ResultsLong5)
ResultsLong5L 		<- split(ResultsLong5,with(ResultsLong5, list(var, Sex, Cohort)))

# do a quick comparison with ThanoEmpirical results:

# function used, given a thano x chrono surface:
get_r <- function(A){
	c(T  = abs(cor(A$pi, A$ta, use = "complete.obs")), 
			A  = abs(cor(A$pi, A$ca, use = "complete.obs")),
			L  = abs(cor(A$pi, A$ta + A$ca, use = "complete.obs")),
			M  = abs(cor(A$pi, A$ca - A$ta, use = "complete.obs"))
	)
}

# --------------------------------------------
# This gets the correlations for each dim, cohort, sex, and span.
Results_r <- do.call(rbind,
		lapply(ResultsLong5L, function(X){
					out     <- X[1:4, ]
					out$Dim <- c("T","A","L","M")
					out$r   <- get_r(X)
					out$pi  <- NULL
					out$ca  <- NULL
					out$ta  <- NULL
					out
				})
)

comparison.tim <- FALSE
if (comparison.tim){
# compare these with the original results
Results_loess_r   	<- local(get(load("Data/Correlations.Rdata")))

# match ordering
Results_r       	<- Results_r[with(Results_r,order(var, Sex, Cohort, Dim)), ]
Results_loess_r   	<- Results_loess_r[with(Results_loess_r,order(var, sex, Cohort, Dim)), ]
Results_r$var   	<- gsub("_", "",Results_r$var)

Results_loess_r.5 	<- Results_loess_r[Results_loess_r$span == "0.5", ]
Results_loess_r.7 	<- Results_loess_r[Results_loess_r$span == "0.7", ]
Results_loess_r.9 	<- Results_loess_r[Results_loess_r$span == "0.9", ]
}
# match cohorts
Results_r       	<- Results_r[Results_r$Cohort %in% unique(Results_loess_r.5$Cohort), ]

# TR: hmmm. positive correlation, but that cloud is fatter than
# I'd like!
#plot(Results_r$r, Results_loess_r.5$r)
#plot(Results_r$r[Results_r$Cohort == 1915], Results_loess_r.5$r[Results_loess_r.5$Cohort == 1915])

# so let's start be redoing Figure 5a, 5b too see if they change much

Hist7      <- Results_r[Results_r$Cohort == 1915, ]
Hist7$rbin <- round(Hist7$r * 100) %/% 10 * 10
Hist7$rbin <- as.factor(Hist7$rbin)
Hist7$R    <- Hist7$r * 100
Hist7$Dim  <- as.factor(Hist7$Dim )

if (make.figs){
# Figures that ought to be mostly comparable with
# Figure 5 in ThanoEmpirical paper. Hopefully
# these histograms have a similar shape!!

pdf("Figure5a_v3.pdf", width = 3, height = 7)
histogram(~R | Dim, 
		data = Hist7[Hist7$Sex == "f", ], 
		col = gray(.4), 
		par.settings = list(strip.background = list(col = gray(.9))), 
		layout = c(1,4), 
		type = "count",
		breaks = seq(0,100, by = 10),
		index.cond = list(c(3, 4, 1, 2)), 
		xlab = list(label = "correlation coef * 100"),
		ylab = list(label = "variable count"),
		ylim = c(0, 47))
dev.off()

pdf("Figure5b_v3.pdf", width = 3, height = 7)
histogram(~R | Dim, 
		data = Hist7[Hist7$Sex == "m", ], 
		col = gray(.4), 
		par.settings = list(strip.background = list(col = gray(.9))), 
		layout = c(1, 4), 
		type = "count",
		breaks = seq(0, 100, by = 10),
		index.cond = list(c(3, 4, 1, 2)), 
		xlab = list(label = "correlation coef * 100"),
		ylab = list(label = "variable count"),
		ylim = c(0, 47))
dev.off()

} # end fig chunk


# TR will want a copy of this output.
save(ResultsLong, file = "ResultsLongBoot.Rdata")
save(Results_r, file = "CorrelationResultsBoot.Rdata")

getwd()

}
