# TR: code to regenerate results of ThanoEmpirical,
#     based on natural splines glm code from Maarten Bijlsma
#     further diagnostics needed on this method.

if (system("hostname",intern=TRUE) %in% c("triffe-N80Vm","tim-ThinkPad-L440")){
	# if I'm on the laptop
	setwd("/home/tim/git/ATL_SexDecomp/ATL_SexDecomp")
} else {
	# in that case I'm on Berkeley system, and other people in the dept can run this too
	setwd(paste0("/data/commons/",system("whoami",intern=TRUE),"/git/ATL_SexDecomp/ATL_SexDecomp"))
}



# TR: modified from MB's script. (under construction)
#setwd('U:/Collaboration/TR AVR')
#load("U:/Collaboration/TR AVR/Data/Data_long.Rdata")
Dat <- local(get(load("Data/Data_long.Rdata")))

library(splines)
library(data.table)
library(reshape2)
library(Epi)

varnames <- c("adl3_", 
		"adl5_", "iadl3_", "iadl5_", "cesd",  "lim_work", "srh", 
		"bmi", "back", "hosp", "hosp_stays", "hosp_nights", "nh", 
		"nh_stays", "nh_nights", "nh_now", "nh_mo", "nh_yr", "nh_days", 
		"doc", "doc_visits", "hhc", "meds", "surg", "dent", "shf", "adl_walk", 
		"adl_dress", "adl_bath", "adl_eat", "adl_bed", "adl_toilet", 
		"iadl_map", "iadl_tel", "iadl_money", "iadl_meds", "iadl_shop", 
		"iadl_meals", "mob", "lg_mus", "gross_mot", "fine_mot", "bp", 
		"diab", "cancer", "lung", "heart", "stroke", "psych", "arth", 
		"cc", "alc_ev", "alc_days", "alc_drinks", "smoke_ev", "smoke_cur", 
		"cesd_depr", "cesd_eff", "cesd_sleep", "cesd_happy", "cesd_lone", 
		"cesd_sad", "cesd_going", "cesd_enjoy", "prob75yo", "alz", "dem", 
		"srm", "pastmem", "ss", "c20b", "name_mo", 
		"name_dmo", "name_yr", "name_dwk", "name_sci", "name_cac", "name_pres", 
		"name_vp", "vocab", "tm", "med_exp", "dwr","twr","iwr",
		"iadl_calc", "mprob", "mprobev", "med_explog") 


# pare down to columns available in current iteration
varnames <- varnames[varnames %in% colnames(Dat)]

# cut data down, define time vars

# no ages under 65
#oops age is in months!
Dat      <- Dat[(Dat$age / 12) >= 65, ]

# birth year must be known
Dat      <- Dat[!is.na(Dat$b_yr), ]

# group to quinquennial
Dat$Coh5 <- Dat$b_yr -  Dat$b_yr %% 5 

# Coh5keep are cohorts kept for fitting
Coh5keep <- c(1900, 1905, 1910, 1915, 1920, 1925, 1930)

# Cog5 are cohorts for predicting
Coh5     <- c(1905, 1910, 1915, 1920, 1925) # i.e. we use the preceding and subsequent cohorts for help fitting

# select only cohorts used for fitting
Dat      <- Dat[Dat$Coh5 %in% Coh5keep, ]

# a final check
all(varnames %in% colnames(Dat))

## integer age representations
Dat$ta_int <- floor(Dat$ta)
Dat$ca_int <- floor(Dat$ca)
Dat$la_int <- floor(Dat$ta + Dat$ca)


# let's remove people with ta = -1
Dat <- Dat[Dat$ta >= 0,]


expit <- function(x){
	exp(x) / (exp(x) + 1)
}
cutla <- function(newdata, year1 = 1992, year2 = 2011, maxl = 100){
	# cut age
	newdata$la <- newdata$ca + newdata$ta
	mini       <- newdata$ca >= (year1 - newdata$b_yr - 1)
	maxi       <- newdata$la < (year2 - newdata$b_yr)
	newdata    <- newdata[mini & maxi, ]
	newdata    <- newdata[newdata$la < maxl, ]
	newdata
}

# TR: work just with this function, have different specifications of ns
# happen inside the same function and save to new output columns.
# consider C + (TAL), with C intercept and two of TAL detrended.
# CTL, or CAL also (currently have CAT). Remember for detrending,
# ns() needs to happen outside the call

fitns <- function(
		varname, 
		Dat, 
		sex = "m",
		t.age = 0:12,
		c.age = 70:100) {
	
	# compose call, explicit knots, artisanally chosen
#	this.call <- paste(varname, "~ 
#ns(b_yr, knots = seq(1902.5,1925.5,by=5)) + 
#ns(ta, knots = c(.5,1,2,4,7.5,10)) +  
#ns(ca, knots = seq(72.5,97.5,by=5))")

    # sex could jsut as easily be filtered before the call...
	Dat  <- Dat[Dat$sex == sex, ]
	
# Tr: include splines in call to make prediction easier.
#	Mcoh <- ns(Dat$b_yr, df = nk + 1)
#	Mta  <- ns(Dat$ta, knots = c(.5,1,2,4,7.5,10))
#	Mca  <- ns(ca, knots = seq(72.5,97.5,by=5))
	
	
	if(max(Dat[, varname], na.rm = TRUE) == 1) {
		# where possible, respect bounded data
		fit        <- glm(Dat[, varname]~ 
						      ns(b_yr, knots = seq(1902.5,1925.5,by=5)) + 
						      ns(ta, knots = c(.5,1,2,4,7.5,10)) +  
						      ns(ca, knots = seq(72.5,97.5,by=5)), 
				          data = Dat, 
						  family = quasibinomial) # removes warning message
				  
	} else {
		fit        <- glm(Dat[, varname]~ 
						ns(b_yr, knots = seq(1902.5,1925.5,by=5)) + 
						ns(ta, knots = c(.5,1,2,4,7.5,10)) +  
						ns(ca, knots = seq(72.5,97.5,by=5)), data = Dat)
	}
	
	# data for prediction, grid
	newdata    <- expand.grid(ta = t.age+.5, 
			                  ca = c.age+.5, 
							  b_yr = min(Dat$b_yr):max(Dat$b_yr))
					  
	# remove extrapolation points for glm prediction
	out        <- cutla(newdata, year1 = 1992, year2 = 2011)
	
	# easier to keep dimensions straight if we predict over rectangular grid, 
	# then throw out values outside range
	out$pi     <- predict(fit, out,type='response')
	
	loessfit <- loess(Dat[, varname] ~ b_yr + ta + ca,
			data = Dat, 
			weights = p_wt2, # this, plus point density both act as weights
			span = 0.5,     # a variable passed in, or smoothness
			# is similiar conceptually to a 1:1:1 aspect ratio. 
			# Everything is in years...
			normalize = FALSE,
			control = loess.control(trace.hat = "approximate"))
	
    # newdata must be a complete grid here, per expand.grid(),
    # otherwise stuff gets out of order.
    lpred       <- melt(predict(loessfit, newdata))
	
	# note that it messes up margin names
	lpred$ta    <- as.numeric(gsub("ta=","",gsub("ta= ","",lpred$ta)))
	lpred$ca    <- as.numeric(gsub("ca=","",gsub("ca= ","",lpred$ca)))
	lpred$b_yr  <- as.numeric(gsub("b_yr=","",lpred$b_yr))
	
	# now cut down to size
	lpred       <- cutla(lpred)
	
	# match ordering
	lpred       <- lpred[with(lpred,order(b_yr,ca,ta)), ]
	out         <- out[with(out,order(b_yr,ca,ta)), ]
	
	# append loess results
	out$piloess <- lpred$value
	
	# return results
    out
}
nsResults <- lapply(varnames, function(varname, Dat){
	Male        <- fitns(varname, Dat, sex = "m")	
	Female      <- fitns(varname, Dat, sex = "f")	
	# add meta vars
	Male$Sex    <- "m"
    Female$Sex  <- "f"
    out         <- rbind(Male, Female)
	out$var     <- varname
	out$ca      <- floor(out$ca)
	out$ta      <- floor(out$ta)
	out$la      <- NULL
	out
}, Dat = Dat)
# needed to generate results for Nov 16 Lab Talk
save(nsResults,file="/home/tim/git/APCTapps/APCTapps/LabTalk/Data/nsResults.Rdata")
#

res <- fitns("adl3_", Dat, "m")
head(res)

# where prediction is a df, as returned by fitns()
pred <- res
# compression trajectory function
CT <- function(pred,value.var = "pi",L = 85, C = 1915, keep = 10){
	ind <- with(pred, la == L & b_yr == C)
	pi <- pred[ind,value.var]
	ta <- pred[ind,"ta"]
	pi <- pi[order(ta)]
	ta <- sort(ta)
	if (length(pi) > keep){
		pi <- pi[1:keep]
		ta <- ta[1:keep]
	}
	list(x = ta, y = pi)
}

library(RColorBrewer)

mycols <- colorRampPalette(brewer.pal(9,"Blues"),space="Lab")(10)

plot(CT(res, L = 95, C = 1915), type = 'l', ylim = c(0,1), xlab = "TTD", ylab = "pi")
lines(CT(res, L = 95, C = 1910))
lines(CT(res, L = 95, C = 1905))

plot(CT(res, L = 85, C = 1915), type = 'l', ylim = c(0,1), xlab = "TTD", ylab = "pi")
lines(CT(res, L = 85, C = 1910))
lines(CT(res, L = 85, C = 1905))
## continue on to iteration. implement bootstrapping.
#nsResults <- lapply(varnames, function(varname, Dat){
#	Male        <- fitns(varname, Dat, sex = "m")	
#	Female      <- fitns(varname, Dat, sex = "f")	
#	# add meta vars
#	Male$Sex    <- "m"
#    Female$Sex  <- "f"
#    out         <- rbind(Male, Female)
#	out$var     <- varname
#	out$ca      <- floor(out$ca)
#	out$ta      <- floor(out$ta)
#	out$la      <- NULL
#	out
#		}, Dat = Dat)
##


#SurfA <- acast(Surfi[Surfi$Sex == "m", ], ta~ca~b_yr,value.var = "pi")
#dim(SurfA)
#range(SurfA, na.rm=TRUE)
#ticks <- seq(0,.75,by=.05)
#dimnames(SurfA)[[3]]
#for (i in dimnames(SurfA)[[3]]){
#	SurfMap(SurfA[,,i],ticks=ticks,bg=TRUE)
#	title(i)
#	locator(1)
#}
#SurfMap(SurfA[,,"1912"],ticks=ticks,bg=TRUE)
#title(1912)
#dev.new()
#SurfMap(SurfA[,,"1917"],ticks=ticks,bg=TRUE)
#title(1917)
#dev.new()
#SurfMap(SurfA[,,"1922"],ticks=ticks,bg=TRUE)
#title(1922)

#nsResultsLong  <- do.call(rbind, nsResults)
# ----------------------------------------------------
# try a more APC-style model, with term for L
# ----------------------------------------------------
# save this object for ThanoEmpirical replication, if desired (at end)
fitns2 <- function(
		varname, 
		Dat, 
		sex = "m",
		t.age = 0:12,
		c.age = 70:100) {
	# sex could just as easily be filtered before the call...
	Dat            <- Dat[Dat$sex == sex, ]
	# compose call, explicit knots, artisanally chosen

	newdata        <- expand.grid(ta = 0:12+.5, ca = 70:100+.5, b_yr = min(Dat$b_yr):max(Dat$b_yr))
	newdata$L      <- newdata$ta + newdata$ca
# remove extrapolation points
	newdata        <- cutla(newdata, year1 = 1992, year2 = 2011)
	newdata$la     <- NULL
	Dat$L          <- Dat$ca + Dat$ta

	this.call <- paste(varname, "~ ns(b_yr, knots = seq(1902.5,1925.5,by=5)) + 
					ns(ta, knots = c(.5,1,2,4,7.5,10)) +  
					ns(L, knots = seq(70,100,by=5))")
	
	if(max(Dat[, varname], na.rm = TRUE) == 1) {
		# where possible, respect bounded data
		fit        <- glm(this.call, 
				          data = Dat, 
				family = binomial) # removes warning message
		
	} else {
		fit        <- glm(this.call, data = Dat)
	}

	newdata$pi     <- predict(fit, newdata, type ="response")
	newdata$L      <- NULL
	newdata
}

# continue on to iteration. implement bootstrapping.
nsResults2 <- lapply(varnames, function(varname, Dat){
			cat(varname,"\n")
			Male        <- fitns2(varname, Dat, sex = "m")	
			Female      <- fitns2(varname, Dat, sex = "f")	
			# add meta vars
			Male$Sex    <- "m"
			Female$Sex  <- "f"
			out         <- rbind(Male, Female)
			out$var     <- varname
			out$ca      <- floor(out$ca)
			out$ta      <- floor(out$ta)
			out$la      <- NULL
			out
		}, Dat = Dat)



#source("R/SurfMap.R")
Loess <- local(get(load("/home/tim/git/ThanoEmpirical/ThanoEmpirical/Data/LoessQuinquenal.Rdata")))
sum(Loess[["adl3__0.7"]]$Male$Surf[,,"1915"], na.rm=TRUE) 
Surfi  <- nsResults[[1]]

Loess[["adl3__0.7"]]$Male$Surf["5","80","1915"]
SurfA["5","80","1914"]
#Surfii <- nsResults2[[1]]
SurfA  <- acast(Surfi[Surfi$Sex == "m", ], ta~ca~b_yr,value.var = "pi")
#SurfB  <- acast(Surfii[Surfii$Sex == "m", ], ta~ca~b_yr,value.var = "pi")

sum(SurfA[,,"1914"], na.rm=TRUE)

sum(SurfA, na.rm=TRUE)
sum(SurfB, na.rm=TRUE)

dim(SurfA)
range(SurfA, na.rm=TRUE)
range(SurfB, na.rm=TRUE)
ticks <- seq(0,.8,by=.05)
for (i in dimnames(SurfA)[[3]]){
	par(mfrow=c(2,1))
	SurfMap(SurfA[,,i],ticks=ticks,bg=TRUE)
	title(i,"no L")
	SurfMap(SurfB[,,i],ticks=ticks,bg=TRUE)
	title(i,"with L")
	locator(1)
}

par(mfrow=c(3,1))
SurfMap(Loess[["adl3__0.7"]]$Male$Surf[,,"1915"],ticks=ticks,bg=TRUE)
title("Loess, 1915-1919")
SurfMap(SurfA[,,"1917"],ticks=ticks,bg=TRUE)
title(i,"no L")
SurfMap(SurfB[,,"1917"],ticks=ticks,bg=TRUE)
title(i,"with L")



replicateThanoEmpirical <- FALSE
if (replicateThanoEmpirical){

	# choose which model to compare with
	#nsResultsLong <- do.call(rbind, nsResults)
	nsResultsLong <- do.call(rbind, lapply(nsResults2, function(X){
						colnames(X)[1] <- "pi"
                        X
						
					} ))
	
nsResultsLong$Cohort <- nsResultsLong$b_yr - nsResultsLong$b_yr %% 5

# take simple means within ca, ta, Cohort, var, Sex
nsResultsLong   	<- data.table(nsResultsLong)
nsResultsLong5  	<- nsResultsLong[, list(pi5 = mean(pi, na.rm=TRUE)), by = list(var, Sex, Cohort, ca, ta )]
nsResultsLong5  	<- as.data.frame(nsResultsLong5)
nsResultsLong5L 	<- split(nsResultsLong5,with(nsResultsLong5, list(var, Sex, Cohort)))

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
nsResults_r <- do.call(rbind,
		lapply(nsResultsLong5L, function(X){
					out     <- X[1:4, ]
					out$Dim <- c("T","A","L","M")
					out$r   <- get_r(X)
					out$pi  <- NULL
					out$ca  <- NULL
					out$ta  <- NULL
					out$pi5  <- NULL
					out
				})
)
#dim(nsResults_r)
#length(unique(nsResults_r$Cohort))
#length(varnames) * 2 * 4 * 7


# compare these with the original results
Results_loess_r <- local(get(load("/home/tim/git/ThanoEmpirical/ThanoEmpirical/Data/Correlations.Rdata")))

head(Results_loess_r)

# match ordering
nsResults_r$var <- gsub("_", "",nsResults_r$var)
nsResults_r     <- nsResults_r[with(nsResults_r,order(var, Sex, Cohort, Dim)), ]
Results_loess_r <- Results_loess_r[with(Results_loess_r,order(var, sex, Cohort, Dim)), ]



head(Results_loess_r)
head(nsResults_r)

Results_loess_r.5 <- Results_loess_r[Results_loess_r$span == "0.5", ]
Results_loess_r.7 <- Results_loess_r[Results_loess_r$span == "0.7", ]
Results_loess_r.9 <- Results_loess_r[Results_loess_r$span == "0.9", ]

# match cohorts
nsResults_r <- nsResults_r[nsResults_r$Cohort %in% unique(Results_loess_r.5$Cohort), ]

head(nsResults_r)
head(Results_loess_r.5)


plot(nsResults_r$r, Results_loess_r.5$r)
abline(a=0,b=1)
plot(nsResults_r$r[nsResults_r$Cohort == 1915], Results_loess_r.7$r[Results_loess_r.7$Cohort == 1915])



Hist7      <- nsResults_r[nsResults_r$Cohort == 1915, ]
Hist7$rbin <- round(Hist7$r * 100) %/% 10 * 10
Hist7$rbin <- as.factor(Hist7$rbin)
Hist7$R    <- Hist7$r * 100
Hist7$Dim  <- as.factor(Hist7$Dim )

#
head(Hist7)
library(lattice)
pdf("Figures/ThanoEmpirical_Figure5a_v2.pdf", width = 3, height = 7)
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

pdf("Figures/ThanoEmpirical_Figure5b_v2.pdf", width = 3, height = 7)
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

getwd()



# redo ThanoEMpirical paper figures
head(nsResultsLong5)

Surf <- acast(nsResultsLong5[with(nsResultsLong5,var == "psych" & Cohort == 1915 & Sex == "m"),], 
		ta~ca, value.var = "pi5")
# 
pdf("Figures/ThanoEmpirical_Figure2a.pdf", width = 10, height = 6)
#dev.new(width = 10, height = 6)
SurfMap(Surf,
		napprox = 9,
		contour = TRUE,
		outline = FALSE,
		bg = TRUE,
		xlab = "chronological age", 
		ylab = "thanatological age")
dev.off()


# ----------------------------
# A: chronological age
Surf <- acast(nsResultsLong5[with(nsResultsLong5,var == "back" & Cohort == 1915 & Sex == "f"),], 
		ta~ca, value.var = "pi5")
pdf("Figures/ThanoEmpirical_Figure2b.pdf", width = 10, height = 6)
#dev.new(width = 10, height = 6)
SurfMap(Surf,
		napprox = 9,
		contour = TRUE,
		outline = FALSE,
		bg = TRUE,
		xlab = "chronological age", 
		ylab = "thanatological age")
dev.off()

# ----------------------------
# L: lifespan:
Surf <- acast(nsResultsLong5[with(nsResultsLong5,var == "smoke_ev" & Cohort == 1915 & Sex == "f"),], 
		ta~ca, value.var = "pi5")
pdf("Figures/ThanoEmpirical_Figure4a.pdf", width = 10, height = 6)
#dev.new(width = 10, height = 6)
SurfMap(Surf,
		napprox = 9,
		contour = TRUE,
		outline = FALSE,
		bg = TRUE,
		xlab = "chronological age", 
		ylab = "thanatological age")
dev.off()

#------------------------------------------------
# M: function of both
Surf <- acast(nsResultsLong5[with(nsResultsLong5,var == "bp" & Cohort == 1915 & Sex == "m"),], 
		ta~ca, value.var = "pi5")
pdf("Figures/ThanoEmpirical_Figure4b.pdf", width = 10, height = 6)
#dev.new(width = 10, height = 6)
SurfMap(Surf,
		napprox = 9,
		contour = TRUE,
		outline = FALSE,
		bg = TRUE,
		xlab = "chronological age", 
		ylab = "thanatological age")
dev.off()
}
# --------------------------------
# now make a flip-book


