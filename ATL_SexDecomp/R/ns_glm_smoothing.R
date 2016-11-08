# CarstensenSmoothing.R

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
newdata <- function(newdata, year1 = 1992, year2 = 2011){
	# cut age
	newdata$la <- newdata$ca + newdata$ta
	mini       <- newdata$la > (year1 - newdata$b_yr - 1)
	maxi       <- newdata$la < (year2 - newdata$b_yr)
	newdata    <- newdata[mini & maxi, ]
	newdata
}

# ideally using the weights is done in a bootstrap
# type step, cause otherwise we use them twice
# and binomial hates them

fitns <- function(
		varname, 
		Dat, 
		sex = "m",
		t.age = 0:12,
		c.age = 70:100) {
	
	# compose call, explicit knots, artisanally chosen
	this.call <- paste(varname, "~ 
ns(b_yr, knots = seq(1902.5,1925.5,by=5)) + 
ns(ta, knots = c(.5,1,2,4,7.5,10)) +  
ns(ca, knots = seq(72.5,97.5,by=5))")

    # sex could jsut as easily be filtered before the call...
	Dat  <- Dat[Dat$sex == sex, ]
	
# Tr: include splines in call to make prediction easier.
#	Mcoh <- ns(Dat$b_yr, df = nk + 1)
#	Mta  <- ns(Dat$ta, knots = c(.5,1,2,4,7.5,10))
#	Mca  <- ns(ca, knots = seq(72.5,97.5,by=5))
	if(max(Dat[, varname], na.rm = TRUE) == 1) {
		# where possible, respect bounded data
		fit        <- glm(this.call, 
				          data = Dat, 
						  family = quasibinomial) # removes warning message
				  
	} else {
		fit        <- glm(this.call, data = Dat)
	}
	#pred           <- predict(fit, Dat, type='response')
	newdata        <- expand.grid(ta = t.age+.5, ca = c.age+.5, b_yr = min(Dat$b_yr):max(Dat$b_yr))
	# remove extrapolation points
	newdata        <- cutla(newdata, year1 = 1992, year2 = 2011)
	#pred           <- predict(fit, Dat, type='response')
	# easier to keep dimensions straight if we predict over rectangular grid, 
	# then throw out values outside range
	#newdata        <- newdata + .5
	logitpi        <- predict(fit, newdata)
	newdata$pi     <- expit(logitpi)
	newdata
}

# continue on to iteration. implement bootstrapping.


