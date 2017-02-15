

# here we check to see if there is problematic sparseness that could
# lead to bias in areas where we estiamte.

# Look here first! this IS the diagnostic parameter set
# define edge cases
# yippie!
YearsFromEdge       <- 3 # attention!
MagnificationFactor <- 5
#varnamechar         <- "adl3_"
varnamechar         <- "psych"
# Achtung!

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

head(Dat)



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

# TR: this line super important 
Dat <- Dat[order(Dat$id), ]


byr <- 1915
2011 - byr - 5
Dat$edgies    <- Dat$la_int > 2011 - Dat$b_yr - YearsFromEdge

# TR: new selection of first cases
get.first <- function(x){
	x       <- sort(x)
	lengths <- rle(x)$lengths
	cumsum(lengths) - lengths + 1
}

select.first <- get.first(Dat$id)
dataid       <- Dat[select.first, ]
# this rescales the draw weight for the edgies,
# and actually edgies should be determined prior to 
# defining dataid so that we don't do it twice.
dataid$drawweight                <- dataid$p_wt2
dataid$drawweight[dataid$edgies] <- dataid$drawweight[dataid$edgies] * MagnificationFactor

# assign person weight of first observation to each subsequent observation
Dat$firstweight                     <- rep(dataid$p_wt2,rle(Dat$id)$lengths)
# given repalcement selection by a draw weight, the remaining observations
# of individuals are thusly reweighted
Dat$rescaleweight                   <- Dat$p_wt2 / Dat$firstweight

# and for those that are also edgies we need to weight in the opposite direction
# same magnitude (but there will be more such people)
Dat$rescaleweight[Dat$edgies]       <- Dat$rescaleweight[Dat$edgies] / MagnificationFactor

# cut down dataid object to two needed columns, id and drawweight
dataid3 <- dataid[,c("id","drawweight")]

#

# copied and pasted
col.index.outc 	<- grep(paste0('^',varnamechar,'$'), colnames(Dat))
col.index.id 	<- grep(paste0('^','id','$'), colnames(Dat))
col.index.b_yr 	<- grep(paste0('^','b_yr','$'), colnames(Dat))
col.index.ta 	<- grep(paste0('^','ta','$'), colnames(Dat))
col.index.ca 	<- grep(paste0('^','ca','$'), colnames(Dat))
col.index.rscw 	<- grep(paste0('^','rescaleweight','$'), colnames(Dat))

Dat.3 <- Dat[,c(col.index.outc,col.index.id,col.index.rscw,
				col.index.b_yr, col.index.ta,col.index.ca)]

#expit <- function(x){
#	exp(x) / (exp(x) + 1)
#}
#cutla2 <- function(newdata, year1 = 1992, year2 = 2011, maxl = 100){
#	# cut age
#	newdata$la <- newdata$ca + newdata$ta
#	mini       <- newdata$ca >= (year1 - newdata$b_yr - 1)
#	maxi       <- newdata$la < (year2 - newdata$b_yr)
#	newdata    <- newdata[mini & maxi, ]
#	newdata    <- newdata[newdata$la < maxl, ]
#	newdata    <- newdata[newdata$ta <= 12, ]
#	newdata
#}

#  end default header here
# ----------------------------------------
# get 2x2 cells:
#unique(Dat$adl3_)
#Dat$adl3_ <- ifelse(Dat$adl3_ > 0, 1, 0)
#count1s <- function(x){
#	sum(x==1,na.rm=TRUE)
#}
#count0s <- function(x){
#	sum(x==0,na.rm=TRUE)
#}
#ones 			<- acast(Dat, tafloor2~cafloor2~Coh5, value.var = "adl3_",count1s)
#zeros 			<- acast(Dat, tafloor2~cafloor2~Coh5, value.var = "adl3_",count0s)
#
#ones 			<- melt(ones, varnames=c("ta","ca","b_yr"))
#zeros 			<- melt(zeros, varnames=c("ta","ca","b_yr"))
#ones        	<- cutla2(ones, year1 = 1992, year2 = 2011)
#zeros        	<- cutla2(zeros, year1 = 1992, year2 = 2011)
#
#ones$zeros 									<- zeros$value
#colnames(ones)[colnames(ones) == "value"] 	<- "ones"
#ones 										<- ones[,c("ta", "ca","la", "b_yr", "ones",  "zeros")]
#ones$N 										<- ones$ones + ones$zeros
#ones 										<- ones[order(ones$N),]
#
#ones


# ------------------------------------------------
# oversample la > year - b_yr - 5
#head(Dat)
#table(2011 - Dat$b_yr)
#Dat$edgies <- Dat$la_int > 2011 - Dat$b_yr - 5
#head(Dat[Dat$edgies & Dat$la >95,])


