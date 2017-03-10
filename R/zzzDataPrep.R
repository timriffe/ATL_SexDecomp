# Author: tim
###############################################################################

# for Tim, this will choke
if (system("hostname",intern=TRUE) %in% c("triffe-N80Vm","tim-ThinkPad-L440")){
	# if I'm on the laptop
	setwd("/home/tim/git/HLETTD")
} else {
	# in that case I'm on Berkeley system, and other people in the dept can run this too
	setwd(paste0("/data/commons/",system("whoami",intern=TRUE),"/git/HLETTD"))
}
getwd()
# install.packages("lubridate")
library(lubridate)
library(data.table)


# cleaning/processing functions
# need to make some codes usable...


# this code imported from ThanoEmpirical, but it is in need of some serious overhaul

convertYN <- function(x){
	xx                  <- rep(NA, length(x))
	xx[grepl("yes", x)] <- 1
	xx[grepl("no", x)]  <- 0
	invisible(xx)
}
convertCI <-  function(x){
	xx                  <- rep(NA, length(x))
	
	xx[x == "1.correct"] <- 0
	xx[x == "2.correct, 1st try"] <- 0
	xx[x == "1.correct, 2nd try"] <- .5
	xx[x == "0.incorrect"]  <- 1
	invisible(as.numeric(xx))
}
convertCESD <- function(x){
	xx <- rep(NA,length(x))
	xx[x == "0.no"]                    <- 0
	xx[x == "4. none or almost none"]  <- 0
	xx[x == "3. some of the time"]     <- .5
	xx[x == "2. most of the time" ]    <- .75
	xx[x ==  "1.yes"]                  <- 1
	xx[x ==  "1. all or almost all"]   <- 1
	xx
}
convertDates <- function(Dat){
	# can't be done with apply because we can't have Date class matrices...
	DateInd       <- grep(pattern="_dt",colnames(Dat))
	for (i in DateInd){
		Dat[,i]    <- as.Date(Dat[,i],origin="1960-1-1")
	}
	invisible(Dat)
}
getThanoAge <- function(Date, DeathDate){
	out <- rep(NA, length(Date))
	Ind <- !is.na(Date)
	out[Ind] <- lubridate::decimal_date(DeathDate[Ind]) - lubridate::decimal_date(Date[Ind])
	out
}
getChronoAge <- function(Date, BirthDate){
	out <- rep(NA, length(Date))
	Ind <- !is.na(Date) & !is.na(BirthDate)
	out[Ind] <- lubridate::decimal_date(Date[Ind]) - lubridate::decimal_date(BirthDate[Ind])
	out
}
imputeWeights <- function(wt,intv_dt){
	if (all(wt == 0)){
		wt2 <- NA * wt
		return(wt2)
	}
	if (sum(wt>0) == 1){
		wt2 <- approx(x = intv_dt[wt>0],
				y = wt[wt>0],
				xout = intv_dt,
				rule = 1:2,
				method = "constant",
				f = .5)$y
	}
	if (sum(wt>0)>=2){
		wt2 <- approx(x = intv_dt[wt>0],
				y = wt[wt>0],
				xout = intv_dt,
				rule = 1:2,
				method = "linear")$y 
	}
	return(wt2)
}

# converts to long format, assumes thano age columns already appended:
#

Dat         <- local(get(load("Data/thanos_long_v3_1.RData")))
nrow(Dat)

# remove missed interviews
Dat         <- Dat[!is.na(Dat$intv_dt), ]
nrow(Dat)
nrow(Dat)/ length(unique(Dat$id)) # avg interviews / id
# change all factors to character (to be later recoded in some instances)
#str(Dat)

#Dat[sapply(Dat, is.factor)] <- lapply(Dat[sapply(Dat, is.factor)], as.character)

# make sex column easier to use:
Dat$sex     <- ifelse(Dat$sex == "1.male","m","f")

# reduce to deceased-only

Dat         <- Dat[Dat$dead == 1, ]
nrow(Dat)
nrow(Dat)/ length(unique(Dat$id))

# convert dates to native R format
Dat         <- convertDates(Dat)

# merge weights:
Dat$nh_wt[is.na(Dat$nh_wt)] <- 0
Dat$p_wt <- Dat$p_wt + Dat$nh_wt

Dat <- data.table(Dat)
# take care of the rest: 

Dat <- Dat[,p_wt2 := imputeWeights(p_wt,intv_dt), by = list(id) ]
#Dat$p_wt2[Dat$p_wt==0]
# all zeros removed
# 2341 observations thrown as leading 0s, affecting 934 ids
# 3227 total observations thrown (including all-0s), 1361 total ids affected
Dat <- Dat[!is.na(Dat$p_wt2),]
nrow(Dat)
nrow(Dat)/ length(unique(Dat$id))
# calculate thanatological age
Dat$ta <- getThanoAge(Dat$intv_dt, Dat$d_dt_r) # TR: d_dt change to d_dt_r

Dat$ca <- getChronoAge(Dat$intv_dt, Dat$b_dt)
# there is one individual with an NA b_dt, and NA age,
# but thano age is known

# convert yes/no responses to 1,0
YNcols <- apply(Dat, 2, function(x){
			xx <- unique(x)
			length(xx) <= 4 & any(grepl("yes",xx))
		})
CIcols <- apply(Dat, 2, function(x){
			xx <- unique(x)
			length(xx) <= 5 & any(grepl("correct",xx))
		}) 

colnames(Dat)[YNcols]

colnames(Dat)[CIcols] 

Dat <- data.frame(Dat)
Dat[YNcols] <- lapply(Dat[YNcols], convertYN)
Dat[CIcols] <- lapply(Dat[CIcols], convertCI)
#################################3
# TR: stopped here
#################################


#head(Dat)

# remove lt, vig, ulc, too inconsistent
Dat$lt        <- NULL
Dat$vig       <- NULL
Dat$ulc       <- NULL
Dat$lt_freq   <- NULL
Dat$mod_freq  <- NULL
Dat$vig_freq  <- NULL
Dat$c86b      <- NULL # only in a couple waves
Dat$dem       <- NULL
Dat$alz       <- NULL
Dat$iadl_calc <- NULL
Dat$prob75yo  <- NULL
Dat$nh_mo     <- NULL
Dat$nh_yr     <- NULL
Dat$nh_days   <- NULL
### med expenditure needs to be removed, even though it has a very clear thano pattern
## mprobev / mprob need to go : too inconsistent
Dat$mprob 		<- NULL
Dat$mprobev 	<- NULL
Dat$med_explog 	<- NULL
Dat$med_exp 	<- NULL
# recode medical expenditure to actual values:

#1=0 to $1,000                 500
#2=~$1000                     1000
#3=$1,001 to 5,000            2500
#4=~$5,000                    5000
#5=$5,001 to $25,000         15000
#6=~$25,000                  25000
#7=$25,001 to $100,000       62500
#8=~$100,000                100000
#9=$100,001 to $500,000     300000
#10=~$500,000               500000
#11=$500,000+              1000000

# med exp thrown out, although it is a very clear pattern.
#rec.vec <- c(500,1000,2500,5000,15000,25000,62500,100000,300000,500000,1000000, NA,NA)
#names(rec.vec)      <- c("1 : 0 to 1000-",
#  "2 : about 1000",
#  "3 : 1001 to 5000-",
#  "4 : about 5000",
#  "5 : 5001 to 25000-",
#  "6 : about 25000",
#  "7 : 25001 to 100000-",
#  "8 : about 100000",
#  "9 : 100001 to 500000-",
#  "10: about 500000",
#  "11: 500000 above" ,
#  "NA" , "")
#
#Dat$med_exp         <- rec.vec[as.character(Dat$med_exp)]
#Dat$med_explog      <- log(Dat$med_exp )
# recode self reported health 1 = excellent - 5 = poor
srhrec              <- c(0:4,NA)
names(srhrec)       <- sort(unique(Dat$srh))
Dat$srh             <- srhrec[Dat$srh] / 4 # now all between 0 and 1. 1 worst.
names(srhrec)       <- sort(unique(Dat$srm))
Dat$srm             <- srhrec[Dat$srm] / 4 

# same, worse, better recode:  (1 bad, 0 good)
pastmem             <- c(0:2,NA)
names(pastmem)      <- sort(unique(Dat$pastmem))
Dat$pastmem         <- pastmem[Dat$pastmem] / 2

# do cesd questions (1 bad, 0 good)
cesdquestions       <- colnames(Dat)[grepl("cesd", colnames(Dat))]
cesdquestions       <- cesdquestions[cesdquestions != "cesd"]
Dat[cesdquestions]  <- lapply(Dat[cesdquestions],convertCESD)

# cesd_enjoy is flipped yet again, because 1 is 'yes I enjoyed life',
# and we want high = bad.
Dat$cesd_enjoy      <- 1 - Dat$cesd_enjoy
Dat$cesd_happy      <- 1 - Dat$cesd_happy
# create a single Total Word Recall variables, twr
#"tr20w"(waves(2-10),"tr40w" (waves1-2)
# i.e. 1 is the worst recall, and 0 is the best recall
Dat$tr20w                   <- 1 - Dat$tr20w / 20
Dat$tr40w                   <- 1 - Dat$tr40w / 40

NAind                       <- is.na(Dat$tr20w) & is.na(Dat$tr40w)
BothInd                     <- !is.na(Dat$tr20w) & !is.na(Dat$tr40w)
Dat$tr20w[is.na(Dat$tr20w)] <- 0
Dat$tr40w[is.na(Dat$tr40w)] <- 0
sum(BothInd) == 0 # (otherwise we'd need to divide these by two after adding)
Dat$twr                     <- Dat$tr20w + Dat$tr40w
Dat$twr[NAind]              <- NA

# vocab: 1 worst 0 best
Dat$vocab <- 1 - Dat$vocab / 10

# total mental: 1 worst, 0 best
Dat$tm    <- 1 - Dat$tm / 15

# delayed word recall
Dat$dr20w                   <- 1 - Dat$dr20w / 20
Dat$dr10w                   <- 1 - Dat$dr10w / 10

NAind                       <- is.na(Dat$dr20w) & is.na(Dat$dr10w)
BothInd                     <- !is.na(Dat$dr20w) & !is.na(Dat$dr10w)
Dat$dr20w[is.na(Dat$dr20w)] <- 0
Dat$dr10w[is.na(Dat$dr10w)] <- 0
sum(BothInd) == 0 # (otherwise we'd need to divide these by two after adding)
Dat$dwr                    <- Dat$dr20w + Dat$dr10w
Dat$dwr[NAind]              <- NA

# immediate word recall
Dat$ir20w                   <- 1 - Dat$ir20w / 20
Dat$ir10w                   <- 1 - Dat$ir10w / 10

NAind                       <- is.na(Dat$ir20w) & is.na(Dat$ir10w)
BothInd                     <- !is.na(Dat$ir20w) & !is.na(Dat$ir10w)
Dat$ir20w[is.na(Dat$ir20w)] <- 0
Dat$ir10w[is.na(Dat$ir10w)] <- 0
sum(BothInd) == 0 # (otherwise we'd need to divide these by two after adding)
Dat$iwr                     <- Dat$ir20w + Dat$ir10w
Dat$iwr[NAind]              <- NA

# memory problem:
#[1] ""                                "0. no"                          
#[3] "1. yes"                          "NA"                             
#[5] "4. disp prev record and no cond"
#mprob <- c(NA,0,1,0,NA)
#names(mprob) <- sort(unique(Dat$mprob))
#Dat$mprob <- mprob[Dat$mprob]
# vocab

# scale to fit btwn 0 and 1
rescale <- function(var,Dat,compelment = FALSE){
	Dat[[var]] <- Dat[[var]] / max(Dat[[var]], na.rm = TRUE)
	if (compelment){
		Dat[[var]] <- 1 - Dat[[var]]
	}
	Dat
}

Dat     <- rescale("mob", Dat, FALSE)
Dat     <- rescale("lg_mus", Dat, FALSE) 
Dat     <- rescale("gross_mot", Dat, FALSE)
Dat     <- rescale("fine_mot", Dat, FALSE)

Dat     <- rescale("ss", Dat, TRUE) # complement because more was better in original
Dat     <- rescale("cc", Dat, FALSE)
Dat     <- rescale("alc_days", Dat, FALSE)

Dat     <- rescale("adl3_", Dat, FALSE)
Dat     <- rescale("adl5_", Dat, FALSE)
Dat     <- rescale("iadl3_", Dat, FALSE)
Dat     <- rescale("iadl5_", Dat, FALSE)
Dat     <- rescale("cesd", Dat, FALSE)

# TODO: does high = bad for these variables?
# "bmi"        "back"       "dent"       "alc_ev"   
# "pastmem"    "dwr"        "twr"        "iwr" 

#checkwaves <- function(var,Dat){
#  table(Dat[[var]],Dat[["wave"]])
#}
#checkwaves("adl3_",Dat)
#checkwaves("adl5_",Dat)
#checkwaves("iadl3_",Dat)
#checkwaves("iadl5_",Dat)
#checkwaves("cesd",Dat)

# -------------------------------------------------------
# for binning purposes, akin to 'completed age'
Dat$tafloor <- floor(Dat$ta)
Dat$cafloor <- floor(Dat$ca)
# I guess actual interview date could be some weeks prior to registered
# interview date? There are two negative thano ages at wave 4 otherwise, but
# still rather close. Likely died shortly after interview.
Dat$tafloor[Dat$tafloor < 0] <- 0

Dat$cafloor2 <- Dat$cafloor - Dat$cafloor %% 2
Dat$tafloor2 <- Dat$tafloor - Dat$tafloor %% 2

Dat$cafloor3 <- Dat$cafloor - Dat$cafloor %% 3
Dat$tafloor3 <- Dat$tafloor - Dat$tafloor %% 3
# save out, so this doesn't need to be re-run every time
save(Dat,file = "Data/Data_long.Rdata")

#Dat <- local(get(load("Data/Data_long.Rdata")))

#apply(Dat[,colnames(Dat)%in%varnames],2,range,na.rm=TRUE)
############################################################################################
# BELOW is just to generate metadata for an appendix:
# these are the ones we keep:
varnames <- c("adl3_", "adl5_", "iadl3_", "iadl5_", "cesd", "lim_work", "srh", 
		"bmi", "back", "hosp", "hosp_stays", "hosp_nights", "nh", "nh_stays", 
		"nh_nights", "nh_now", "doc", "doc_visits", "hhc", "meds", "surg", 
		"dent", "shf", "adl_walk", "adl_dress", "adl_bath", "adl_eat", 
		"adl_bed", "adl_toilet", "iadl_map", "iadl_tel", "iadl_money", 
		"iadl_meds", "iadl_shop", "iadl_meals", "mob", "lg_mus", "gross_mot", 
		"fine_mot", "bp", "diab", "cancer", "lung", "heart", "stroke", 
		"psych", "arth", "cc", "alc_ev", "alc_days", "alc_drinks", "smoke_ev", 
		"smoke_cur", "cesd_depr", "cesd_eff", "cesd_sleep", "cesd_happy", 
		"cesd_lone", "cesd_sad", "cesd_going", "cesd_enjoy", "srm", "pastmem", 
		"ss", "c20b", "name_mo", "name_dmo", "name_yr", "name_dwk", "name_sci", 
		"name_cac", "name_pres", "name_vp", "vocab", "tm", "med_exp", 
		"dwr", "twr", "iwr", "mprob", "mprobev", "med_explog")


varnames <- varnames[varnames %in% colnames(Dat)]


# 1) how many NAs are there for each of these questions? Important to know,
# because we're about to do some major imputing!!

NApre <- sapply(varnames, function(vn, Dat){
			sum(is.na(Dat[[vn]]))
		}, Dat = Dat) # save this object to compare!

# for each of these varnames, lets interpolate to fill missings for waves in
# which the person was interviewed but not asked.
imputeSkippedQuestions <- function(vn,intv_dt){
	nas <- is.na(vn)
	if (all(nas)){
		return(vn)
	}
	if (sum(!nas)==1){
		vn <- approx(x = intv_dt,
				y = vn,
				xout = intv_dt,
				rule = 1:2,
				method = "constant")$y
		return(vn)
	}
	
	if (sum(!nas > 1)){
		
		vn <- approx(x = intv_dt,
				y = vn,
				xout = intv_dt,
				rule = 1:2,
				method = "linear")$y
		return(vn)
	}
	
}

# which varnames will need to be NA'd for entire waves afterwards?:
NAout <- sapply(varnames, function(vn, Dat){
			tapply(Dat[[vn]],Dat$wave,function(x){
						all(is.na(x))
					})
		}, Dat = Dat)
# this object will help re-impute NAs where necessary

Dat  <- data.table(Dat)
{
	############## yikes, it was faster to write all this than to figure out how to 
# do it elegantly in data.table....... sorry dear reader! it's more complicated
# than a simple column apply because intv_dt is needed as well
	Dat[,adl3_:=imputeSkippedQuestions(adl3_,intv_dt), by = list(id) ]
	Dat[,adl5_:=imputeSkippedQuestions(adl5_,intv_dt), by = list(id) ]
	Dat[,iadl3_:=imputeSkippedQuestions(iadl3_,intv_dt), by = list(id) ]
	Dat[,iadl5_:=imputeSkippedQuestions(iadl5_,intv_dt), by = list(id) ]
	Dat[,cesd:=imputeSkippedQuestions(cesd,intv_dt), by = list(id) ]
	Dat[,lim_work:=imputeSkippedQuestions(lim_work,intv_dt), by = list(id) ]
	Dat[,srh:=imputeSkippedQuestions(srh,intv_dt), by = list(id) ]
	Dat[,bmi:=imputeSkippedQuestions(bmi,intv_dt), by = list(id) ]
	Dat[,back:=imputeSkippedQuestions(back,intv_dt), by = list(id) ]
	Dat[,hosp:=imputeSkippedQuestions(hosp,intv_dt), by = list(id) ]
	
	Dat$hosp_stays <- as.numeric(Dat$hosp_stays) 
	Dat[,hosp_stays:=imputeSkippedQuestions(hosp_stays,intv_dt), by = list(id) ] # Error, investigate
	
	Dat[,hosp_nights:=imputeSkippedQuestions(hosp_nights,intv_dt), by = list(id) ]
	Dat[,nh:=imputeSkippedQuestions(nh,intv_dt), by = list(id) ]
	
	Dat$nh_stays <- as.numeric(Dat$nh_stays) 
	Dat[,nh_stays:=imputeSkippedQuestions(nh_stays,intv_dt), by = list(id) ] # Error, investigate
	Dat$nh_nights <- as.numeric(Dat$nh_nights) 
	Dat[,nh_nights:=imputeSkippedQuestions(nh_nights,intv_dt), by = list(id) ] # Error, investigate
	
	Dat[,nh_now:=imputeSkippedQuestions(nh_now,intv_dt), by = list(id) ]
	Dat[,doc:=imputeSkippedQuestions(doc,intv_dt), by = list(id) ]
	Dat[,hhc:=imputeSkippedQuestions(hhc,intv_dt), by = list(id) ]
	Dat[,meds:=imputeSkippedQuestions(meds,intv_dt), by = list(id) ]
	Dat[,surg:=imputeSkippedQuestions(surg,intv_dt), by = list(id) ]
	Dat[,dent:=imputeSkippedQuestions(dent,intv_dt), by = list(id) ]
	Dat[,shf:=imputeSkippedQuestions(shf,intv_dt), by = list(id) ]
	Dat[,adl_walk:=imputeSkippedQuestions(adl_walk,intv_dt), by = list(id) ]
	Dat[,adl_dress:=imputeSkippedQuestions(adl_dress,intv_dt), by = list(id) ]
	Dat[,adl_bath:=imputeSkippedQuestions(adl_bath,intv_dt), by = list(id) ]
	Dat[,adl_eat:=imputeSkippedQuestions(adl_eat,intv_dt), by = list(id) ]
	Dat[,adl_bed:=imputeSkippedQuestions(adl_bed,intv_dt), by = list(id) ]
	Dat[,adl_toilet:=imputeSkippedQuestions(adl_toilet,intv_dt), by = list(id) ]
	Dat[,iadl_map:=imputeSkippedQuestions(iadl_map,intv_dt), by = list(id) ]
	Dat[,iadl_tel:=imputeSkippedQuestions(iadl_tel,intv_dt), by = list(id) ]
	Dat[,iadl_money:=imputeSkippedQuestions(iadl_money,intv_dt), by = list(id) ]
	Dat[,iadl_meds:=imputeSkippedQuestions(iadl_meds,intv_dt), by = list(id) ]
	Dat[,iadl_shop:=imputeSkippedQuestions(iadl_shop,intv_dt), by = list(id) ]
	Dat[,iadl_meals:=imputeSkippedQuestions(iadl_meals,intv_dt), by = list(id) ]
	Dat[,mob:=imputeSkippedQuestions(mob,intv_dt), by = list(id) ]
	Dat[,lg_mus:=imputeSkippedQuestions(lg_mus,intv_dt), by = list(id) ]
	Dat[,gross_mot:=imputeSkippedQuestions(gross_mot,intv_dt), by = list(id) ]
	Dat[,fine_mot:=imputeSkippedQuestions(fine_mot,intv_dt), by = list(id) ]
	Dat[,bp:=imputeSkippedQuestions(bp,intv_dt), by = list(id) ]
	Dat[,diab:=imputeSkippedQuestions(diab,intv_dt), by = list(id) ]
	Dat[,cancer:=imputeSkippedQuestions(cancer,intv_dt), by = list(id) ]
	Dat[,lung:=imputeSkippedQuestions(lung,intv_dt), by = list(id) ]
	Dat[,heart:=imputeSkippedQuestions(heart,intv_dt), by = list(id) ]
	Dat[,stroke:=imputeSkippedQuestions(stroke,intv_dt), by = list(id) ]
	Dat[,psych:=imputeSkippedQuestions(psych,intv_dt), by = list(id) ]
	Dat[,arth:=imputeSkippedQuestions(arth,intv_dt), by = list(id) ]
	Dat[,cc:=imputeSkippedQuestions(cc,intv_dt), by = list(id) ]
	Dat[,alc_ev:=imputeSkippedQuestions(alc_ev,intv_dt), by = list(id) ]
	Dat[,alc_days:=imputeSkippedQuestions(alc_days,intv_dt), by = list(id) ]
	
	Dat$alc_drinks <- as.numeric(Dat$alc_drinks) # data.table needs consistent classes...
	Dat[,alc_drinks:=imputeSkippedQuestions(alc_drinks,intv_dt), by = list(id) ] # Error, investigate (see line above)
	
	Dat[,smoke_ev:=imputeSkippedQuestions(smoke_ev,intv_dt), by = list(id) ]
	Dat[,smoke_cur:=imputeSkippedQuestions(smoke_cur,intv_dt), by = list(id) ]
	Dat[,cesd_depr:=imputeSkippedQuestions(cesd_depr,intv_dt), by = list(id) ]
	Dat[,cesd_eff:=imputeSkippedQuestions(cesd_eff,intv_dt), by = list(id) ]
	Dat[,cesd_sleep:=imputeSkippedQuestions(cesd_sleep,intv_dt), by = list(id) ]
	Dat[,cesd_happy:=imputeSkippedQuestions(cesd_happy,intv_dt), by = list(id) ]
	Dat[,cesd_lone:=imputeSkippedQuestions(cesd_lone,intv_dt), by = list(id) ]
	Dat[,cesd_sad:=imputeSkippedQuestions(cesd_sad,intv_dt), by = list(id) ]
	Dat[,cesd_going:=imputeSkippedQuestions(cesd_going,intv_dt), by = list(id) ]
	Dat[,cesd_enjoy:=imputeSkippedQuestions(cesd_enjoy,intv_dt), by = list(id) ]
	Dat[,srm:=imputeSkippedQuestions(srm,intv_dt), by = list(id) ]
	Dat[,pastmem:=imputeSkippedQuestions(pastmem,intv_dt), by = list(id) ]
	Dat[,ss:=imputeSkippedQuestions(ss,intv_dt), by = list(id) ]
	Dat[,c20b:=imputeSkippedQuestions(c20b,intv_dt), by = list(id) ]
	Dat[,name_mo:=imputeSkippedQuestions(name_mo,intv_dt), by = list(id) ]
	Dat[,name_dmo:=imputeSkippedQuestions(name_dmo,intv_dt), by = list(id) ]
	Dat[,name_yr:=imputeSkippedQuestions(name_yr,intv_dt), by = list(id) ]
	Dat[,name_dwk:=imputeSkippedQuestions(name_dwk,intv_dt), by = list(id) ]
	Dat[,name_sci:=imputeSkippedQuestions(name_sci,intv_dt), by = list(id) ] 
	Dat[,name_cac:=imputeSkippedQuestions(name_cac,intv_dt), by = list(id) ]
	Dat[,name_pres:=imputeSkippedQuestions(name_pres,intv_dt), by = list(id) ]
	Dat[,name_pres:=imputeSkippedQuestions(name_pres,intv_dt), by = list(id) ]
	Dat[,vocab:=imputeSkippedQuestions(vocab,intv_dt), by = list(id) ]
	Dat[,tm:=imputeSkippedQuestions(tm,intv_dt), by = list(id) ]
	Dat[,dwr:=imputeSkippedQuestions(dwr,intv_dt), by = list(id) ]
	Dat[,twr:=imputeSkippedQuestions(twr,intv_dt), by = list(id) ]
	Dat[,iwr:=imputeSkippedQuestions(iwr,intv_dt), by = list(id) ]
	### again, sorry this was insanely bad coding.
}
# now re-insert NAs for waves that simply didn't include question X:


# this picks up almost everything...
#ImputeThese <- sapply(varnames, function(vn, Dat){
#    checkImpute <- any(tapply(Dat[[vn]],Dat$wave,function(x){
#    any(is.na(x)) & any(!is.na(x))
#  }))
#  },Dat=Dat)

imputeNAsInTheseVars <- colnames(NAout)[colSums(NAout) > 0]
NAout <- NAout[,imputeNAsInTheseVars]
waves <- 1:10
vn <-"name_cac"
for (vn in imputeNAsInTheseVars){
	wavesi <- waves[NAout[,vn]]
	Dat[[vn]][Dat$wave %in% wavesi] <- NA
}

# compare with NApre
NApost <- sapply(varnames, function(vn, Dat){
			sum(is.na(Dat[[vn]]))
		}, Dat = Dat) 

plot(NApost,NApre,asp=1)
abline(a=0,b=1)
hist(NApost / NApre) # OK so this wasn't pointless.
mean((NApre - NApost) / NApre, na.rm=TRUE)
save(Dat,file = "Data/Data_long_imputed.Rdata")


# now we compare before and after values for these variables.
# this is just for the sake of a variable appendix.


#DatIn         <- local(get(load("Data/thanos_long_v2_2.gz")))
## remove missed interviews
#DatIn         <- DatIn[!is.na(DatIn$intv_dt), ]
## reduce to deceased-only
#DatIn         <- DatIn[DatIn$dead == 1, ] # cut down size to reduce character searching in a moment
#DatIn         <- convertDates(DatIn)
#DatFinal <- local(get(load("Data/Data_long.Rdata")))
#SurfaceList <- local(get(load("Data/SurfaceList.Rdata")))
#DatIn$dwr                    <- DatIn$dr20w + DatIn$dr10w
#DatIn$iwr                     <- DatIn$ir20w + DatIn$ir10w
#
#DatFinal      <- DatFinal[DatFinal$age >= 65, ]
#DatFinal      <- DatFinal[!is.na(Dat$b_yr), ]
#DatFinal$Coh5 <- DatFinal$b_yr -  DatFinal$b_yr %% 5 
#Coh5keep <- c(1900, 1905, 1910, 1915, 1920, 1925, 1930)
#DatFinal      <- DatFinal[DatFinal$Coh5 %in% Coh5keep, ]
#
#KeepVec <- paste(DatFinal$id,DatFinal$intv_dt)
#rownames(DatIn) <- paste(DatIn$id,DatIn$intv_dt)
##all(KeepVec %in% rownames(DatIn)) TRUE
#DatIn <- DatIn[KeepVec, ]
#DatIn$Coh5 <- DatIn$b_yr -  DatIn$b_yr %% 5 
#nrow(DatIn);nrow(DatFinal) # Needed to do this for accurate tabulations...
# ------------------------paste(Dat$id,Dat$intv_dt)

# now, for each variable we find the first unique individual for 'each' unique reponse within each variable.
# then we find the resulting response in the DatOut df, producing a list of before-after responses. eek
# just use for short-long names...
#Meta <- read.csv( "Data/PercentThano.csv",stringsAsFactors=FALSE)
#Meta <- Meta[,c("Short","Long")]
#varnames<- varnames[varnames != "med_explog"]
#Variables <- list()
#for (vn in varnames){
#  
#  if (length(unique(DatFinal[[vn]])) < 13){
#    uniqueTab   <- suppressWarnings(table(DatIn[[vn]],exclude=c("NA",NA),useNA="no"))
#    unique1915  <- suppressWarnings(table(DatIn[[vn]][DatIn$Coh5==1915],exclude=c("NA",NA),useNA="no"))
#    responsesIn <- names(uniqueTab)
#    
#    firstIDs    <- sapply(responsesIn, function(res, DatIn){
#        DatIn$id[DatIn[[vn]] == res][1]
#      },DatIn=DatIn)
#    
#    # now iterate over IDs to get corresponding numeric equivalent:
#    responsesOut <- sapply(firstIDs, function(id, .vn, DatFinal){
#        DatOut[[vn]][DatOut$id == id][1]
#      }, .vn = vn, DatFinal=DatFinal)
#    dfout <- data.frame(Original = responsesIn, 
#      Recode = responsesOut, 
#      Count = as.integer(uniqueTab), 
#      Count1915 = as.integer(unique1915),
#      stringsAsFactors = FALSE)
#    dfout <- dfout[order(dfout$Recode), ]
#    Variables[[vn]] <- dfout
#  }
#}
#
#xtable
#
#varnames[!varnames %in% names(Variables)]
#library(Hmisc)
#DatIn[,varnames[varnames %in% colnames(DatIn)]]
# install.packages("Hmisc")
#latex(describe(DatIn[,varnames[varnames %in% colnames(DatIn)]]),file="")

# check for each variable if it's asked in every wave...
#
#Problems <- sapply(varnames,function(vn,Dat){
#			any(unlist(tapply(Dat[[vn]],Dat[["wave"]],function(x){
#						all(is.na(x))
#					})))
#		}, Dat = Dat)
#sum(Problems)
#Problems[Problems]
#vn <- "srh"
#Problems2 <- sapply(varnames,function(vn,Dat){
#			any(unlist(tapply(Dat[[vn]],Dat[["wave"]],function(x){
#										all(is.na(x))
#									}))[-1])
#		}, Dat = Dat)
#
#sum(Problems2)
#Problems2[Problems2]
#
#Problems3 <- sapply(varnames,function(vn,Dat){
#			sum(unlist(tapply(Dat[[vn]],Dat[["wave"]],function(x){
#										all(is.na(x))
#									})))
#		}, Dat = Dat)
#sort(Problems3)
#
#### med expenditure needs to be removed, even though it has a very clear thano pattern
### mprobev / mprob need to go : too inconsistent
#med_explog
#med_exp

# the others are either missing just the first or the first two interviews. 
# That's OK, will only affect left side.
# can use this info to trim, since it's effectively extrapolation.

###########################################################################
# provisional: make matrices of means to take a look.
###########################################################################


wmean <- function(x,w=rep(1,length(x))){
	if (length(x)==0){
		return(NA)
	}
	sum(x * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
}

#Dat <- local(get(load("Data/Data_long.Rdata")))
Dat <- local(get(load("Data/Data_long_imputed.Rdata")))
Dat <- data.table(Dat)

## for binning purposes, akin to 'completed age'
Dat$tafloor <- floor(Dat$ta)
Dat$cafloor <- floor(Dat$ca)
Dat$tafloor[Dat$tafloor<0]<-0
## I guess actual interview date could be some weeks prior to registered
## interview date? There are two negative thano ages at wave 4 otherwise, but
## still rather close. Likely died shortly after interview.
#Dat$tafloor[Dat$tafloor < 0] <- 0

# greater than 15 we lose precision, cut off:
Dat <- Dat[Dat$tafloor <= 15, ]
Dat <- Dat[Dat$cafloor >= 60, ]


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
		"iadl_calc", "mprob", "mprobev", "med_explog") # this is bigger than the final list...

varnames <- varnames[varnames %in% colnames(Dat)]

Dat      <- Dat[Dat$age >= 65, ]
Dat      <- Dat[!is.na(Dat$b_yr), ]
Dat$Coh5 <- Dat$b_yr -  Dat$b_yr %% 5 
Coh5keep <- c(1900, 1905, 1910, 1915, 1920, 1925, 1930)
Coh5     <- c(1905, 1910, 1915, 1920, 1925) # i.e. we use the preceding and subsequent cohorts for help fitting
Dat      <- Dat[Dat$Coh5 %in% Coh5keep, ]
# This is a sloppy old-school way this
Dat         <- data.frame(Dat)
SurfaceList <- list()
#Dat$mod_freq
varname <- "diab"
head(Dati)
library(data.table)
library(reshape2)
Dat <- as.data.frame(Dat)

for (varname in varnames){
	
	Dati <- Dat[, c("sex","tafloor","cafloor","Coh5","p_wt2",varname)]
	colnames(Dati)[ncol(Dati)] <- "V1"
	
	Mean <- 
			data.table(Dati)[,  list(V1 = wmean(V1,p_wt2)),
					by = list(sex,tafloor,cafloor,Coh5)]
	Mean <- data.frame(Mean)   
	
	MeanM <- acast(Mean[Mean$sex == "m", ],tafloor~cafloor~Coh5,value.var="V1" )
	MeanF <- acast(Mean[Mean$sex == "f", ],tafloor~cafloor~Coh5,value.var="V1" )
	
	#setnames(Dat,cols="V1",value=varname)
	
	SurfaceList[[varname]] <- list(Male = MeanM,Female = MeanF)
}
source("R/SurfMap.R")

# ------------------------------------------------------        
save(SurfaceList, file = "Data/SurfaceList.Rdata")

