# --------------------------------------
# TR
# This script recodes most variables, makes some new ones.
# This is the step after extracting from the RAND version M
# file in HRS_Rand_extract.R
# the next script to run after this one is CreateMatrices.R
# --------------------------------------

# Sets working directory for Tim's machines
if (system("hostname",intern=TRUE) %in% c("triffe-N80Vm","tim-ThinkPad-L440")){
  # if I'm on the laptop
  setwd("/home/tim/git/HLETTD")
} else {
	if (system("hostname",intern=TRUE) == "PC-403478"){
		# on MPIDR PC
		setwd("U://git//HLETTD")
	} else {
		# in that case I'm on Berkeley system, and other people in the dept can run this too
		setwd(paste0("/data/commons/",system("whoami",intern=TRUE),"/git/HLETTD"))
	}
}
getwd()

# get data into long format (takes some minutes)
if (!file.exists("Data/RAND_p_long.Rdata")){
	source("R/0_RAND_reshape.R")
}



library(lubridate)
library(data.table)

#---------------------------------------------------------
# start utility function preamble
#----------------------------------------------------------
# convert yes no coded questions into binary
convertYN <- function(x){
    xx                  <- rep(NA, length(x))
    xx[grepl("yes", x)] <- 1
    xx[grepl("no", x)]  <- 0
    invisible(xx)
}
#-----------------------------------------------------
# convert odd binary with first and second try into single binary
# TR: changed Aug 4, 2016. No more intermediate values: percent incorrect.
convertCI <-  function(x){
  xx                  <- rep(NA, length(x))
  
  xx[x == "1.correct"] <- 0
  xx[x == "2.correct, 1st try"] <- 0
  #xx[x == "1.correct, 2nd try"] <- .5
  xx[x == "1.correct, 2nd try"] <- 0
  xx[x == "0.incorrect"]  <- 1
  invisible(as.numeric(xx))
}
#-----------------------------------------------------
# convert CESD variables into binary
convertCESD <- function(x){
  xx <- rep(NA,length(x))
  xx[x == "0.no"]                    <- 0
  xx[x == "4. none or almost none"]  <- 0
  # TR: changed Aug 4, 2016. No more intermediate values: percent incorrect.
#  xx[x == "3. some of the time"]     <- .5
#  xx[x == "2. most of the time" ]    <- .75
  xx[x == "3. some of the time"]     <- 0
  xx[x == "2. most of the time" ]    <- 1
  xx[x ==  "1.yes"]                  <- 1
  xx[x ==  "1. all or almost all"]   <- 1
  xx
}

#-----------------------------------------------------
# convert dates to R internal format
convertDates <- function(Dat){
    # can't be done with apply because we can't have Date class matrices.
	# all date columns can be detected with the pattern _dt, it turns out.
    DateInd       <- grep(pattern = "_dt",colnames(Dat))
    for (i in DateInd){
        Dat[,i]    <- as.Date(Dat[, i], origin = "1960-1-1")
    }
    invisible(Dat)
}

#-----------------------------------------------------
# two functions to get exact years lived and left
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


# -------------------------------#
# Weight imputation function     #
# see code for annotation        #
# -------------------------------#
imputeWeights <- function(wt,intv_dt){
	# positive weights, also used for indexing
	ind <- wt > 0
	# if all weights 0, replace w NA
    if (all(wt == 0)){
        wt2 <- NA * wt
        return(wt2)
    }
	# if only one valid weight, all later
	# observations will keep that weight.
    if (sum(ind) == 1){
        wt2 <- approx(x = intv_dt[ind],
                y = wt[ind],
                xout = intv_dt,
                rule = 1:2,
                method = "constant",
                f = .5)$y
    }
	# if at least two valid observations, we
	# interpolate linearly for any missing,
	# but extrapolate (rightward) with constant
    if (sum(ind)>=2){
        wt2 <- approx(x = intv_dt[ind],
                y = wt[ind],
                xout = intv_dt,
                rule = 1:2,
                method = "linear")$y 
    }
    return(wt2)
}

# end utility function preamble
#----------------------------------------------------------

#----------------------------------------------------------
# load in long files from PHC
Dat         <- local(get(load("Data/RAND_p_long.Rdata")))
(start_cases <- nrow(Dat))
#--------------------------------------------------------------------------#
# we could expand to newer data, but would require a lot more work.        #
# would rather wait until more waves are out, linked, and integrated       #
# into the RAND data. Especially since mortality linking has been changing #
# and some changes will have retrospective impacts. Could cause bad        #
# headache to get into now.                                                #
#--------------------------------------------------------------------------#
#Dat         <- local(get(load("Data/thanos_long_v3_1.RData")))

# remove missed interviews
Dat         <- Dat[!is.na(Dat$intv_dt), ]

# TR: factors not present in later revisions.
# change all factors to character (to be later recoded in some instances)
## Dat[sapply(Dat, is.factor)] <- lapply(Dat[sapply(Dat, is.factor)], as.character)

# make sex column easier to use:
Dat$sex     <- ifelse(Dat$sex == "1.male","m","f")


# reduce to deceased-only
Dat$dead    <- ifelse(is.na(Dat$d_dt), 0, 1) 
Dat         <- Dat[Dat$dead == 1, ]
(dead_cases <- nrow(Dat))
# stats for paper
nrow(Dat) 
dead_cases / length(unique(Dat$id)) # 4.822015 / person on average
#hist(rle(sort(Dat$id))$lengths)

# convert dates to native R format
Dat         <- convertDates(Dat)

# --------------------------------------------------#
# merge weights (big assumption here:               #
# weights in institutions are valid and             #
# comparable with weights outside institutions.     #
# soooo annoying ppl in institutions don't have     #
# comparable weights.                               #
# --------------------------------------------------#
Dat$nh_wt[is.na(Dat$nh_wt)] <- 0
Dat$p_wt 	<- Dat$p_wt + Dat$nh_wt

# --------------------------------------------------#
# now we do weight interpolation/extrapolation      #
# --------------------------------------------------#
Dat 		<- data.table(Dat)
Dat 		<- Dat[, p_wt2 := imputeWeights(p_wt,intv_dt), by = list(id) ]
Dat 		<- Dat[!is.na(Dat$p_wt2),]

(cases_with_weights <- nrow(Dat))
cases_with_weights / length(unique(Dat$id)) # 4.683211

# --------------------------------------------------#
# calculate thanatological ages                     #
# --------------------------------------------------#
Dat$ta 		<- getThanoAge(Dat$intv_dt, Dat$d_dt)
Dat$ca 		<- getChronoAge(Dat$intv_dt, Dat$b_dt)
Dat$la_int 	<- floor(Dat$ta + Dat$ca)
# there is one individual with an NA b_dt, and NA age,
# but thano age is known
# --------------------------------------------------#
# locate yes/no, correct/incorrect columns          #
# --------------------------------------------------#
YNcols <- apply(Dat, 2, function(x){
        xx <- unique(x)
        length(xx) <= 4 & any(grepl("yes",xx))
        })
CIcols <- apply(Dat, 2, function(x){
          xx <- unique(x)
          length(xx) <= 5 & any(grepl("correct",xx))
        }) 

# which columns are these anyway?
colnames(Dat)[YNcols]
colnames(Dat)[CIcols] 

# convert to binary
Dat         <- data.frame(Dat)
Dat[YNcols] <- lapply(Dat[YNcols], convertYN)
Dat[CIcols] <- lapply(Dat[CIcols], convertCI)
#head(Dat)

# --------------------------------------------------#
# save off final names                              #
# --------------------------------------------------#
varnames_check <- local(get(load("Data/varnames.Rdata")))
all(varnames_check %in% colnames(Dat))
varnames_check[!varnames_check %in% colnames(Dat)]


names(varnames_check)      <- varnames_check
varnames_check[ "adl3_"]   <- "adl3"
varnames_check[ "adl5_"]   <- "adl5"
varnames_check[ "iadl3_"]  <- "iadl3"
varnames_check[ "iadl5_"]  <- "iadl5"

save(varnames_check, file = "Data/varnamesP.Rdata")

# mprob  = "pastmem" 
# iwr, twr, dwr all created below.
#grep("dwr",colnames(Dat))
#--------------------------------------------------------
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
# recode medical expenditure to mid-range values:


# recode self reported health to binary:
# excellent to good = 0, fair, poor = 1.
Dat$srh             <- as.character(Dat$srh)
srhrec              <- c(0,0,0,1,1,NA)
names(srhrec)       <- sort(unique(Dat$srh))

Dat$srm             <- as.character(Dat$srm)
names(srhrec)       <- sort(unique(Dat$srm))
Dat$srm             <- srhrec[Dat$srm] 

# now move to binary

# same, worse, better recode:  0 betterm 0 same 1 worse
Dat$pastmem         <- as.character(Dat$pastmem)
pastmem             <- c(0,0,1,NA)
names(pastmem)      <- sort(unique(Dat$pastmem))
Dat$pastmem         <- pastmem[Dat$pastmem] 

# TODO
# left off here June 18, 2017. Beyond this point still not examined!!!!
#

# do cesd questions (1 bad, 0 good)
cesdquestions       <- colnames(Dat)[grepl("cesd", colnames(Dat))]
cesdquestions       <- cesdquestions[cesdquestions != "cesd"]
Dat[cesdquestions]  <- lapply(Dat[cesdquestions],convertCESD)

# cesd_enjoy is flipped yet again, because 1 is 'yes I enjoyed life',
# and we want high = bad.
Dat$cesd_enjoy      <- 1 - Dat$cesd_enjoy
Dat$cesd_happy      <- 1 - Dat$cesd_happy

# ---------------------------------------------------------------
# create a single Total Word Recall variables, twr
#"tr20w"(waves(2-10),"tr40w" (waves1-2)
# i.e. 1 is the worst recall, and 0 is the best recall. This
# will operate the same as a binary var, but I don't want to make
# it binary because I wouldn't know where to set the breakpoint.
# plus I doubt it would affect the aggregate pattern conclusions anyway.
# same story for vocab, total memory, delayed word recall, 
# immediate word recall
# ---------------------------------------------------------------
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
Dat$dwr                     <- Dat$dr20w + Dat$dr10w
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
#colnames(Dat)
#names(mprob) <- sort(unique(Dat$mprob))
#Dat$mprob <- mprob[Dat$mprob]
# vocab

# ------------------------------------------------------
# TR: now change to binary coding, using ad hoc breaks
# scale to fit btwn 0 and 1
#rescale <- function(var,Dat,complement = FALSE){
#  Dat[[var]] <- Dat[[var]] / max(Dat[[var]], na.rm = TRUE)
#  if (compelment){
#    Dat[[var]] <- 1 - Dat[[var]]
#  }
#  Dat
#}
# ------------------------------------------------------

# TR: mob through cesd all change to binary, with ad hoc breaks
# -------------------
# mob 
# Dat     <- rescale("mob", Dat, FALSE)
hist(Dat$mob) # break at > 1
Dat$mob <- ifelse(is.na(Dat$mob), NA, ifelse(Dat$mob > 1, 1, 0) )

# -------------------
# lg_mus
#Dat     <- rescale("lg_mus", Dat, FALSE) 
hist(Dat$lg_mus) # break at > 1
Dat$lg_mus <- ifelse(is.na(Dat$lg_mus), NA, ifelse(Dat$lg_mus > 1, 1, 0) )

# -------------------
# gross_mot
#Dat     <- rescale("gross_mot", Dat, FALSE)
hist(Dat$gross_mot) # > 1
Dat$gross_mot <- ifelse(is.na(Dat$gross_mot), NA, ifelse(Dat$gross_mot > 1, 1, 0) )

# -------------------
# gross_mot
#Dat     <- rescale("fine_mot", Dat, FALSE)
hist(Dat$fine_mot) # > 0
Dat$fine_mot <- ifelse(is.na(Dat$fine_mot), NA, ifelse(Dat$fine_mot > 0, 1, 0) )

# -------------------
# ss
#Dat     <- rescale("ss", Dat, TRUE) # complement because more was better in original
hist(Dat$ss) # < 4
Dat$ss <- ifelse(is.na(Dat$ss), NA, ifelse(Dat$ss < 4, 1, 0) )

# -------------------
# cc nr chronic cond
#Dat     <- rescale("cc", Dat, FALSE)
hist(Dat$cc) # > 2?
Dat$cc <- ifelse(is.na(Dat$cc), NA, ifelse(Dat$cc > 2, 1, 0) )

# -------------------
# alc_days
#Dat     <- rescale("alc_days", Dat, FALSE)
hist(Dat$alc_days) # > 1
Dat$alc_days <- ifelse(is.na(Dat$alc_days), NA, ifelse(Dat$alc_days > 1, 1, 0) )

# -------------------
# adl3_
#Dat     <- rescale("adl3_", Dat, FALSE)
hist(Dat$adl3_) # > 0
Dat$adl3_ <- ifelse(is.na(Dat$adl3_), NA, ifelse(Dat$adl3_ > 0, 1, 0) )

# -------------------
# adl5_
#Dat     <- rescale("adl5_", Dat, FALSE)
hist(Dat$adl5_) # > 0
Dat$adl5_ <- ifelse(is.na(Dat$adl5_), NA, ifelse(Dat$adl5_ > 0, 1, 0) )

# -------------------
# iadl3_
# Dat     <- rescale("iadl3_", Dat, FALSE)
hist(Dat$iadl3_) # > 0
Dat$iadl3_ <- ifelse(is.na(Dat$iadl3_), NA, ifelse(Dat$iadl3_ > 0, 1, 0) )

# -------------------
# iadl5_
#Dat     <- rescale("iadl5_", Dat, FALSE)
hist(Dat$iadl5_) # > 0
Dat$iadl5_ <- ifelse(is.na(Dat$iadl5_), NA, ifelse(Dat$iadl5_ > 0, 1, 0) )

# -------------------
# cesd
#Dat     <- rescale("cesd", Dat, FALSE)
hist(Dat$cesd) # > 2
Dat$cesd <- ifelse(is.na(Dat$cesd), NA, ifelse(Dat$cesd > 2, 1, 0) )

# -------------------
# TODO: does high = bad for these variables?
# it doesn't matter if all we want to do is measure
# which direction things correlate in, but it's nice
# to know for visual interpretation of surfaces. 
# "bmi"        "back"       "dent"       "alc_ev"   
# "pastmem"    "dwr"        "twr"        "iwr" 
# -------------------


# -------------------
# check cases by wave ( tapering in recent waves because selected down to deaths..)
#checkwaves <- function(var,Dat){
#  table(Dat[[var]],Dat[["wave"]])
#}
#checkwaves("adl3_",Dat)
#checkwaves("adl5_",Dat)
#checkwaves("iadl3_",Dat)
#checkwaves("iadl5_",Dat)
#checkwaves("cesd",Dat)
# -------------------

# -------------------------------------------------------
# for binning purposes, akin to 'completed age'
Dat$tafloor <- floor(Dat$ta)
Dat$cafloor <- floor(Dat$ca)

# I guess actual interview date could be some weeks prior to registered
# interview date? There are two negative thano ages at wave 4 otherwise, but
# still rather close. Likely died shortly after interview.
Dat$tafloor[Dat$tafloor < 0] <- 0

# We use higher bin widths fur purposes of visualizing raw data,
# Just diagnostics. larger widths help cancel out noise. This
# Such binning can be done as an alternative to the loess smoothing,
# where we take weighted means in cells. It'd probably make sense
# to keep the final year of life in a single year width, but the 
# general pattern ought to remain visible.
Dat$cafloor2 <- Dat$cafloor - Dat$cafloor %% 2
Dat$tafloor2 <- Dat$tafloor - Dat$tafloor %% 2

Dat$cafloor3 <- Dat$cafloor - Dat$cafloor %% 3
Dat$tafloor3 <- Dat$tafloor - Dat$tafloor %% 3
#----------------------------------------------
# save out, so this doesn't need to be re-run every time
save(Dat,file = "Data/Data_long.Rdata")
# next step would be CreateMatrices.R, usually

#----------------------------------------------


#Dat <- local(get(load("Data/Data_long.Rdata")))

#apply(Dat[,colnames(Dat)%in%varnames],2,range,na.rm=TRUE)
#----------------------------------------------
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

# ---------------------------------------------
# optional imputing of questions. Asked stats geeks, and
# multiple imputation is not advised. Within individual 
# trajectories are probably innocuous, as well, just
# increases number of points. Let's generate downstream 
# results using the non-imputed version. Ergo, the following
# code chunk is not used for the published analysis
# ---------------------------------------------

skip <- TRUE
if (!skip){
	
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
}

# end


# ------------------------------------------------------
# DEPRECATED. MAY NO LONGER PERTAIN TO ABOVE CODE
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
