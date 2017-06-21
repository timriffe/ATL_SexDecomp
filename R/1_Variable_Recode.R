# --------------------------------------
# TR
# This script recodes most variables, makes some new ones.
# This is the step after extracting from the RAND version M
# file in HRS_Rand_extract.R
# the next script to run after this one is CreateMatrices.R
# --------------------------------------
# source("R/1_Variable_Recode.R")
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
	x <- as.character(x)
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
  
  x <- as.character(x)
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
  if (class(x) == "numeric"){
		return(x)
  }
	
  x <- as.character(x)
  xx <- rep(NA,length(x))
  xx[x == "0.no"]                    <- 0
  xx[x == "4. none or almost none"]  <- 0
  xx[x == "4.none or almost none"]  <- 0
  # TR: changed Aug 4, 2016. No more intermediate values: percent incorrect.
#  xx[x == "3. some of the time"]     <- .5
#  xx[x == "2. most of the time" ]    <- .75
  xx[x == "3. some of the time"]     <- 0
  xx[x == "3.some of the time"]     <- 0
  xx[x == "2. most of the time" ]    <- 1
  xx[x == "2.most of the time" ]    <- 1
  xx[x ==  "1.yes"]                  <- 1
  xx[x ==  "1. all or almost all"]   <- 1
  xx[x ==  "1.all or almost all"]   <- 1
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

#Dat[Dat$id==53901020,] # investigate this silly case
#(start_cases <- nrow(Dat)) #
# we could expand to newer data, but would require a lot more work.        #
# would rather wait until more waves are out, linked, and integrated       #
# into the RAND data. Especially since mortality linking has been changing #
# and some changes will have retrospective impacts. Could cause bad        #
# headache to get into now.                                                #
#--------------------------------------------------------------------------#
#Dat         <- local(get(load("Data/thanos_long_v3_1.RData")))

# remove missed interviews
Dat         <- Dat[!is.na(Dat$intv_dt), ]
#Dat[Dat$intv == "1.resp,alive",]
# TR: factors not present in later revisions.
# change all factors to character (to be later recoded in some instances)
## Dat[sapply(Dat, is.factor)] <- lapply(Dat[sapply(Dat, is.factor)], as.character)

# make sex column easier to use:
Dat$sex     <- ifelse(Dat$sex == "1.male","m","f")


# reduce to deceased-only
Dat$dead    <- ifelse(is.na(Dat$d_dt), 0, 1) 
Dat         <- Dat[Dat$dead == 1, ]
#(dead_cases <- nrow(Dat))
# stats for paper
#nrow(Dat) 
#dead_cases / length(unique(Dat$id)) # 4.822015 / person on average
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

#(cases_with_weights <- nrow(Dat))
#cases_with_weights / length(unique(Dat$id)) # 4.683211

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
varnames <- local(get(load("Data/varnamesP.Rdata")))

YNcols <- apply(Dat, 2, function(x){
        xx <- unique(x)
        length(xx) <= 4 & any(grepl("yes",xx))
        })
CIcols <- apply(Dat, 2, function(x){
          xx <- unique(x)
          length(xx) <= 5 & any(grepl("correct",xx))
        }) 

# which columns are these anyway?
#colnames(Dat)[YNcols]
#colnames(Dat)[CIcols] 

# convert to binary
Dat         <- data.frame(Dat)
Dat[YNcols] <- lapply(Dat[YNcols], convertYN)
Dat[CIcols] <- lapply(Dat[CIcols], convertCI)
#head(Dat)




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
Dat$cesd_happy   <- NULL # too few non-NA cases
Dat$iadl_tel  <- NULL    # too few non-NA cases
### med expenditure needs to be removed, even though it has a very clear thano pattern
## mprobev / mprob need to go : too inconsistent
Dat$mprob 		<- NULL
Dat$mprobev 	<- NULL
Dat$med_explog 	<- NULL
Dat$med_exp 	<- NULL
# recode medical expenditure to mid-range values:
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
varnames_check <- varnames_check[varnames_check %in% colnames(Dat)]
save(varnames_check, file = "Data/varnamesP.Rdata")

#pdf("Figures/histograms.pdf")
#lapply(varnames_check,function(vname,Dat){
#			x <- Dat[[vname]]
#			if (length(unique(x))>3){
#				hist(x,breaks=unique(as.integer(pretty(x,n=25))),main=vname)	
#			}
#			},Dat=Dat)
#dev.off()
# found this out: should be integer...
Dat$hosp_nights <- floor(Dat$hosp_nights)
binaries <- unlist(lapply(Dat[varnames_check],function(x){
			all(x[!is.na(x)] %in% c(0,1))
		}))
int <- unlist(lapply(Dat[varnames_check],function(x){
					x <- x[!is.na(x)]
					all(x == as.integer(x))
				}))
#range(Dat$hosp_nights-floor(Dat$hosp_nights),na.rm=TRUE)
#write.csv(data.frame(Morbidity=varnames_check,
#				binary=binaries,
#				count = int & !binaries, 
#				other = !int & !binaries),
#		file = "Data/variabletypesP_untransformed.csv",row.names=FALSE)
#save(Dat,file="Data/Praw.Rdata")

# recode self reported health to binary:
# excellent to good = 0, fair, poor = 1.
Dat$srh             <- as.character(Dat$srh)
srhrec              <- c(0,0,0,1,1,NA)
names(srhrec)       <- sort(unique(Dat$srh))
Dat$srh 			<- srhrec[Dat$srh]

Dat$srm             <- as.character(Dat$srm)
names(srhrec)       <- sort(unique(Dat$srm))
Dat$srm             <- srhrec[Dat$srm] 

# now move to binary

# same, worse, better recode:  0 betterm 0 same 1 worse
Dat$pastmem         <- as.character(Dat$pastmem)
pastmem             <- c(0,0,1,NA)
names(pastmem)      <- sort(unique(Dat$pastmem))
Dat$pastmem         <- pastmem[Dat$pastmem] 



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
# immediate word recall. A quasibinom will work here
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
# alc_days (any)
#Dat     <- rescale("alc_days", Dat, FALSE)
hist(Dat$alc_days) # > 1
Dat$alc_days <- ifelse(is.na(Dat$alc_days), NA, ifelse(Dat$alc_days > 1, 1, 0) )

# -------------------
# adl3 (any)
#Dat     <- rescale("adl3_", Dat, FALSE)
hist(Dat$adl3) # > 0
Dat$adl3 <- ifelse(is.na(Dat$adl3), NA, ifelse(Dat$adl3 > 0, 1, 0) )

# -------------------
# adl5 (any)
#Dat     <- rescale("adl5_", Dat, FALSE)
hist(Dat$adl5) # > 0
Dat$adl5 <- ifelse(is.na(Dat$adl5), NA, ifelse(Dat$adl5 > 0, 1, 0) )

# -------------------
# iadl3 (any)
# Dat     <- rescale("iadl3_", Dat, FALSE)
hist(Dat$iadl3) # > 0
Dat$iadl3 <- ifelse(is.na(Dat$iadl3), NA, ifelse(Dat$iadl3 > 0, 1, 0) )
# -------------------
# iadl5 (any)
#Dat     <- rescale("iadl5_", Dat, FALSE)
hist(Dat$iadl5) # > 0
Dat$iadl5 <- ifelse(is.na(Dat$iadl5), NA, ifelse(Dat$iadl5 > 0, 1, 0) )

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
hist(Dat$ta)
# I guess actual interview date could be some weeks prior to registered
# interview date? There are two negative thano ages at wave 4 otherwise, but
# still rather close. Likely died shortly after interview.
(Dat[Dat$tafloor < 0, ])
# there is one individual with an erroneous death date (or id!), throwing out.
Dat                          <- Dat[Dat$ta > -1, ]
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
save(Dat,file = "Data/Data_longP.Rdata")
# next step would be CreateMatrices.R, usually

#lapply(Dat[varnames_check],unique)
# how many NAs per variable?
#varnames_check[!varnames_check%in%colnames(Dat)]
#NAs <- lapply(Dat[varnames_check],function(x){
#			sum(is.na(x))
#		})
#sort(unlist(NAs) / nrow(Dat))
# end
