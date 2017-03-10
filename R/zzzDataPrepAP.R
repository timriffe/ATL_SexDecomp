

# this jsut needs to be repeated because the other dataset omits alive people...

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

Dat           <- local(get(load("Data/thanos_long_v3_1.RData")))
# remove missed interviews
Dat           <- Dat[!is.na(Dat$intv_dt), ]
# make sex column easier to use:
Dat$sex       <- ifelse(Dat$sex == "1.male","m","f")
Dat           <- convertDates(Dat)
Dat           <- data.table(Dat)
Dat$nh_wt[is.na(Dat$nh_wt)] <- 0
Dat$p_wt      <- Dat$p_wt + Dat$nh_wt
Dat           <- Dat[, p_wt2 := imputeWeights(p_wt, intv_dt), by = list(id) ]
Dat           <- Dat[!is.na(Dat$p_wt2),]
Dat$age       <- getChronoAge(Dat$intv_dt, Dat$b_dt)





rescale <- function(var,Dat,compelment = FALSE){
	Dat[[var]] <- Dat[[var]] / max(Dat[[var]], na.rm = TRUE)
	if (compelment){
		Dat[[var]] <- 1 - Dat[[var]]
	}
	Dat
}
Dat     <- rescale("adl3_", Dat, FALSE)
Dat     <- rescale("iadl3_", Dat, FALSE)
Dat     <- rescale("adl5_", Dat, FALSE)
Dat     <- rescale("iadl5_", Dat, FALSE)

setnames(Dat,"adl5_","adl5")
setnames(Dat,"iadl5_","iadl5")
setnames(Dat,"adl3_","adl3")
setnames(Dat,"iadl3_","iadl3")
## let's just smooth these in AP.

Dat$DateDec <- lubridate::decimal_date(Dat$intv_dt)
colnames(Dat)
Dat <- Dat[Dat$age > 65, ]

Males       <- as.data.frame(Dat)[Dat$sex == "m", ]
Females     <-  as.data.frame(Dat)[Dat$sex == "f", ]


fitAPloess <- function(DatSex, varname = "adl5"){
	lubridate::decimal_date(Males$intv_dt[1])
	mod   <- loess(paste0(varname,'~age * DateDec') ,
			data = DatSex, 
			weights = p_wt2, # this, plus point density both act as weights
			span = .5,       # a variable passed in, or smoothness
			normalize = FALSE,
			family = "gaussian",
			control = loess.control(trace.hat="approximate")
	)
	# predict 1992 to 2013
	newdata        <- expand.grid(age = 65:105 + .5, DateDec = 1992:2013 + .5)
	Surf           <- predict(mod, newdata)
	dimnames(Surf) <- list(65:105,1992:2013)
	Surf
}

iadlm5 <- fitAPloess(Males,"iadl5")
iadlf5 <- fitAPloess(Females,"iadl5")
adlm5 <- fitAPloess(Males,"adl5")
adlf5 <- fitAPloess(Females,"adl5")

adlm3 <- fitAPloess(Males,"adl3")
adlf3 <- fitAPloess(Females,"adl3")

save(adlm3, file = "Data/APadlm3.Rdata")
save(adlf3, file = "Data/APadlf3.Rdata")
save(adlm5, file = "Data/APadlm5.Rdata")
save(adlf5, file = "Data/APadlf5.Rdata")







