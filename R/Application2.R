

# for TR computers...otherwise you need to set wd yourself
if (system("hostname",intern=TRUE) %in% c("triffe-N80Vm","tim-ThinkPad-L440")){
	# if I'm on the laptop
	setwd("/home/tim/git/ATL_SexDecomp/ATL_SexDecomp")
} 
getwd()
# install.packages("lubridate")
library(lubridate)
library(data.table)


# cleaning/processing functions
# need to make some codes usable...


# this code imported from ThanoEmpirical, but it is in need of some serious overhaul
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
	out[Ind] <- decimal_date(Date[Ind]) - decimal_date(BirthDate[Ind])
	out
}
getThanoAge <- function(Date, DeathDate){
	out <- rep(NA, length(Date))
	Ind <- !is.na(Date) & !is.na(DeathDate)
	out[Ind] <- decimal_date(DeathDate[Ind]) - decimal_date(Date[Ind])
	out
}

# TR: This is sketchy. To be replaced one day
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
# dead only = full coordinates:
Dat           <- Dat[Dat$dead == 1, ]

# remove missed interviews
Dat           <- Dat[!is.na(Dat$intv_dt), ]
# make sex column easier to use:
Dat$sex       <- ifelse(Dat$sex == "1.male","m","f")

# make date R-friendly
Dat           <- convertDates(Dat)

# data.table for easier processing of some stuff
Dat           <- data.table(Dat)

# for some years there are separate weights for those in care homes vs
# society at large. In this case, they add, because the default weights
# are set to zero
Dat$nh_wt[is.na(Dat$nh_wt)] <- 0
Dat$p_wt      <- Dat$p_wt + Dat$nh_wt

# but for other years the care home pop gets zero weights
# but we don't want to lose them, ergo: sketchy imputation 
# For some reason HRS gives them zero weight, 
# but they are clearly in-universe...
Dat           <- Dat[, p_wt2 := imputeWeights(p_wt, intv_dt), by = list(id) ]
Dat           <- Dat[!is.na(Dat$p_wt2),]


# decimal values for chrono-thano age
Dat$age       <- getChronoAge(Dat$intv_dt, Dat$b_dt)
Dat$ttd       <- getThanoAge(Dat$intv_dt, Dat$d_dt)

# single age classes (not used)
Dat$ca        <- floor(Dat$age)
Dat$ta        <- floor(Dat$ttd)

# make 5 year cohorts
Dat$coh5      <- Dat$b_yr - Dat$b_yr %% 5

# Coh5keep inlcudes neighboring outside cohorts for help fitting
Coh5keep <- c(1900, 1905, 1910, 1915, 1920, 1925, 1930)

# these are the cohorts we predict for
Coh5     <- c(1905, 1910, 1915, 1920, 1925) 
Dat      <- Dat[Dat$coh5 %in% Coh5keep, ]

unique(Dat$adl3_)
unique(Dat$adl5_)
Dat$adl1            <- as.integer(Dat$adl5_ >= 1)
Dat$adl2            <- as.integer(Dat$adl5_ >= 2)
Dat$adl3            <- as.integer(Dat$adl5_ >= 3)

colnames(Dat)
Dat$iadl1            <- as.integer(Dat$iadl5_ >= 1)
Dat$iadl2            <- as.integer(Dat$iadl5_ >= 2)
Dat$iadl3            <- as.integer(Dat$iadl5_ >= 3)


Dat$srhpoor         <- ifelse(Dat$srh == "5. poor",1,ifelse(Dat$srh == "NA",NA,0))


#Dat[, srhpoor := imputeSkippedQuestions(srhpoor,intv_dt), by = list(id) ]


#rescale <- function(var,Dat,compelment = FALSE){
#	Dat[[var]] <- Dat[[var]] / max(Dat[[var]], na.rm = TRUE)
#	if (compelment){
#		Dat[[var]] <- 1 - Dat[[var]]
#	}
#	Dat
#}
#Dat     <- rescale("adl3_", Dat, FALSE)
#Dat     <- rescale("iadl3_", Dat, FALSE)
#Dat     <- rescale("adl5_", Dat, FALSE)
#Dat     <- rescale("iadl5_", Dat, FALSE)
#
#setnames(Dat,"adl5_","adl5")
#setnames(Dat,"iadl5_","iadl5")
#setnames(Dat,"adl3_","adl3")
#setnames(Dat,"iadl3_","iadl3")
## let's just smooth these in AP.

Dat$DateDec <- lubridate::decimal_date(Dat$intv_dt)

# we want to be clearly in Gompertzlandia
Dat <- Dat[Dat$age > 65, ]
# now we do a simplistic smooth of the cloud of points in the space.
FitLoess <- function(varname, 
		Dat, 
		sex,
		t.age = 0:12,    # some 5-year cohorts simply don't have 15 years, cut it lower
		c.age = 70:100,  # standard matrix size, though we may NA certain unobserved cells
		span = .5, # will vary
		.Coh5){
	# conservative here to cut tails
	maxL  <- 100
	minL  <- 70
	Dats <- Dat[Dat$sex == sex, ]

	# multiplicative give the most freedom.
	mod   <- loess(paste0(varname,'~coh5 * ta * ca') ,
			data = Dats, 
			weights = p_wt2, # this, plus point density both act as weights
			span = span,     # a variable passed in, or smoothness
			# is similiar conceptually to a 1:1:1 aspect ratio. Everything is in years...
			normalize = FALSE,
			control = loess.control(trace.hat="approximate")
	)
	
	newdata        <- expand.grid(ta = t.age, ca = c.age, coh5 = .Coh5)
	# easier to keep dimensions straight if we predict over rectangular grid, 
	# then throw out values outside range
	#newdata        <- newdata + .5
	Surf           <- predict(mod, newdata)
	
	dimnames(Surf) <- list(floor(t.age),floor(c.age), .Coh5)
	
	# need to trim on the left side where applicable, since some questions didn't enter until
	# wave 2 or 3. There are many such cases, so check and make sure we dont' extrapolate.
	# sex <- "m"
	# varname <- "cesd"
	MissingWaves <- tapply(Dat[Dat$sex==sex,varname],Dat[Dat$sex==sex,"wave"], function(x){
				all(is.na(x))
			})
	RightYear <- 2011; LeftYear <- 1992
#	# this reduces extrapolation outside of data points 
	for (i in 1:dim(Surf)[3]){
		#maxL  <- 2011 - Coh5[i] - 1
		#maxt  <- tamax[as.character(Coh5[i])]
		#keept <- as.integer(rownames(Surf)) <= maxt
		A     <- Surf[,,i]
		MaxL <- RightYear - .Coh5[i] - 1
		A[ col(A) - 1 + 70 + row(A) - 1 > MaxL] <- NA
# possibly need to trim lower left corner too: dimnames(A)
		MinL <- LeftYear - (.Coh5[i] + 5)
		A[col(A) + 70 - 1 < MinL] <- NA
		#A[!keept, ] <- NA 
		Surf[,,i] <- A
	}
	list(Surf = Surf, span = span, sex = sex, varname = varname, Cohorts = Coh5)
}

Dat <- as.data.frame(Dat)

Female <- FitLoess(varname = "srhpoor", 
		Dat = Dat, 
		sex = "f",
		t.age = 0:12,    # some 5-year cohorts simply don't have 15 years, cut it lower
		c.age = 70:100,  # standard matrix size, though we may NA certain unobserved cells
		span =  .5, 
		.Coh5 = Coh5)

Male <- FitLoess(varname = "srhpoor", 
		Dat = Dat, 
		sex = "m",
		t.age = 0:12,    # some 5-year cohorts simply don't have 15 years, cut it lower
		c.age = 70:100,  # standard matrix size, though we may NA certain unobserved cells
		span =  .5, 
		.Coh5 = Coh5)


plot(rowMeans(Female$Surf, na.rm=TRUE), ylim=c(0,.4), xlim=c(0,15), type='l',col="red")
lines(rowMeans(Male$Surf, na.rm=TRUE),col="blue")

ADL1f <- FitLoess(varname = "adl1", 
		Dat = Dat, 
		sex = "f",
		t.age = 0:12,    # some 5-year cohorts simply don't have 15 years, cut it lower
		c.age = 70:100,  # standard matrix size, though we may NA certain unobserved cells
		span =  .5, 
		.Coh5 = Coh5)

ADL1m <- FitLoess(varname = "adl1", 
		Dat = Dat, 
		sex = "m",
		t.age = 0:12,    # some 5-year cohorts simply don't have 15 years, cut it lower
		c.age = 70:100,  # standard matrix size, though we may NA certain unobserved cells
		span =  .5, 
		.Coh5 = Coh5)
ADL2f <- FitLoess(varname = "adl2", 
		Dat = Dat, 
		sex = "f",
		t.age = 0:12,    # some 5-year cohorts simply don't have 15 years, cut it lower
		c.age = 70:100,  # standard matrix size, though we may NA certain unobserved cells
		span =  .5, 
		.Coh5 = Coh5)

ADL2m <- FitLoess(varname = "adl2", 
		Dat = Dat, 
		sex = "m",
		t.age = 0:12,    # some 5-year cohorts simply don't have 15 years, cut it lower
		c.age = 70:100,  # standard matrix size, though we may NA certain unobserved cells
		span =  .5, 
		.Coh5 = Coh5)
ADL3f <- FitLoess(varname = "adl3", 
		Dat = Dat, 
		sex = "f",
		t.age = 0:12,    # some 5-year cohorts simply don't have 15 years, cut it lower
		c.age = 70:100,  # standard matrix size, though we may NA certain unobserved cells
		span =  .5, 
		.Coh5 = Coh5)

ADL3m <- FitLoess(varname = "adl3", 
		Dat = Dat, 
		sex = "m",
		t.age = 0:12,    # some 5-year cohorts simply don't have 15 years, cut it lower
		c.age = 70:100,  # standard matrix size, though we may NA certain unobserved cells
		span =  .5, 
		.Coh5 = Coh5)

# IADL
IADL1f <- FitLoess(varname = "iadl1", 
		Dat = Dat, 
		sex = "f",
		t.age = 0:12,    # some 5-year cohorts simply don't have 15 years, cut it lower
		c.age = 70:100,  # standard matrix size, though we may NA certain unobserved cells
		span =  .5, 
		.Coh5 = Coh5)

IADL1m <- FitLoess(varname = "iadl1", 
		Dat = Dat, 
		sex = "m",
		t.age = 0:12,    # some 5-year cohorts simply don't have 15 years, cut it lower
		c.age = 70:100,  # standard matrix size, though we may NA certain unobserved cells
		span =  .5, 
		.Coh5 = Coh5)
IADL2f <- FitLoess(varname = "iadl2", 
		Dat = Dat, 
		sex = "f",
		t.age = 0:12,    # some 5-year cohorts simply don't have 15 years, cut it lower
		c.age = 70:100,  # standard matrix size, though we may NA certain unobserved cells
		span =  .5, 
		.Coh5 = Coh5)

IADL2m <- FitLoess(varname = "iadl2", 
		Dat = Dat, 
		sex = "m",
		t.age = 0:12,    # some 5-year cohorts simply don't have 15 years, cut it lower
		c.age = 70:100,  # standard matrix size, though we may NA certain unobserved cells
		span =  .5, 
		.Coh5 = Coh5)
IADL3f <- FitLoess(varname = "iadl3", 
		Dat = Dat, 
		sex = "f",
		t.age = 0:12,    # some 5-year cohorts simply don't have 15 years, cut it lower
		c.age = 70:100,  # standard matrix size, though we may NA certain unobserved cells
		span =  .5, 
		.Coh5 = Coh5)

IADL3m <- FitLoess(varname = "iadl3", 
		Dat = Dat, 
		sex = "m",
		t.age = 0:12,    # some 5-year cohorts simply don't have 15 years, cut it lower
		c.age = 70:100,  # standard matrix size, though we may NA certain unobserved cells
		span =  .5, 
		.Coh5 = Coh5)
graphics.off()
plot(rowMeans(ADL1f$Surf, na.rm=TRUE), ylim=c(0,.8), xlim=c(0,15),col="red",type='l')
lines(rowMeans(ADL1m$Surf, na.rm=TRUE),col="blue")
lines(rowMeans(ADL2m$Surf, na.rm=TRUE),col="blue",lty=2)
lines(rowMeans(ADL2f$Surf, na.rm=TRUE),col="red",lty=2)
lines(rowMeans(ADL3m$Surf, na.rm=TRUE),col="blue",lty=4)
lines(rowMeans(ADL3f$Surf, na.rm=TRUE),col="red",lty=4)

plot(rowMeans(IADL1f$Surf, na.rm=TRUE), ylim=c(0,.8), xlim=c(0,15),col="red",type='l')
lines(rowMeans(IADL1m$Surf, na.rm=TRUE),col="blue")
lines(rowMeans(IADL2m$Surf, na.rm=TRUE),col="blue",lty=2)
lines(rowMeans(IADL2f$Surf, na.rm=TRUE),col="red",lty=2)
lines(rowMeans(IADL3m$Surf, na.rm=TRUE),col="blue",lty=4)
lines(rowMeans(IADL3f$Surf, na.rm=TRUE),col="red",lty=4)

DatForAlyson <- list()
DatForAlyson[["SRHpoor"]] <- list(Males = Male, Females = Female)
DatForAlyson[["IADL1"]] <- list(Males = IADL1m, Females = IADL1f)
DatForAlyson[["IADL2"]] <- list(Males = IADL2m, Females = IADL2f)
DatForAlyson[["IADL3"]] <- list(Males = IADL3m, Females = IADL3f)
DatForAlyson[["ADL1"]] <- list(Males = ADL1m, Females = ADL1f)
DatForAlyson[["ADL2"]] <- list(Males = ADL2m, Females = ADL2f)
DatForAlyson[["ADL3"]] <- list(Males = ADL3m, Females = ADL3f)

save(DatForAlyson, file = "/home/tim/git/ATL_SexDecomp/ATL_SexDecomp/Data/SmArrays.Rdata")

Dat <- local(get(load("/home/tim/git/ATL_SexDecomp/ATL_SexDecomp/Data/SmArrays.Rdata")))

library(RColorBrewer)
Surf<- Dat$IADL1$Males$Surf[,,"1915"]
Ramp <- colorRampPalette(brewer.pal(9,"Reds"),space="Lab")
getwd()

ttd <- 12
age <- 30

pdf("/home/tim/git/ATL_SexDecomp/ATL_SexDecomp/PAApresentation/Figures/IADLttdlines.pdf",width=ttd/3+1.2,height=5)
par(xpd=TRUE,xaxs='i',yaxs='i',mai=c(.8,.8,.3,.3))
matplot(0:12,Surf, type= 'l', lty = 1, col = Ramp(ncol(Surf)+5)[-c(1:5)], lwd=2, axes=FALSE,xlab="",ylab="")
segments(0,0,12,0)
segments(c(0,5,10),0,c(0,5,10),-.01)
text(c(0,5,10),-.01,c(0,5,10),pos=1)
segments(0,0,0,.8)
segments(0,seq(0,.8,by=.2),-.15,seq(0,.8,by=.2))
text(0,seq(0,.8,by=.2),seq(0,.8,by=.2),pos=2)
text(6,-.08,"time-to-death",cex=1.5)
text(-1.2,.45,"Prevalence",srt=90,cex=1.5,pos=2)
dev.off()

#----------------------------------
pdf("/home/tim/git/ATL_SexDecomp/ATL_SexDecomp/PAApresentation/Figures/IADLagelines.pdf",width=age/3+1.2,height=5)
ages <- as.integer(colnames(Surf))
par(xpd=TRUE,xaxs='i',yaxs='i',mai=c(.8,.8,.3,.3))
matplot(as.integer(colnames(Surf)),t(Surf), type= 'l', lty = 1, 
		col = rev(Ramp(nrow(Surf)+5)[-c(1:5)]), lwd=2, axes=FALSE,xlab="",ylab="")
segments(min(ages),0,max(ages),0)
segments(ages[ages%%5==0],0,ages[ages%%5==0],-.01)
text(ages[ages%%5==0],-.01,ages[ages%%5==0],pos=1)
segments(min(ages),0,min(ages),.8)
segments(min(ages),seq(0,.8,by=.2),min(ages)-.2,seq(0,.8,by=.2))
text(min(ages)-.2,seq(0,.8,by=.2),seq(0,.8,by=.2),pos=2)
text(85,-.08,"Age",cex=1.5)
text(68,.4,"Prevalence",srt=90,cex=1.5)
dev.off()


