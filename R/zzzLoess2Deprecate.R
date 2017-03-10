
# Author: tim
###############################################################################

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
	# hist(Dat$p_wt2)
	# multiplicative give the most freedom.
	Dats <- Dat[Dat$sex == sex, ]
	
	mod   <- loess(paste0(varname,'~b_yr * ta * ca') ,
			data = Dats, 
			weights = p_wt2, # this, plus point density both act as weights
			span = .5,     # a variable passed in, or smoothness
			# is similiar conceptually to a 1:1:1 aspect ratio. Everything is in years...
			normalize = FALSE,
			family = "gaussian",
			control = loess.control(trace.hat="approximate")
	)
	
	newdata        <- expand.grid(ta = t.age+.5, ca = c.age+.5, b_yr = .Coh5)
	# easier to keep dimensions straight if we predict over rectangular grid, 
	# then throw out values outside range
	#newdata        <- newdata + .5
	Surf           <- predict(mod, newdata)
	
	dimnames(Surf) <- list(floor(t.age),floor(c.age), .Coh5)
	
	# need to trim on the left side where applicable, since some questions didn't enter until
	# wave 2 or 3. There are many such cases, so check and make sure we dont' extrapolate.
	# sex <- "m"
	# varname <- "cesd"
	MissingWaves <- tapply(Dats[Dats$sex==sex,varname],Dats[Dats$sex==sex,"wave"], function(x){
				all(is.na(x))
			})
	LeftYear <- 1992
	RightYear <- 2013
	if (any(MissingWaves)){
		Waves <- which(MissingWaves)
		# if it's 1 or two, trim, if it's more, then return NULL for now. Be conservative.
		if (any(Waves %in% 3:9)){
			return(NULL) # takes care of gaps. Also takes care of late missing waves.
		}
		if (any(Waves %in% c(1,2))){
			WaveGetYr <- max(Waves[Waves < 3]) + 1
			LeftYear  <- round(mean(as.numeric(format(Dat$intv_dt[Dat$wave == WaveGetYr], "%Y"))))
		}
		if (any(Waves == 10)){
			RightYear  <- round(mean(as.numeric(format(Dat$intv_dt[Dat$wave == 9], "%Y"))))
		}
		
	}
	# this reduces extrapolation outside of data points 
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


FitLoessAP <- function(varname, 
		Dat, 
		sex,
		c.age = 70:100,  # standard matrix size, though we may NA certain unobserved cells
		span = .5, # will vary
		.Coh5){
	# conservative here to cut tails
	maxL  <- 100
	minL  <- 70
	# hist(Dat$p_wt2)
	# multiplicative give the most freedom.
	Dats <- Dat[Dat$sex == sex, ]
	
	mod   <- loess(paste0(varname,'~b_yr * ta * ca') ,
			data = Dats, 
			weights = p_wt2, # this, plus point density both act as weights
			span = .5,     # a variable passed in, or smoothness
			# is similiar conceptually to a 1:1:1 aspect ratio. Everything is in years...
			normalize = FALSE,
			family = "gaussian",
			control = loess.control(trace.hat="approximate")
	)
	
	newdata        <- expand.grid(ta = t.age+.5, ca = c.age+.5, b_yr = .Coh5)
	# easier to keep dimensions straight if we predict over rectangular grid, 
	# then throw out values outside range
	#newdata        <- newdata + .5
	Surf           <- predict(mod, newdata)
	
	dimnames(Surf) <- list(floor(t.age),floor(c.age), .Coh5)
	
	# need to trim on the left side where applicable, since some questions didn't enter until
	# wave 2 or 3. There are many such cases, so check and make sure we dont' extrapolate.
	# sex <- "m"
	# varname <- "cesd"
	MissingWaves <- tapply(Dats[Dats$sex==sex,varname],Dats[Dats$sex==sex,"wave"], function(x){
				all(is.na(x))
			})
	LeftYear <- 1992
	RightYear <- 2013
	if (any(MissingWaves)){
		Waves <- which(MissingWaves)
		# if it's 1 or two, trim, if it's more, then return NULL for now. Be conservative.
		if (any(Waves %in% 3:9)){
			return(NULL) # takes care of gaps. Also takes care of late missing waves.
		}
		if (any(Waves %in% c(1,2))){
			WaveGetYr <- max(Waves[Waves < 3]) + 1
			LeftYear  <- round(mean(as.numeric(format(Dat$intv_dt[Dat$wave == WaveGetYr], "%Y"))))
		}
		if (any(Waves == 10)){
			RightYear  <- round(mean(as.numeric(format(Dat$intv_dt[Dat$wave == 9], "%Y"))))
		}
		
	}
	# this reduces extrapolation outside of data points 
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

