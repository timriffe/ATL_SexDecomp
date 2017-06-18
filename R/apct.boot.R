
# Author: tim
# based on script from Maarten J. Bijlsma
###############################################################################
library(splines)
library(reshape2)
# TR: new selection of first cases

# utility
cutla <- function(newdata, year1 = 1992, year2 = 2011, maxl = 100){
	# cut age
	newdata$la <- newdata$ca + newdata$ta
	mini       <- newdata$ca >= (year1 - newdata$b_yr - 1)
	maxi       <- newdata$la < (year2 - newdata$b_yr)
	newdata    <- newdata[mini & maxi, ]
	newdata    <- newdata[newdata$la < maxl, ]
	newdata
}

# This is a general boot function for this particular dataset and method,
# takes care of some data management too. Edge oversampling is by default turned
# off, but can be controlled with YearsFromEdge and MagFactor

apct.boot <- function(
		data,
		varname = "adl3_",
		nboot = 250,
		t.age = 0:12,     
		c.age = 70:100,
		b_yr_range = min(data$b_yr):max(data$b_yr),
		YearsFromEdge = 0,          # for optional edge diagnostics
		MagFactor = 1) {
	
	# test:
	if (!"la_int" %in% colnames(data)){
		data$la_int <- floor(data$ta + data$ca)
	}
	
	# the most recent wave year
	year2 <- max(as.integer(format(as.Date(data$intv_dt, format="%d/%m/%Y"),"%Y")))
	
	
    # data is the HRS object, selected for dead people
	# and whatever other data prep was required, e.g. to
	# make vars binary if wished
	# varname is 
    
	# this line super important 
	data            <- data[order(data$id), ]
	
	# determine potential edge cases
	data$edgies   <- data$la_int > year2 - data$b_yr - YearsFromEdge
	if (YearsFromEdge == 0){
		data$edgies <- FALSE
	}
	
	# select first cases, for sampling
	ids          <- data$id
	lengths      <- rle(ids)$lengths
	select.first <- cumsum(lengths) - lengths + 1
	dataid       <- data[select.first, ]

	# another object useful for resampling
	rows         <- 1:nrow(data)
	rowlist      <- split(rows, ids)
		
	dataid$drawweight                	<- dataid$p_wt2
	
	# this either will or won't have an effect depending on edgies and mag factor
	# by default no edge oversampling
	dataid$drawweight[dataid$edgies] 	<- dataid$drawweight[dataid$edgies] * MagFactor
	
    # assign person weight of first observation to each subsequent observation
	data$firstweight                    <- rep(dataid$drawweight, lengths)
    # given repalcement selection by a draw weight, the remaining observations
	# of individuals are thusly reweighted
	data$rescaleweight                  <- data$p_wt2 / data$firstweight
	
	# reduce size of object to just required columns
	col.index.outc 	<- grep(paste0('^', varname, '$'), colnames(data))
	col.index.id 	<- grep(paste0('^', 'id', '$'), colnames(data))
	col.index.b_yr 	<- grep(paste0('^', 'b_yr', '$'), colnames(data))
	col.index.ta 	<- grep(paste0('^', 'ta', '$'), colnames(data))
	col.index.ca 	<- grep(paste0('^', 'ca', '$'), colnames(data))
	col.index.rscw 	<- grep(paste0('^', 'rescaleweight', '$'), colnames(data))
	colselect       <- c(col.index.outc, col.index.id, col.index.rscw, col.index.b_yr, 
			             col.index.ta, col.index.ca, col.index.rscw)
	data            <- data[, colselect]
	gc()
	
	# and for those that are also edgies we need to weight in the opposite direction
	# same magnitude (but there will be more such people)
	# this line redundant, since we now divide by drawweight (by default equal to p_wt2)
    # Dat$rescaleweight[Dat$edgies]       <- Dat$rescaleweight[Dat$edgies] / MagnificationFactor
	
	# cut down dataid object to two needed columns, id and drawweight
	dataid     <- dataid[,c("id","drawweight")]
	
	
	# data for prediction, grid
	newdata    <- expand.grid(ta = t.age+.5, 
			                  ca = c.age+.5, 
			                  b_yr = b_yr_range)
	
	# remove extrapolation points for glm prediction
	out        <- cutla(newdata, year1 = 1992, year2 = year2)
	
	# number of cells in the jacked up Lexis surface that we actually need estimates for
	ncell      <- nrow(out)
	
	# matrix in which we save our estimates
	bootsave   <- matrix(ncol = ncell, nrow = nboot)
	idlength   <- nrow(dataid)

	for(b in 1:nboot) {
		# draw IDs with replacement and with weight
		selectid    <- sample(dataid[, 1],
				           size = idlength,
						   replace = TRUE,
				           prob = dataid[,2])
		
		select.rows <- unlist(rowlist[as.character(selectid)])
		data.boot   <- data[select.rows, ]
	
		## perform the ns code
		col.index   <- grep(paste0('^', varname,'$'), colnames(data.boot))

		if(max(data.boot[, col.index], na.rm = TRUE) == 1) {
			# determine if data is binary or not
			# if binary do (quasi)binomial glm
			fit        <- glm(data.boot[, col.index] ~ 
								ns(b_yr, knots = seq(1902.5, 1925.5, by = 5)) + 
								ns(ta, knots = c(.5, 1, 2, 4, 7.5, 10)) +  
								ns(ca, knots = seq(72.5, 97.5, by = 5)), 
							  data = data.boot,
							  weights = rescaleweight,
							  family = quasibinomial)
			
		} else { # if not binary do lr glm
			fit        <- glm(data.boot[, col.index] ~ 
								ns(b_yr, knots = seq(1902.5, 1925.5, by = 5)) + 
								ns(ta, knots = c(.5, 1, 2, 4, 7.5, 10)) +  
								ns(ca, knots = seq(72.5, 97.5, by = 5)), 
							  data = data.boot,
							  weights = rescaleweight)
		}
		
		# easier to keep dimensions straight if we predict over rectangular grid, 
		# then throw out values outside range
		# output so that bootstrap function can ...bootstrap it
		bootsave[b, ] <- predict(fit, out, type = 'response')
		
		print(b)
		
	}
	# don't need last estimate
	out$pi      <- NULL
	# change centroids to lower bounds for easier plotting
	out$ta     	<- floor(out$ta)
	out$ca     	<- floor(out$ca)
	
	# return both vec and it's named dimensions
	return(list(boot.est = bootsave,
					dims = out))
	
}
# utility, called in get.goods()
boot.ci <- function(bootdata,conf.level) {
	# bootdata is a dataframe or matrix like bootsave
	# conf.level is the confidence level we want
	# e.g. a 95% CI, a 99% CI, etc.
	
	# CI bounds
	lower <- (1-conf.level)/2
	upper <- conf.level+lower
	
	return(apply(bootdata,2,quantile,probs=c(lower,upper)))
	
}
# mean, median, 2.5% and 97.5% quantiles
get.booty <- function(boot.list){
	# separate parts
	dims        <- boot.list$dims
	boot.est    <- boot.list$boot.est
	dims$mean   <- apply(boot.est, 2, median, na.rm = TRUE)
	dims$median <- apply(boot.est, 2, mean, na.rm = TRUE)
	dims        <- cbind(dims, t(boot.ci(boot.est, .95)))
	dims$Width  <- c(diff(t(dims[,c("2.5%","97.5%")])))
	dims
}
# turn any of the above columns into a TTD by age matrix for a selected cohort
get.mat <- function(goods, column = "median", cohort = 1915){
	acast(goods[goods$b_yr == cohort,], ta ~ ca, value.var = column)
}
# order of call = apct.boot, get.booty, get.array??
get.array <- function(goods, column = "median"){
	acast(goods, ta ~ ca ~ b_yr, value.var = column)
}


apct.boot.wrapper <- function(
		Dat,
		varname = "adl3",
		sex = "f",
		t.age = 0:12,     
		c.age = 70:100,
		b_yr_range = min(data$b_yr):max(data$b_yr),
		nboot = 250){
	
	step1 <- apct.boot(data = Dat[sex == sex,],
			varname = varname,
			t.age = t.age,     
			c.age = c.age,
			b_yr_range = min(data$b_yr):max(data$b_yr),
			nboot = nboot,
			YearsFromEdge = 0,          
			MagFactor = 1)
	# creates an array of the median estimate.
	# we throw away all the rest. Could adapt
	# this to create an array for a set of quantiles..
	Surf <- get.array(get.booty(step1))		
	
	list(Surf = Surf, sex = sex, varname = varname)
}



