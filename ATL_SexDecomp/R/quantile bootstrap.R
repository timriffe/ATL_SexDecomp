
#if (system("hostname",intern=TRUE) %in% c("triffe-N80Vm", "tim-ThinkPad-L440")){
#	# if I'm on the laptop
#	setwd("/home/tim/git/APCT/APCT")
#} else {
#	if (system("hostname",intern=TRUE) == "PC-403478"){
#		# on MPIDR PC
#		setwd("U://git//APCT//APCT")
#	} else {
#		# in that case I'm on Berkeley system, and other people in the dept can run this too
#		setwd(paste0("/data/commons/",system("whoami",intern=TRUE),"/git/APCT/APCT"))
#	}
#}	
# select this first!!
#varnamechar <-"adl3_"
varnamechar <-"psych"


setwd("/home/tim/git/ATL_SexDecomp/ATL_SexDecomp")
# load data
# and other data prep
#{setwd('U:/Collaboration/TR AVR')
 load("Data/Data_long.Rdata")
 
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
               "iadl_calc", "mprob", "mprobev", "med_explog") # this is bigger than the final list...
 # pare down to columns available in current iteration
 varnames <- varnames[varnames %in% colnames(Dat)]
 
 # cut data down, define time vars
 
 # no ages under 65
 Dat      <- Dat[Dat$age >= 65, ]
 
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
 
 ## loess curves
 Dat$ta_int <- floor(Dat$ta)
 Dat$ca_int <- floor(Dat$ca)
 Dat$la_int <- floor(Dat$ta + Dat$ca)
 
 # let's remove people with ta = -1
 Dat <- Dat[Dat$ta_int > -1,]
 
 # crude but effective function for taking means of residuals
 # and plotting them
 # mostly just to save space in the code
 Datm <- Dat[Dat$sex=='m',]
 Datf <- Dat[Dat$sex=='f',]
 
 # set knots
 nk <- 7
 
 # give row numbers (rows are already ordered chronologically 
 # for this data)
 # and give row weights
 Dat$rownr <- NA
 Dat$rescaleweight <- 0
 for(i in unique(Dat$id)) {
   
   rownr                <- sum(Dat$id==i)
   Dat$rownr[Dat$id==i] <- 1:rownr
   
   drawweight <- Dat$p_wt2[Dat$id==i & Dat$rownr==1]
   Dat$rescaleweight[Dat$id==i] <- Dat$p_wt2[Dat$id==i]/drawweight
   
   #print(i)
 }
 
 # determine weights
 # draw weights
 sumdraw <- sum(Dat$p_wt2[Dat$rownr==1])
 Dat$drawprob <- 0
 Dat$drawprob[Dat$rownr==1] <- Dat$p_wt2[Dat$rownr==1]/sumdraw
 
 cutla <- function(newdata, year1 = 1992, year2 = 2011, maxl = 100){
   # cut age
   newdata$la <- newdata$ca + newdata$ta
   mini       <- newdata$ca >= (year1 - newdata$b_yr - 1)
   maxi       <- newdata$la < (year2 - newdata$b_yr)
   newdata    <- newdata[mini & maxi, ]
   newdata    <- newdata[newdata$la < maxl, ]
   newdata
 }
 


id.dat <- data.frame(cbind(Dat$id[Dat$rownr==1],Dat$drawprob[Dat$rownr==1]))
# id.dat = dataid
# Dat = dataset

t.age <- 0:12
c.age <- 70:100
b_yr_range <- min(Dat$b_yr):max(Dat$b_yr)

# choose variable
varnamechar <- varnames[1]

# make a smaller dataset for the loop/bootstrap
col.index.outc 	<- grep(paste0('^',varnamechar,'$'), colnames(Dat))
col.index.id 	<- grep(paste0('^','id','$'), colnames(Dat))
col.index.b_yr 	<- grep(paste0('^','b_yr','$'), colnames(Dat))
col.index.ta 	<- grep(paste0('^','ta','$'), colnames(Dat))
col.index.ca 	<- grep(paste0('^','ca','$'), colnames(Dat))
col.index.rscw 	<- grep(paste0('^','rescaleweight','$'), colnames(Dat))

Dat.2 <- Dat[,c(col.index.outc,col.index.id,col.index.rscw,
                col.index.b_yr, col.index.ta,col.index.ca)]

# premake prediction data:
newdata    <- expand.grid(ta = t.age+.5, 
		ca = c.age+.5, 
		b_yr = b_yr_range)

# remove extrapolation points for glm prediction
out        <- cutla(newdata, year1 = 1992, year2 = 2011)

# make a matrix for saving bootstrap estimates in
# in this case, we happen to know that we need
# 3950 columns (1 for each estimate cell in the Lexis diagram 
# and bootstrap size number of rows

# bootstrap function
apct.boot <- function(
		dataid,
		data,
		bootlength,
		t.age = 0:12,       # TR: added these args that are otherwise only coincidentally visible
		c.age = 70:100,
		b_yr_range = min(data$b_yr):max(data$b_yr)) {
  # dataid is just a dataframe with ids and draw weights
  # we bootstrap this way so that we effectively
  # do block resampling
  # data should be the dataset we are bootstraping
  # ideally it only contains the necessary columns so that the
  # bootstrap is more efficient
  # those are: 
  # 1. The outcome variable of interest
  # 2. ID variables (with draw weights)
  # 3. the APCT variables of interest
  # 4. longitudinal weights
  # bootlength is a number telling us how many bootstrap
  # iterations it should do
  
  # data for prediction, grid
  newdata    <- expand.grid(ta = t.age+.5, 
			ca = c.age+.5, 
			b_yr = b_yr_range)
	
  # remove extrapolation points for glm prediction
  out        <- cutla(newdata, year1 = 1992, year2 = 2011)
	
  # number of cells in the jacked up Lexis surface that we actually need estimates for
  ncell    <- nrow(out)
  # matrix in which we save our estimates
  bootsave <- matrix(rep(NA,ncell*bootlength),nrow=bootlength)
  idlength <- nrow(dataid)
  # b <- 1
  for(b in 1:bootlength) {
    
    # draw IDs with replacement and with weight
    selectid <- sample(dataid[,1],size=idlength,replace=TRUE,
                       prob=dataid[,2])
    
    # see how often an ID appears
    idtab <- sort(table(selectid))
    
    # save the names of the IDs (read: their position)
    idtabname <- as.numeric(names(idtab))
    
    # save the data of the people who are selected once
    id.once <- idtabname[idtab==1]
    
    data.boot <- data[which(data$id %in% id.once),]
    
    ## then save data for people who are selected more than once
    # how many times are people selected?
    selecttimes <- unique(idtab)
    selecttimes <- selecttimes[selecttimes != 1]
    

    for(st in selecttimes) {
      
      # select the ids that were selected randomly
      # 'st number of times
      id.mult <- idtabname[idtab==st]
      
      # save the data of these people
      temp <- data[which(data$id %in% id.mult),]
      
      # and copy it 'st' number of times
      temp <- temp[rep(1:length(temp$id),each=st),]
      
      # and append it to data.boot
      data.boot <- rbind(data.boot,temp)
      
      # just to be sure
      temp <- NULL
      
    }
    
    ## perform the ns code
    col.index <- grep(paste0('^',varnamechar,'$'), colnames(data.boot))
    
    # remove all columns that we don't need
    # first we find the indexes of the variables we do need
    # (ideally this is done outside the loop for more efficiency)
    
    if(max(data.boot[, col.index], na.rm = TRUE) == 1) {
      # determine if data is binary or not
      # if binary do (quasi)binomial glm
      fit        <- glm(data.boot[, col.index] ~ 
                          ns(b_yr, knots = seq(1902.5,1925.5,by=5)) + 
                          ns(ta, knots = c(.5,1,2,4,7.5,10)) +  
                          ns(ca, knots = seq(72.5,97.5,by=5)), 
                        data = data.boot,
                        weights = rescaleweight,
                        family = quasibinomial)
      
    } else { # if not binary do lr glm
      fit        <- glm(data.boot[, col.index] ~ 
                          ns(b_yr, knots = seq(1902.5,1925.5,by=5)) + 
                          ns(ta, knots = c(.5,1,2,4,7.5,10)) +  
                          ns(ca, knots = seq(72.5,97.5,by=5)), 
                        data = data.boot,
                        weights = rescaleweight)
    }
    

    
    # easier to keep dimensions straight if we predict over rectangular grid, 
    # then throw out values outside range
    out$pi     <- predict(fit, out,type='response')
    # length(out$pi); ncol(bootsave)
    # output so that bootstrap function can ...bootstrap it
    bootsave[b,] <- out$pi
    
    print(b)
    
  }
  # don't need last estimate
  out$pi <- NULL
  
  # return both vec and it's named dimensions
  return(list(boot.est = bootsave,
			  dims = out))

}

bootlength  <- 500
boot.list   <- apct.boot(id.dat,Dat.2,bootlength)
boot.list3  <- apct.boot(dataid3,Dat.3,bootlength)


hist(boot.list3$boot.est[,which(with(dims, b_yr==1910 & ta == 10 & la == 98))],breaks=20,freq=F)
hist(boot.list$boot.est[,which(with(dims, b_yr==1910 & ta == 10 & la == 98))],freq=F,breaks=20,add=TRUE,col = "#00000050")



head(dims)
boot.est    <- boot.list$boot.est
dims        <- boot.list$dims

Median      	<- apply(boot.list$boot.est,2,median,na.rm=T)
Median3      	<- apply(boot.list3$boot.est,2,median,na.rm=T)
# what are the coefficients?
# get the mean and medians
Mean        	<- apply(boot.list$boot.est,2,mean,na.rm=T)
Mean3        	<- apply(boot.list3$boot.est,2,mean,na.rm=T)
colnames(Dat)
# get the comparison default thingy:
fit         	<- glm(adl3_ ~ 
						ns(b_yr, knots = seq(1902.5,1925.5,by=5)) + 
						ns(ta, knots = c(.5,1,2,4,7.5,10)) +  
						ns(ca, knots = seq(72.5,97.5,by=5)), 
				       weights = p_wt2,
					   data = Dat,
					   family = quasibinomial)
Comparison     	<- predict(fit, out,type='response')

# assign them to the dims
dims$Mean   	<- Mean
dims$Median 	<- Median
dims$Mean3 	    <- Mean3
dims$Median3 	<- Median3
dims$Comparison <- Comparison

# adjust the dims a little bit
dims$ta     	<- floor(dims$ta)
dims$ca     	<- floor(dims$ca)

# make matrices for plotting:
CompareMat <- acast(dims[dims$b_yr == 1915,], ta~ca, value.var = "Comparison")
MeanMat    <- acast(dims[dims$b_yr == 1915,], ta~ca, value.var = "Mean")
MedMat     <- acast(dims[dims$b_yr == 1915,], ta~ca, value.var = "Median")
MedMat3     <- acast(dims[dims$b_yr == 1915,], ta~ca, value.var = "Median3")
library(reshape2)

args(SurfMap)
source("R/SurfMap.R")
par(mfrow=c(1,2))
#SurfMap(MeanMat)
SurfMap(MedMat)
title("Boot median")
SurfMap(MedMat3)
title("Edgy boot median")
SurfMap(CompareMat)
title("default glm")

range(MedMat3-MedMat,na.rm=TRUE)


library(LexisUtils)
range(MeanMat - MedMat, na.rm=TRUE)
range(CompareMat - MedMat, na.rm=TRUE)
range(CompareMat - MeanMat, na.rm=TRUE)
range(MedMat3 - MeanMat, na.rm=TRUE)
range(MedMat3 - CompareMat, na.rm=TRUE)
range(booty3_97 - booty_97, na.rm=TRUE)

range(t(MedMat3 - CompareMat)/t(CompareMat),na.rm=TRUE)
mean(MedMat3 - CompareMat, na.rm=TRUE)


booty3_97    <- acast(dims[dims$b_yr == 1915,], ta~ca, value.var = "booty3_97.5%")
booty_97    <- acast(dims[dims$b_yr == 1915,], ta~ca, value.var = "booty_97.5%")
range(t(booty3_97 - booty_97)/t(booty_97),na.rm=TRUE)


library(RColorBrewer)
breaks <- seq(-.1,.1,length=11)
breaks <- seq(-.07,.07,length=11)
cols   <- colorRampPalette(brewer.pal(9,"RdBu"),space="Lab")
image(t(booty3_97 - booty_97)/t(booty_97), breaks = breaks, col = rev(cols(length(breaks)-1)))
contour(t(booty3_97 - booty_97)/t(booty_97),breaks=breaks,add=TRUE)
# comparison glm()



sqrt(36^2 + 48^2)
unique(Dat$psych)

image(CompareMat)



hist(boot.est[,2000])

(50 * 10) / 60

# quantile CI function
boot.ci <- function(bootdata,conf.level) {
  # bootdata is a dataframe or matrix like bootsave
  # conf.level is the confidence level we want
  # e.g. a 95% CI, a 99% CI, etc.
  
  # CI bounds
  lower <- (1-conf.level)/2
  upper <- conf.level+lower
  
  return(apply(bootdata,2,quantile,probs=c(lower,upper)))
  
}
dims <- cbind(dims, t(boot.ci(boot.list$boot.est,0.95)))
dims <- cbind(dims, t(boot.ci(boot.list3$boot.est,0.95)))
#dput(colnames(dims))
colnames(dims) <- c("ta", "ca", "b_yr", "la", "Mean", "Median", "Mean3", "Median3", 
		"Comparison", "booty_2.5%", "booty_97.5%", "booty3_2.5%", "booty3_97.5%")


# test:


devtools::load_all("/home/tim/git/TR1/TR1/HMDHFDplus")
flt <- readHMDweb("USA","fltper_1x1",username = us, password = pw)
mxf <- acast(flt, Age~Year, value.var = "mx")
LexisUtils::LexisMap(mxf)
abline(v=1995)
plot(apply(mxf,2,which.min) - 1)

plot(0:110, mxf[,"1995"], log = 'y',type = 'l')
lines(0:110, mxf[,"1996"], col = "red")
matplot(0:110,mxf[,as.character(1992:1997)], type = 'l', 
		log = 'y', col = RColorBrewer::brewer.pal(8,"Blues")[-c(1,2)], lty=1)
Diffs <- t(diff(t(mxf[,as.character(1992:1997)])))
matplot(0:110,Diffs, type = 'l', 
		
		col = RColorBrewer::brewer.pal(8,"Blues")[-c(1,2,3)], lty=1, ylim = c(-.001,.001))
abline(h=0,col="red")
lines(0:110, Diffs[,"1995"], col = "magenta")
abline(v=10)



matplot(0:110,Diffs/mxf[,as.character(1992:1996)], type = 'l', 
		col = RColorBrewer::brewer.pal(8,"Blues")[-c(1,2,3)], lty=1)
abline(h=0,col="red")
lines(0:110, Diffs[,"1995"]/mxf[,as.character(1995)], col = "magenta")
abline(v=10)
