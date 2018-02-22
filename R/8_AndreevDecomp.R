# written by Alyson van Raalte
############################################
# TR: Note, this script uses almost the entire HMD, and it is hard coded to 
# run from the unzipped database, which you can download from here:
# http://www.mortality.org/cgi-bin/hmd/hmd_download.php,
# 'All HMD statistics'
# unzipping will create a folder call /hmd_statistics/.
# create a path like "C:\\hmd_statistics" to give to the decomposition function defined below
# in order to 
# Code for Figure 5 - Decomp results

library(plyr)
library(dplyr)
options(scipen=4,stringsAsFactors=FALSE)


# TR: this can be reworked to make use of object TTDprev, created in 
# the script 3_ApplicationPrep.R
TTDdf = read.csv("C:\\Dropbox\\sex decomp Tim\\postWP\\TTDdf.csv",head=T)


Ctry="SWE"
Year1=1970
Year2=1980
sex="m"
Morbtype="adl5_2"

# Decomposition using Andreev et al. formula
Decomp.DLY <- function(Ctry,Year1,Year2,sex,Morbtype,hmdpath = "C:\\hmd_statistics"){
  # Ctry is HMD short form, sex is 'm' or 'f', Morbtypes can be seen unique(meanprev$Morbidity)
  Morbdat <- filter(TTDdf,varname==Morbtype,sex==sex)
  
  if(sex=="m"){
	  this.path <- file.path(hmdpath,"lt_male","mltper_1x1",paste0(Ctry,".mltper_1x1.txt"))
	 
  } else {
	  this.path <- file.path(hmdpath,"lt_female","fltper_1x1",paste0(Ctry,".fltper_1x1.txt"))
  }
  LT <- read.table(this.path, head = TRUE, skip = 2, na.strings = ".")
  
  LT$Age <- rep(0:110,length(unique(LT$Year)))
  
  LT1 <- filter(LT,Year==Year1,Age>=60)
  LT2 <- filter(LT,Year==Year2,Age>=60)
  
  # vector of remaining lifespans  
  ages <- 60:110
  rls1 <- list()   
  rls2 <- list()
  for (i in 1:length(ages)){
    rls1[[i]] <- (LT1$dx[LT1$Age>=ages[i]])/sum(LT1$dx[LT1$Age>=ages[i]])  
    rls2[[i]] <- (LT2$dx[LT2$Age>=ages[i]])/sum(LT2$dx[LT2$Age>=ages[i]])  
  }
  
  # adding NAs to make each vector the same length
  rlswNA1 <- lapply(rls1, function(x) {c(x, rep(NA, length(ages) - length(x)))})
  rlswNA2 <- lapply(rls2, function(x) {c(x, rep(NA, length(ages) - length(x)))})
  
  # putting into matrix format
  rlsmat1 <- do.call(rbind, rlswNA1)
  rlsmat2 <- do.call(rbind, rlswNA2)
  
  # disability prevalence by age (pie i)
  DP1 <- colSums(t(rlsmat1[,1:12])*Morbdat$prev[1:12],na.rm=T)
  DP2 <- colSums(t(rlsmat2[,1:12])*Morbdat$prev[1:12],na.rm=T)
  # changing zeros at oldest ages (caused by zero deaths) to max DP
  DP1[DP1==0] <- max(DP1)
  DP2[DP2==0] <- max(DP2)
  
  # proportions disabled at each age
  pi1 <- DP1
  pi2 <- DP2
  
  # important life table columns
  Lx1 = LT1$Lx ; Lx2=LT2$Lx
  lx1 = LT1$lx ; lx2=LT2$lx
  qx1 = LT1$qx ; qx2=LT2$qx
  
  # changing radix to 1
  Lx1=Lx1/lx1[1]
  Lx2=Lx2/lx2[1]
  lx1=lx1/lx1[1]
  lx2=lx2/lx2[1]
  
  # disability life years
  h1 <- rev(cumsum(rev(Lx1*pi1)))/lx1
  h2 <- rev(cumsum(rev(Lx2*pi2)))/lx2
  h1[is.na(h1)==T] <- 0
  h2[is.na(h2)==T] <- 0
  
  # Person-years divided by survivors (P in Andreev et al.)
  P1 <- Lx1/lx1
  P2 <- Lx2/lx2
  P1[is.na(P1)==T] <- 1 #doesn't really matter what I put here since it is at ages beyond max age at death
  P2[is.na(P2)==T] <- 1
  
  # Formulas 10 & 11 from Andreev et al.
  # Health and mortality components
  HC <- MC <- rep(0,length(ages))
  for(i in 1:(length(HC)-1)){
    MC[i] <- 0.25*(lx1[i]+lx2[i])*(P2[i]-P1[i])*(pi1[i]+pi2[i])+
      0.5*(h1[i+1]*lx2[i]+h2[i+1]*lx1[i])*(qx1[i]-qx2[i])
    HC[i] <- 0.25*(lx1[i]+lx2[i])*(P1[i]+P2[i])*(pi2[i]-pi1[i])
  }
  sum(MC+HC)
  sum(h2[1]-h1[1])
  
  res <- data.frame(Ctry=Ctry,Year1=Year1,Year2=Year2,total.diff=h2[1]-h1[1],mortality.comp=sum(MC),
                    disability.comp=sum(HC),
                    e60.T1=LT1$ex[1],e60.T2=LT2$ex[1],e60.change=LT2$ex[1]-LT1$ex[1],
                    DLY.T1=h1[1],DLY.T2=h2[1],HLY.T1=LT1$ex[1]-h1[1],HLY.T2=LT2$ex[1]-h2[1],
                    Sex=sex,Morbtype=Morbtype,error=sum(MC+HC)-sum(h2[1]-h1[1]))
  return(res)
}




#--------------running on all ctries

# Ctries <- read.csv("C:\\Dropbox\\sex decomp Tim\\Data\\country_codes_HMD.csv",head=T)
Ctries <- read.csv("Data/country_codes_HMD.csv",head=T)
Ctry <- Ctries$Ctry
# TR: redone here to remove hard-coded file. This list includes Greece, which is a new
# HMD country. Others may join in the future, exercise still valid as such.
#Ctry <- HMDHFDplus::getHMDcountries()

# within country comparison
Ctrydf <- list()
for(c in 1:length(Ctry)){
  Ctrydf[[c]] <- data.frame(Ctry=Ctries$Ctry[c],
                            Year1=seq(Ctries$fyear10[c],Ctries$lyear10[c]-10,10),
                            Year2=seq(Ctries$fyear10[c]+10,Ctries$lyear10[c],10))
}
Ctrydf <- rbind.fill(Ctrydf)
# TR: appears to be from:
# https://stackoverflow.com/questions/11693599/alternative-to-expand-grid-for-data-frames
expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))
varchoose <- expand.grid.df(Ctrydf,data.frame(Sex=c("f","m")),
                            data.frame(morbtype=c("adl5_1","adl5_2","adl5_3","iadl5_1","iadl5_2","iadl5_3","nh_now","name_mo","srhpoor")))
varchoose <- filter(varchoose,Year1>=1950)

do.this <- FALSE
if (do.this){
	Bigres <- list()
	for (i in 1:dim(varchoose)[1]){
		Bigres[[i]] <- Decomp.DLY(Ctry=varchoose[i,1],Year1=varchoose[i,2],Year2=varchoose[i,3],
				sex=varchoose[i,4],Morbtype=varchoose[i,5])
	}
	Bigres <- rbind.fill(Bigres)
	
	write.csv(Bigres,file="Data/Within_country_decomp.csv",row.names=F)
	
}


#Bigres <- list()
#for (i in 1:15){
#  Bigres[[i]] <- Decomp.DLY(Ctry=varchoose[i,1],Year1=varchoose[i,2],Year2=varchoose[i,3],
#                            sex=varchoose[i,4],Morbtype=varchoose[i,5])
#}
#Bigres <- rbind.fill(Bigres)
#Bigres
