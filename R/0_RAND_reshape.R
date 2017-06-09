
# Author: tim
###############################################################################

setwd("/home/tim/git/HLETTD")

rm(list = ls())
gc()
data.path <- "/home/tim/Data/HRS/HRSvP"
# location of RAND data product
library(foreign)
library(stringr)

# read in full HRS rand data
temp                <- read.dta(file.path(data.path,"Data/Original/RAND/rndhrs_p.dta") )
# column names manually grabbed by MG for PHC's script
colsRAND            <- read.csv(file.path(data.path,"Data/namechangesHRS.csv"))

vary_want           <- c("inw", colsRAND$HRSold)
# rm bwc86
vary_want           <- vary_want[vary_want != "bwc86"]
# how many waves in present RAND version?
waves               <- 1:12
waves2              <- sprintf("%02d", waves)
# the time-invariant variables that we keep perhaps for no reason
invariant_keep		<- c("hhidpn","hacohort","rawtsamp","rabmonth",
		"rabyear","rabdate","radmonth","radyear",
		"raddate","ragender","raracem","rahispan",
		"raedyrs","rabplace","rameduc","rafeduc")
# some overlap...
vary_want           <- vary_want[!vary_want %in% invariant_keep]

# ------------------------------- #
# preamble functions              #
# ------------------------------- #

# intuit time-varying
varies <- function(x,waves=1:12){
	xx <- strsplit(x, "[^0-9]+")[[1]]
	if (length(xx) >= 2){
		dig <- as.integer(xx[2])
		return(dig %in% waves)
	}
	FALSE
}

# --------------------------------
# first character of root names, must be r
varies_r <- function(x,waves=1:12){
	varies(x,waves) & substring(x,1,1) == "r"
}

# --------------------------------
# get root name of a time-varying column name
get_root <- function(x){
	from <- stringr::str_locate_all(pattern ="[^0-9]", x)[[1]][2,1]
	to   <- nchar(x)
	substring(x,from,to)
}

# --------------------------------
# rearrange time-varying column names to something
# reshape will correctly guess the meaning of
# steps: 1 remove leading character, r
#        2 separate first digits
#        3

rearrange_varying <- function(x){
	# wave digit
	dig <- strsplit(x, "[^0-9]+")[[1]][2]
	# make two-digit, leading 0 if necessary
	dig <- ifelse(nchar(dig)==1,paste0(0,dig),dig)
	# compose palatable name:
	paste(get_root(x), dig, sep = ".")
}

# -------------------------------- #
# column name manipulations        #
# -------------------------------- #

# get all names
x                <- colnames(temp)
# first let's rename inwn to rninw ... to be consistent w string processing...
inws             <- x[grepl(x,pattern="inw")]
inw_new          <- paste0("r",waves,"inw")
colnames(temp)[grepl(x,pattern="inw")] <- inw_new
# get all names again...
x                <- colnames(temp)

# intuit all varying columns whose names start with r
varying          <- sapply(x,varies_r)

# some invariant columns to remove
invariant_ind    <- x %in% invariant_keep
# sum(invariant_all);length(invariant_keep) # 16

# get all root names of time-varying columns
varying_all      <- x[varying]
roots_all        <- sapply(varying_all,get_root)
# esoteric, yes, but we want the reference to be the column
# positions of the full data object, not the subset..., and can
# only match these on integer position...
varying_keep     <- 1:ncol(temp) %in% which(varying)[roots_all %in%	vary_want]

# -------------------------------- #
# cut columns down                 #
# -------------------------------- #
temp             <- temp[, invariant_ind | varying_keep]

# start again clean
# rearrange names of the varying columns. Kinda tricky
y                <- colnames(temp)
varying          <- sapply(y,varies)
# overwrite time-varying colnames, where applicable
y[varying]       <- sapply(y[varying], rearrange_varying)
# overwrite colnames of data object
colnames(temp)   <- y
# rename, just because
wide             <- temp
rm(temp);gc()
# -------------------------------- #
# fill in missing columns          #
# -------------------------------- #
# let's take care of colnumn names containing other numbers. Examples include 20,40,75,86.
# replace w names...
colnames(wide)   <- gsub(colnames(wide),pattern="20\\.",replacement="twenty\\.")
colnames(wide)   <- gsub(colnames(wide),pattern="40\\.",replacement="fourty\\.")
colnames(wide)   <- gsub(colnames(wide),pattern="75\\.",replacement="sevfiv\\.")
colnames(wide)   <- gsub(colnames(wide),pattern="10\\.",replacement="ten\\.")
# now get back to business
z                <- colnames(wide)
varying          <- sapply(colnames(wide),varies)
# unique root names
roots_all        <- unique(unlist(lapply(strsplit(z[varying],split="\\."),"[[",1)))
# expansion thereof
expansion        <- c(t(outer(roots_all,waves2,function(r,w){
							paste(r,w,sep=".")
						})))
# add decoy columns
varying_decoy    <- expansion[!expansion %in% z[varying]]
wide[varying_decoy]		<- NA
# reuse z..
z                <- colnames(wide)
varying          <- grepl(z,pattern="\\.")

# put time-invariant in front:
wide_inv         <- wide[,!varying]
wide_varying     <- wide[,varying]
# alphabetize varying: (sticks root names together, but integers not handled properly
# this is why we stuck 0s in front of waves, so 10,11,12 aren't at the front..
wide_varying     <- wide_varying[,sort(colnames(wide_varying))]
# stick back together:
wide             <- cbind(wide_inv, wide_varying)
# now remove those pesky 0s that were necessary for alphabetizing...
colnames(wide)   <- gsub(colnames(wide),pattern="\\.0",replacement="\\.")


# -------------------------------- #
# reshape to long                  #
# -------------------------------- #

# fix the varying varnames too
expansion <- sort(expansion)
expansion <- gsub(expansion,pattern="\\.0",replacement="\\.")

# the reshape call
long <- reshape(
		wide,
		idvar = "hhidpn",
		direction = "long",
		varying = expansion,
		times = list(waves))
# order by waves within individual
long <- long[order(long$hhidpn,long$time), ]
# NA weights have no time-varying variables
long <- long[!is.na(long$wtresp), ]

# -------------------------------- #
# Some final details               #
# -------------------------------- #
# dput(colnames(long)[!colnames(long) %in% colsRAND[,2]])
# dput(colsRAND[,2][!colsRAND[,2] %in% colnames(long) ])
have <- c("time", "adlrten", "aimrten", "atrtwenty", "bwctwenty", "dlrctwenty", 
		"hdlrtwenty", "himrtwenty", "htrfourty", "imrctwenty", "livsevfiv", 
		"trfourty", "trtwenty")
need <- c("wave", "adlr10", "aimr10", "atr20", "bwc20", "dlrc20", 
		"hdlr20", "himr20", "htr40", "imrc20", "liv75", 
		"tr40", "tr20")
colnames(long)[colnames(long) %in% have] <- need

# ----------------------------------------------------- #
# now rename once more to legacy (ThanoEmpirical) names #
# ----------------------------------------------------- #
x                   <- colnames(long)                   

# a vector of column names, whose names are identical 
# to values, for named selection
xnew                <- x
names(xnew)         <- x
# a recode vector, based on MG's table
recvec              <- colsRAND[,3]
names(recvec)       <- colsRAND[,2]
recvec              <- recvec[names(recvec) %in% x]
# now change names, maintaining order
xnew[names(recvec)] <- recvec
#cbind(x,xnew) # quick check ordering maintained

# now rename columns to downstream used names.
colnames(long)      <- xnew

# save out:
save(long,file="Data/RAND_p_long.Rdata")
