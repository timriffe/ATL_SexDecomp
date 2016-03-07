
# Author: tim
###############################################################################


# make new ATL surface plot.

Dat <- local(get(load("/home/tim/git/ATL_SexDecomp/ATL_SexDecomp/Data/SmArrays.Rdata")))
names(Dat)
library(RColorBrewer)
ramp <-
names()
X <- Dat[["SRHpoor"]]$Females$Surf[,,3]

source("/home/tim/git/ATL_SexDecomp/ATL_SexDecomp/R/SurfMap.R")
ticks <- seq(0,.8,by=.1)
graphics.off()
dev.new(width=7,height=4.2)
pdf("Figures/SRHpoor_Females.pdf",width=7,height=4.2)
SurfMap(Dat[["IADL1"]]$Females$Surf[,,3], ticks = ticks, bg=TRUE,ylab="time-to-death",xlab="age",
		mai=c(.6,.5,.4,1.2))
text(103,16,"prevalence")
dev.off()
pdf("Figures/IADL1_Males.pdf",width=7,height=4.2)
SurfMap(Dat[["IADL1"]]$Males$Surf[,,3], ticks = ticks, bg=TRUE,ylab="time-to-death",xlab="age",
		mai=c(.6,.5,.4,1.2))
text(103,16,"prevalence")
dev.off()

names(Dat)