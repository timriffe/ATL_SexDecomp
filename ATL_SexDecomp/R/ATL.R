
# Author: tim
###############################################################################
setwd("/home/tim/git/ATL_SexDecomp/ATL_SexDecomp")

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

X <- Dat[["IADL1"]]$Males$Surf[,,3]
#pdf("Figures/IADL1_Males.pdf",width=7,height=4.2)
graphics.off()
#dev.new(width=7,height=4)
pdf("Figures/IADL1_Males_Surf_a.pdf",width=7,height=4)
SurfMap(X, ticks = ticks, bg=TRUE,ylab="",xlab="",
		mai=c(.6,.6,.1,.1),legnd = FALSE)
text(85,-2.5,"Age")
text(68, 7.5, "TTD",srt=90)
dev.off()
#display.brewer.all()
agecols <- colorRampPalette(brewer.pal(9,"PuBu")[-c(1:3)])
X2 <- X[,colSums(!is.na(X))>0]
#dev.new(width=7,height=4)
pdf("Figures/IADL1_Males_Age_c.pdf",width=7,height=4)
par(mai=c(.6,.6,.1,.1), xaxs="i",yaxs="i",xpd=TRUE)
plot(NULL, type = "n", xlim = c(70,101),ylim=c(0,.8),xlab="",ylab="",axes=FALSE)
matplot(as.integer(colnames(X2)),t(X2),
		add = TRUE,
		type='l',
		lty = 1,
		col = rev(agecols(ncol(X2))),
		lwd = seq(1.5,2.5,length=ncol(X2)))
# manual axes to match
segments(70,0,70,.8);segments(70,0,100,0)
segments(70,seq(0,.8,by=.2),69.6,seq(0,.8,by=.2))
segments(seq(70,100,by=5),0,seq(70,100,by=5),-.03)
text(69.6,seq(0,.8,by=.2),seq(0,.8,by=.2),pos=2,cex=.8)
text(seq(70,100,by=5),-.04,seq(70,100,by=5),pos=1,cex=.8)

# label y locations
yvals <- ncol(X2) - col(X2) + 1 == row(X2)
X2[yvals]
TTDs <- c(12:1,"TTD = 0")
keep <- c(TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,rep(TRUE,6))
text((83:96)[keep],X2[yvals][keep]-.01,TTDs[keep],pos=4,cex=.8)
text(85,-.12,"Age",xpd=TRUE)
text(67.5,.4,"Prevalence",srt=90)
dev.off()
# trick to get x dims to match.

ttdcols <- colorRampPalette(brewer.pal(9,"Greens")[-c(1:3)])
#dev.new(width=((7 - .7) / diff(c(70,101))) * 13 + .7,height=4)
pdf("Figures/IADL1_Males_TTD_b.pdf",width=((7 - .7) / diff(c(70,101))) * 13 + .7,height=4)
par(mai=c(.6,.6,.1,.1),xaxs="i",yaxs="i",xpd=TRUE)
plot(NULL, type = "n", xlim = c(0,13),ylim=c(0,.8),xlab="",ylab="",axes=FALSE)
matplot(as.integer(rownames(X2)),X2, 		
		type = 'l',
		lty = 1,
		add=TRUE,
		col = gray(.6),
		lwd = seq(1))
ages <- as.integer(colnames(X2))
keep <- ages %% 5 == 0
matplot(as.integer(rownames(X2)),X2[,keep], 		
		type = 'l',
		lty = 1,
		add=TRUE,
		col = ttdcols(sum(keep)),
		lwd = seq(3,2,length=sum(keep)))
segments(0,0,0,.8);segments(0,0,13,0)
segments(0,seq(0,.8,by=.2),-.4,seq(0,.8,by=.2))
segments(seq(0,13,by=5),0,seq(0,13,by=5),-.03)
text(-.4,seq(0,.8,by=.2),seq(0,.8,by=.2),pos=2,cex=.8)
text(seq(0,13,by=5),-.04,seq(0,13,by=5),pos=1,cex=.8)

yvals <- X2[,keep][c(13,26,37,45)]
xvals <- c(12,12,10,5)
text(xvals,yvals,c(75,80,85,"Age = 90"),pos=4,cex=.8)
text(6, -.12, "TTD")
text(-2.6,.4,"Prevalence",srt=90)
dev.off()

pdf("Figures/blankholder.pdf",width=((7 - .7) / diff(c(70,101))) * 13 + .7,height=4)
plot(NULL, type="n",axes=FALSE,xlim=c(0,2),ylim=c(0,2),xlab="",ylab="")
dev.off()
c(7,((7 - .7) / diff(c(70,101))) * 13 + .7) * 1.2 / 10