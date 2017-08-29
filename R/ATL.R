
# Author: tim
###############################################################################
# for TR computers...otherwise you need to set wd yourself
if (system("hostname",intern=TRUE) %in% c("triffe-N80Vm","tim-ThinkPad-L440")){
	# if I'm on the laptop
	setwd("/home/tim/git/HLETTD")
} 

library(RColorBrewer)
source("R/SurfMap.R")

# load list of smoothed variables
Results        <- local(get(load("Data/resultsP.Rdata")))
varnames       <- sapply(Results, "[[", "var")
names(Results) <- varnames
ResultsLong    <- local(get(load("Data/ResultsPLong.Rdata")))


# Male cesddepr 1915, Female srhpoor 1920


#X <- apply(Results[["cesddepr"]]$Male$Surf[,,c("1915","1916","1917","1918","1919")],
#		c(1,2),
#		function(x){
#			N <- sum(!is.na(x))
#			ifelse(N >= 2, mean(x, na.rm = TRUE), NA )
#		})
X <- apply(Results[["srhpoor"]]$Female$Surf[,,c("1920","1921","1922","1923","1924")],
		c(1,2),
		function(x){
			N <- sum(!is.na(x))
			ifelse(N >= 2, mean(x, na.rm = TRUE), NA )
		})
#X <- Results[["iadl5_1"]]$Male$Surf[,,"1915"]
ticks <- seq(0,.45,by=.05)
#pdf("Figures/IADL1_Males.pdf",width=7,height=4.2)
graphics.off()
#dev.new(width=7,height=5)
XS <- X
XS <- XS[rowSums(!is.na(XS)) > 0, colSums(!is.na(XS)) > 0]
pdf("Figures/srhpoor_f_Surf_a.pdf",width=7,height=5)
par(mai=c(.6, .6, .1, .1), xaxs = "i", yaxs = "i", xpd = TRUE)
plot(NULL, type = "n",xlim = c(70,91),ylim=c(0,13),axes=FALSE, xlab = "", ylab = "")
SurfMap(XS, ticks = ticks, bg=TRUE,ylab="",xlab="",
		mai=c(.6,.6,.1,.1),legnd = FALSE,add=TRUE)
text(80,-1.5,"Age", cex = 1.2)
text(68.5, 7, "TTD", cex = 1.2,srt=90)
dev.off()


# TODO: work out so x ticks match below and right.

# make the TTD over age lines
agecols <- colorRampPalette(brewer.pal(9,"PuBu")[-c(1:3)])
ages    <- as.integer(colnames(XS))
ymax    <- max(pretty(XS))
#dev.new(width=7,height=5)
pdf("Figures/srhpoor_f_Age_c.pdf", width = 7, height = 5)
par(mai=c(.6, .6, .1, .1), xaxs = "i", yaxs = "i", xpd = TRUE)
plot(NULL, 
		type = "n", 
		xlim = range(ages) + c(0, 1), 
		ylim = c(0, ymax), 
		xlab = "", 
		ylab = "", 
		axes = FALSE)
matplot(ages, 
		t(XS),
		add = TRUE,
		type='l',
		lty = 1,
		col = rev(agecols(ncol(XS))),
		lwd = seq(1.5, 2.5, length = ncol(XS)))
# manual axes to match
segments(ages[1], 0, ages[1], ymax)
segments(ages[1], 0, max(ages) + 1, 0)
segments(ages[1],pretty(XS), ages[1] - .4, pretty(XS))
segments(ages[ages %% 5 == 0], 0, ages[ages %% 5 == 0], -.02)
text(ages[1] - .4,pretty(XS),pretty(XS),pos=2,cex=.8)
text(ages[ages %% 5 == 0],-.02,ages[ages %% 5 == 0],pos=1,cex=.8)

# label y locations
yvals <- ncol(XS) - col(XS) + 1 == row(XS)
TTDs  <- c(12:1,"TTD = 0")
keep  <- c(TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE)
text((83:96)[keep]-5,XS[yvals][keep]-.01,TTDs[keep],pos=c(4,4,4,4,1),cex=.8)

text(80,-.058,"Age", cex = 1.2)
text(68.25, .25, "Prevalence", cex = 1.2,srt=90)
dev.off()
# trick to get x dims to match.


# make the age over TTD lines
ttdcols <- colorRampPalette(brewer.pal(9,"Greens")[-c(1:3)])
#dev.new(width=((7 - .7) / diff(c(70,max(ages)+1))) * 13 + .7,height=5)
XSS <- XS[, as.integer(colnames(XS))%%5 == 0]
pdf("Figures/srhpoor_f_TTD_b.pdf",
		width=((7 - .7) / diff(c(70,max(ages)+1))) * 13 + .7,
		height=5)
par(mai=c(.6,.6,.1,.1),xaxs="i",yaxs="i",xpd=TRUE)
plot(NULL, type = "n", xlim = c(0,13),ylim=c(0,.55),xlab="",ylab="",axes=FALSE)
matplot(as.integer(rownames(XS)),XS, 		
		type = 'l',
		lty = 1,
		add=TRUE,
		col = gray(.8),
		lwd = seq(1))
ages <- as.integer(colnames(XS))
keep <- ages %% 5 == 0
matplot(as.integer(rownames(XS)),XS[,keep], 		
		type = 'l',
		lty = 1,
		add=TRUE,
		col = ttdcols(sum(keep)),
		lwd = seq(3,2,length=sum(keep)))
segments(0,0,0,.45);segments(0,0,13,0)
segments(0,seq(0,.8,by=.2),-.4,seq(0,.8,by=.2))
segments(seq(0,13,by=5),0,seq(0,13,by=5),-.015)
text(-.4,seq(0,.8,by=.2),seq(0,.8,by=.2),pos=2,cex=.8)
text(seq(0,13,by=5),-.015,seq(0,13,by=5),pos=1,cex=.8)

yvals <- XS[,keep][c(13,26,37,45)]
xvals <- c(12,12,10,5)
text(xvals,yvals,c(75,80,85,"Age = 90"),pos=c(4,1,1,1),cex=.8)
text(6, -.045, "TTD")
text(-2.6,.4,"Prevalence",srt=90)
dev.off()

pdf(
		"Figures/blankholder.pdf", 
		width = ((7 - .7) / diff(c(70, 101))) * 13 + .7,
		height = 4)
plot(
		NULL, 
		type = "n", 
		axes = FALSE, 
		xlim = c(0, 2), 
		ylim = c(0, 2), 
		xlab = "", 
		ylab = "")
dev.off()
#c(7,((7 - .7) / diff(c(70,101))) * 13 + .7) * 1.2 / 10

#ticks <- seq(0,.8,by=.1)
#graphics.off()
#dev.new(width=7,height=4.2)
#pdf("Figures/SRHpoor_Females.pdf",width=7,height=4.2)
#SurfMap(Dat[["IADL1"]]$Females$Surf[,,3], ticks = ticks, bg=TRUE,ylab="time-to-death",xlab="age",
#		mai=c(.6,.5,.4,1.2))
#text(103,16,"prevalence")
#dev.off()
#pdf("Figures/IADL1_Males.pdf",width=7,height=4.2)
#SurfMap(Dat[["IADL1"]]$Males$Surf[,,3], ticks = ticks, bg=TRUE,ylab="time-to-death",xlab="age",
#		mai=c(.6,.5,.4,1.2))
#text(103,16,"prevalence")
#dev.off()
#
#names(Dat)
