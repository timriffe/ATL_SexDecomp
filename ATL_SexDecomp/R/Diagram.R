# Author: tim
###############################################################################
library(reshape2)
ECtri2census <- function(Deaths){
	Deaths              <- Deaths[Deaths$Cohort != Deaths$Period,]
	PC 					<- acast(Deaths, 
			Period~Cohort, 
			sum, 
			value.var = "Deaths")
	
	CPcs 				<- apply(PC, 1, cumsum)
	
	CPcs[t(PC == 0)] 	<- NA
	Census 				<- melt(CPcs,
			varnames = c("Cohort","Period"), 
			value.name = "Population")
	Census 				<- Census[!is.na(Census$Population), ]
	Census$Age 			<- Census$Period - Census$Cohort - 1 
	
	Census
}

# underlying mortality:
stationary <- matrix(c(20,20,20,20,20,20,20,20,20,20,
				60,60,60,60,60,60,60,60,60,60,
				35,35,35,35,35,35,35,35,35,35), ncol = 10, nrow = 3, byrow = TRUE,
		dimnames = list(0:2,2000:2009))
change <- matrix(c(0, 0, 0, 0,    -10,  -10,-15,-15,-15,-15,
				   0,0,0,0,       5,   5, 10,10,10,10,
				   0,0,0,0,       0, 5, 5, 5, 5, 5 ), ncol = 10, nrow = 3, byrow = TRUE,
		dimnames = list(0:2,2000:2009))


upper <- stationary + change
lower <- upper

upper[1, ] <- upper[1, ] + 5
upper[3, ] <- upper[3, ] - 5
lower[1,6] <- lower[1,6] - 5

upper[3, ] <- upper[3, ] - 25
lower[3, ] <- lower[3, ] - 25
#lower[2, 1:4] <- lower[2,1:4] + 5
#upper[2, 6:10] <- upper[2,6:10] + 5


Upper 			<- melt(upper, varnames = c("Age", "Period"), value.name = "Deaths")
Lower 			<- melt(lower, varnames = c("Age", "Period"), value.name = "Deaths")
Upper$Cohort 	<- Upper$Period - Upper$Age - 1
Lower$Cohort 	<- Lower$Period - Lower$Age
Deaths 			<- rbind(Upper, Lower)

# get Census counts!

Census <- acast(ECtri2census(Deaths), Age~Period, value.var = "Population")
Census <- cbind(Census,Census[,ncol(Census)])
colnames(Census)[11] <- 2010
#plot(2000:2010,colSums(Census))
# now get birthdays!
Birthdays <- Census[,-1] + lower
Birthdays <- rbind(Birthdays,rep(0,10))
rownames(Birthdays)[4] <- 3

TriCol <- gray(.7)

#draw.lower.tri <- function()


LTtri <- function(Census, Birthdays){
	
	pl <- Census[,-1] / Birthdays[-4,]
	pu <- Birthdays[-1,] / Census[,-ncol(Census)]
	
	ax <- 1/3
	pinter <- rbind(pl[1, ], pu[1, ], pl[2, ], pu[2, ], pl[3, ], pu[3, ])
	
	linter <- apply(rbind(1, pinter), 2, cumprod)
	#matplot(linter, type = 'l')
	Linter <- (linter[-1, ] + linter[-nrow(linter), ]) / 4
	dinter <- -apply(linter,2,diff)
	minter <- dinter / Linter
	#matplot(minter, log='y',type = 'l')
	Tinter <- apply(Linter,2,function(x){
				rev(cumsum(rev(x)))
			})
	einter <- Tinter / linter[-nrow(linter), ]
}


#N <- 1e6
#x <- runif(N)
## birthdays.
#surv <- runif(N)
#
#keep <- (x + surv) < 1
#x <- x[keep]
#surv <- surv[keep]
#mean(surv)

# 1) calculate exposures for triangles.

# 2) calculate AP and AC

# 3) figure out how to apply g(y) to the death counts: by triangle?

# ya, for triangles, assuming uniformity we have 1/3, 2/3, 1 1/3, 1 2/3, etc.. applied to the triangles.
# so this will give proportions in each vertical segment. Good.
graphics.off()
y <- c(1/3,2/3,1+1/3,1+2/3,2+1/3,2+2/3)
gy <- c(.9,.6,.2,.1,.05,.02)
plot(y,gy)
# then for periods we can average G(a) on either left or right?
		
# get proportions of G for each census slice:

Census

lower

upper
lowermult <- c(1/3,1+1/3,2+1/3)
uppermult <- c(2/3,1+2/3,2+2/3)

lowermult <- c(.9,.2,.05)
uppermult <- c(.6,.1,.02)
G0APl <- lower * lowermult 
G0APu <- upper * uppermult
G0APl <- melt(G0APl,varnames=c("Age","Period"),value.name = "G0")
G0APu <- melt(G0APu,varnames=c("Age","Period"),value.name = "G0")

G0APl$Cohort <- G0APl$Period - G0APl$Age 
G0APu$Cohort <- G0APu$Period - G0APu$Age -1

G0APC <- rbind(G0APl, G0APu)
G0APC <- G0APC[G0APC$Cohort >= 2000, ]
G0 <- tapply(G0APC$G0, G0APC$Cohort, sum)
G0[c("2007","2008","2009")] <- G0["2006"]

# now G1: this is sooo manual!!

lowermult <- c(.9,.2)
uppermult <- c(.6,.1)
G1APl <- lower[-1,] * lowermult 
G1APu <- upper[-1,] * uppermult
G1APl <- melt(G1APl,varnames=c("Age","Period"),value.name = "G1")
G1APu <- melt(G1APu,varnames=c("Age","Period"),value.name = "G1")

G1APl$Cohort <- G1APl$Period - G1APl$Age 
G1APu$Cohort <- G1APu$Period - G1APu$Age -1

G1APC <- rbind(G1APl, G1APu)
G1APC <- G1APC[G1APC$Cohort >= 1999, ]
G1 <- tapply(G1APC$G1, G1APC$Cohort, sum)
G1[c("2007","2008")] <- G1["2006"]
G1 <- G1[1:10]

lowermult <- c(.9)
uppermult <- c(.6)
G2APl <- lower[-c(1,2),] * lowermult 
G2APu <- upper[-c(1,2),] * uppermult
G2APl <- matrix(G2APl,nrow=1,byrow=TRUE,dimnames=list(2,names(G2APl)))
G2APu <- matrix(G2APu,nrow=1,byrow=TRUE,dimnames=list(2,names(G2APu)))
G2APl <- melt(G2APl,varnames=c("Age","Period"),value.name = "G2")
G2APu <- melt(G2APu,varnames=c("Age","Period"),value.name = "G2")

G2APl$Cohort <- G2APl$Period - G2APl$Age 
G2APu$Cohort <- G2APu$Period - G2APu$Age -1

G2APC <- rbind(G2APl, G2APu)
G2APC <- G2APC[G2APC$Cohort >= 1998, ]
G2 <- tapply(G2APC$G2, G2APC$Cohort, sum)
G2["2007"] <- G2["2006"]
G2 <- G2[1:10]

Gbirthdays <- matrix(c(rbind(G0,G1,G2)),nrow=3,dimnames=list(0:2,2000:2009))

############################################

print(xtable::xtable(cbind(y,gy)),row.names=FALSE)

#\begin{table}[ht]
#\centering
#\begin{tabular}{rr}
#\hline
#& TTD & $g(y)$ \\ 
#\hline
#$\frac{1}{3}$ & 0.90 \\ 
#$\frac{2}{3}$ & 0.60 \\ 
#$1\frac{1}{3}$ & 0.20 \\ 
#$1\frac{2}{3}$ & 0.10 \\ 
#$2\frac{1}{3}$ & 0.05 \\ 
#$2\frac{2}{3}$ & 0.02 \\ 
#\hline
#\end{tabular}
#\end{table}
##############################################
# now we need exposures AP:

# exposures:
(pop1 + pop2) / 2 + (dl - du) / 6

Exposure <- (Census[,1:10] + Census[,-1]) /2 + (lower - upper) / 6
Mx <- (lower + upper) / Exposure
ax <- (lower * 1/3 + upper * 2/3) / (upper + lower)
qx <-  Mx / (1 + (1 - ax) * Mx)
qx[qx>1] <-1

px <- 1-qx
lx <- rbind(1,apply(px,2,cumprod))
dx <- -apply(lx,2,diff)
Lx <- lx[-4,] - ax * dx 
e0 <- colSums(Lx)


Gab <- Gbirthdays / Birthdays[-4,]
GAP <- (Gab + rbind(Gab[-1, ], 0)) / 2
eUAP <- colSums(GAP * Lx)

plot(2000:2009, eUAP/e0, type = 'l')
plot(2000:2009, e0 - eUAP, type = 'l')

period <- data.frame(Year = 2000:2009, e0 = e0, eH = e0 - eUAP)

# now cohort EH

#mx 				<- Dx / Exp																						# Eq 50  MPv5
#qx 				<- Dx / (pop + dl)																		# Eq 70  MPv5 (same in V6)
## ax 				<- ((1 / 3) * dl + (2 / 3) * du) / Dx  						  # Eq 71  MPv5
## prefer solve for ax from identity formula- more flexible given version
#ax        <- ((mx + 1) * qx - mx) / (mx * qx)    

colnames(Birthdays) <- 2000:2009
BD <- melt(Birthdays, varnames = c("Age","Period"),value.name = "Birthdays")
BD$Cohort <- BD$Period - BD$Age
BD <- BD[BD$Cohort %in% c(2000:2006),]
BAC <- acast(BD, Age~Cohort, value.var = "Birthdays")

pxc <- BAC[-1,] / BAC[-4,] 

lxc <- apply(rbind(1,pxc),2,cumprod)
Lxc <- (lxc[-1, ] + lxc[-4, ]) / 2

e0c <- colSums(Lxc)

# now get Gac
GD <- melt(Gab, varnames = c("Age","Period"),value.name = "Gb")
GD$Cohort <- GD$Period - GD$Age
GD <- GD[GD$Cohort %in% c(2000:2006),]
Gac <- acast(GD, Age~Cohort, value.var = "Gb")
GAC <- (Gac + rbind(Gac[-1, ], 0)) / 2

eUAC <- colSums(GAC * Lxc)

cohort <- data.frame(Year = 2000:2006, e0 = e0c, eH = e0c - eUAC)

plot(period$Year, period$eH, type = 'l')
lines(cohort$Year, cohort$eH, col = "blue")
# try the odd LT method: exiting divided by entering square:

pdf("Figures/DiagramLexis.pdf",height=5,width=10)
par(xaxs="i",yaxs = "i", xpd=FALSE)
plot(NULL, type = "n", asp = 1, xlim = c(2000, 2010), ylim = c(0, 3), axes = FALSE, xlab = "", ylab = "")
segments(1995:2010, 0, 1998:2013, 3, col = TriCol)
segments(2000:2010, 0, 2000:2010, 3, col = TriCol)
segments(2000, 0:3, 2010, 0:3, col = TriCol)

# years
text(2000:2010,-.3,2000:2010,srt=90,xpd=TRUE)
# ages
text(2000, 0:3, 0:3, pos = 2, xpd=TRUE)
text(col(lower) - 1 / 3 + 2000, row(lower) + 1 / 3 - 1, lower, col = "red")
text(col(upper) - 2 / 3 + 2000, row(upper) + 2 / 3 - 1, upper, col = "red")

text(col(Census) - 1 + 2000, row(Census) - .5, Census, xpd = TRUE)

text(col(Birthdays) - .8 + 2000, row(Birthdays) - 1, Birthdays, col = "blue")
text(col(Gbirthdays) - .2 + 2000, row(Gbirthdays) - 1, round(Gbirthdays), col = "forestgreen")

text(2003,-.8,"Year",cex=2,xpd=TRUE)
text(1999.5,3.5,"Age",cex=2,xpd=TRUE)

#legend
text(2006,-.7,"Census Pop",pos=4,xpd=TRUE)
text(2006,-.9,"Birthdays",pos=4,xpd=TRUE,col="blue")
text(2006,-1.1,"Unhealthy",pos=4,xpd=TRUE,col="forestgreen")
text(2006,-1.3,"Deaths (tri)",pos=4,xpd=TRUE,col="red")
dev.off()

pdf("Figures/TTDgy.pdf",height=5,width=10)
plot(c(0,y), c(1,gy), type = 'l', xlab = "", ylab = "", axes = FALSE, 
		ylim=c(0,1),xlim=c(0,3), col = "forestgreen", lwd = 2)
points(y,gy,pch=19,col="forestgreen")
segments(0,0,3,0)
segments(0,0,0,1)
segments(seq(0,3,by=1),0,seq(0,3,by=1),-.04,xpd=TRUE)
segments(0,seq(0,1,by=.5),-.06,seq(0,1,by=.5),xpd=TRUE)
text(seq(0,3,by=1),-.04,seq(0,3,by=1),pos=1,xpd=TRUE)
text(-.06,seq(0,1,by=.5),seq(0,1,by=.5),pos=2,xpd=TRUE)
text(-.2,1.2,"g(y)",xpd=TRUE,cex=2)
text(1.5,-.2,"time to death",xpd=TRUE,cex=2)
dev.off()

pdf("Figures/e0eHtoy.pdf",height=5,width=8)
par(mai=c(.5,.5,.5,.5))
plot(2000:2009, period$e0, type = 'l', ylim = c(0,1.6),lwd=2,axes=FALSE, xlab = "",ylab = "")
lines(2000:2006, cohort$e0, col = "blue",lwd=2)
#polygon(c(2000:2009,2009:2000),c(period$e0,rev(period$eH)), col = "#00000020", border= FALSE)
#polygon(c(2000:2006,2006:2000),c(cohort$e0,rev(cohort$eH)), col = "#0000FF20", border= FALSE)
lines(2000:2009, period$eH,lwd=2, lty = 2)
lines(2000:2006, cohort$eH, col = "blue", lwd=2, lty = 2)

segments(2000,0,2009,0)
segments(2000,0,2000,1.6)
segments(seq(2000,2009,by=1),0,seq(2000,2009,by=1),-.04,xpd=TRUE)
segments(2000,seq(0,1.6,by=.5),2000-.15,seq(0,1.6,by=.5),xpd=TRUE)
text(seq(2000,2009,by=1),-.04,seq(2000,2009,by=1),pos=1,xpd=TRUE)
text(2000-.15,seq(0,1.6,by=.5),seq(0,1.6,by=.5),pos=2,xpd=TRUE)

text(2001.5,1.5,"cohort e(0)",col="blue",cex=1.2)
text(2004,1.3,"period e(0)",cex=1.2)

text(2001.5,.85,"cohort\nhealthy e(0)",col="blue",cex=1.2,pos=4)
text(2004,.5,"period healthy e(0)",cex=1.2)
dev.off()

matplot(0:2,GAP, type = 'l', col = "#00000050", lty = 1)


