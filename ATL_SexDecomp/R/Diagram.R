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

text(col(Birthdays) - .5 + 2000, row(Birthdays) - 1, Birthdays, col = "blue")
#text(2000:2009 + .5, 3, 0, col = "blue")

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


