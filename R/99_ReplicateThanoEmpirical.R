
# Author: tim
###############################################################################

# this is a side endeavor. Since the smoothing procedure changed between the present
# paper are earlier version from
# Riffe et al (2017) "Time-to-death patterns in markers of age and disability" VYPR v14
# I'd like to remake the  box plot distributions in the final figure of that paper to see if same

Results               <- local(get(load("Data/resultsP.Rdata")))
varnames              <- sapply(Results, "[[", "var")
names(Results)        <- varnames

library(reshape2)
ResultsLong           <- melt(Results)
colnames(ResultsLong) <- c("ta","ca","Cohort","pi","Delete","Sex","var")
ResultsLong$Delete    <- NULL
ResultsLong$pi        <- suppressWarnings(as.numeric(ResultsLong$pi))
ResultsLong           <- ResultsLong[!is.na(ResultsLong$pi), ]

# cut down to same time reference
ResultsLong           <- ResultsLong[ResultsLong$Cohort >= 1905 & ResultsLong$Cohort < 1930, ]
ResultsLong$Cohort5   <- ResultsLong$Cohort - ResultsLong$Cohort %% 5  
# take simple means within ca, ta, Cohort, var, Sex
ResultsLong   		  <- data.table(ResultsLong)

# choosing median here, could just as well be mean.
ResultsLong5  		  <- ResultsLong[, list(pi = median(pi, na.rm=TRUE)), 
		                             by = list(var, Sex, ca, ta, Cohort5)]
ResultsLong5  		  <- as.data.frame(ResultsLong5)
colnames(ResultsLong5)[5] <- "Cohort"
ResultsLong5L 		  <- split(ResultsLong5,with(ResultsLong5, list(var, Sex, Cohort)))

# do a quick comparison with ThanoEmpirical results:

# function used, given a thano x chrono surface:
get_r <- function(A){
	c(T  = abs(cor(A$pi, A$ta, use = "complete.obs")), 
			A  = abs(cor(A$pi, A$ca, use = "complete.obs")),
			L  = abs(cor(A$pi, A$ta + A$ca, use = "complete.obs")),
			M  = abs(cor(A$pi, A$ca - A$ta, use = "complete.obs"))
	)
}
	
# --------------------------------------------
# This gets the correlations for each dim, cohort, sex, and span.
	Results_r <- do.call(rbind,
			lapply(ResultsLong5L, function(X){
						out     <- X[1:4, ]
						out$Dim <- c("T","A","L","M")
						out$r   <- get_r(X)
						out$pi  <- NULL
						out$ca  <- NULL
						out$ta  <- NULL
						out
					})
	)
	
ind <- Results_r$r[Results_r$Dim=="T"] - Results_r$r[Results_r$Dim=="A"] > .8
# Male cesddepr 1915, Female srhpoor 1920

#Results_r[Results_r$Dim=="T",][ind, ]
	
	comparison.tim <- FALSE
	if (comparison.tim){
# compare these with the original results
		Results_loess_r   	<- local(get(load("/home/tim/git/ThanoEmpirical/ThanoEmpirical/Data/Correlations.Rdata")))
	
# match ordering
		Results_r       	<- Results_r[with(Results_r,order(var, Sex, Cohort, Dim)), ]
		Results_loess_r   	<- Results_loess_r[with(Results_loess_r,order(var, sex, Cohort, Dim)), ]
		Results_r$var   	<- gsub("_", "",Results_r$var)
		
		Results_loess_r.5 	<- Results_loess_r[Results_loess_r$span == "0.5", ]
		Results_loess_r.7 	<- Results_loess_r[Results_loess_r$span == "0.7", ]
		Results_loess_r.9 	<- Results_loess_r[Results_loess_r$span == "0.9", ]
	}
# match cohorts
	Results_r       	<- Results_r[Results_r$Cohort %in% unique(Results_loess_r.5$Cohort), ]
	
# TR: hmmm. positive correlation, but that cloud is fatter than
# I'd like!
#plot(Results_r$r, Results_loess_r.5$r)
#plot(Results_r$r[Results_r$Cohort == 1915], Results_loess_r.5$r[Results_loess_r.5$Cohort == 1915])
	
# so let's start be redoing Figure 5a, 5b too see if they change much
	
	Hist7      <- Results_r[Results_r$Cohort == 1915, ]
	Hist7$rbin <- round(Hist7$r * 100) %/% 10 * 10
	Hist7$rbin <- as.factor(Hist7$rbin)
	Hist7$R    <- Hist7$r * 100
	Hist7$Dim  <- as.factor(Hist7$Dim )

	if (make.figs){
# Figures that ought to be mostly comparable with
# Figure 5 in ThanoEmpirical paper. Hopefully
# these histograms have a similar shape!!
		library(lattice)
		pdf("Figures/TE_Figure5a_v3.pdf", width = 3, height = 7)
		histogram(~R | Dim, 
				data = Hist7[Hist7$Sex == "Female", ], 
				col = gray(.4), 
				par.settings = list(strip.background = list(col = gray(.9))), 
				layout = c(1,4), 
				type = "count",
				breaks = seq(0,100, by = 10),
				index.cond = list(c(3, 4, 1, 2)), 
				xlab = list(label = "correlation coef * 100"),
				ylab = list(label = "variable count"),
				ylim = c(0, 47))
		dev.off()
		
		pdf("Figures/TE_Figure5b_v3.pdf", width = 3, height = 7)
		histogram(~R | Dim, 
				data = Hist7[Hist7$Sex == "Male", ], 
				col = gray(.4), 
				par.settings = list(strip.background = list(col = gray(.9))), 
				layout = c(1, 4), 
				type = "count",
				breaks = seq(0, 100, by = 10),
				index.cond = list(c(3, 4, 1, 2)), 
				xlab = list(label = "correlation coef * 100"),
				ylab = list(label = "variable count"),
				ylim = c(0, 47))
		dev.off()
		
	} # end fig chunk
	
	
# TR will want a copy of this output.
	#save(ResultsLong, file = "Data/ResultsLongBoot.Rdata")
	save(Results_r, file = "Data/CorrelationResultsBoot.Rdata")
	
	getwd()


