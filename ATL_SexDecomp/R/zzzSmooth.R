
if (system("hostname",intern=TRUE) %in% c("triffe-N80Vm","tim-ThinkPad-L440")){
	# if I'm on the laptop
	setwd("/home/tim/git/ATL_SexDecomp/ATL_SexDecomp")
} else {
	# in that case I'm on Berkeley system, and other people in the dept can run this too
	setwd(paste0("/data/commons/",system("whoami",intern=TRUE),"/git/ATL_SexDecomp/ATL_SexDecomp"))
}
source("R/zzzLoess2Deprecate.R")
source("R/SurfMap.R")
Dat 	<- local(get(load("Data/Data_long_imputed.Rdata")))
Dat 	<- as.data.frame(Dat)
Dat     <- Dat[!is.na(Dat$b_yr), ]
Dat     <- Dat[Dat$age >= 65, ]
Dat <- Dat[Dat$b_yr >= ]
unique(Dat$sex)
"adl5_"
coh5 <- 1915
table(Dat$b_yr[Dat$sex == "f"])
table(Dat$b_yr[Dat$sex == "m"])
plot(table(Dat$b_yr[Dat$sex == "f"]))
colnames(Dat)


getcoh <- function(varname = "iadl5_",coh5 = 1915,sex = "f", Dat){
	FitLoess(varname = varname, 
			Dat[Dat$b_yr >= (coh5 - 2) & Dat$b_yr <= (coh5 + 2), ], 
			sex = sex,
			t.age = 0:12,    # some 5-year cohorts simply don't have 15 years, cut it lower
			c.age = 70:100,  # standard matrix size, though we may NA certain unobserved cells
			span = .7, # will vary
			.Coh5 = coh5)
}

cohs <- c(1910,1915,1920)
sexes <- c("f","m")
vars <- c("iadl5_","adl5_")

Results <- list()
for (sex in sexes){
	Cohorts <- list()
	for (coh in cohs){
		varsi <- list()
		for (varname in vars){
			varsi[[varname]] <- getcoh(
					varname = varname,
					coh5 = coh,
					sex = sex, 
					Dat)$Surf[,,1]
		}
		Cohorts[[as.character(coh)]] <- varsi
	}
	Results[[sex]] <- Cohorts
}


names(Results$f[[1]])

ticks <- seq(0,1,by=.1)

Results <- local(get(load("Data/ResultsIADL_ADL.Rdata")))

# IADL5
graphics.off()
# females
dev.new(height=10,width=5)
par(mfrow=c(3,1))
SurfMap(Results[["f"]][["1910"]][[1]], ticks = ticks)
SurfMap(Results[["f"]][["1915"]][[1]], ticks = ticks)
SurfMap(Results[["f"]][["1920"]][[1]], ticks = ticks)
# males
dev.new(height=10,width=5)
par(mfrow=c(3,1))
SurfMap(Results[["m"]][["1910"]][[1]], ticks = ticks)
SurfMap(Results[["m"]][["1915"]][[1]], ticks = ticks)
SurfMap(Results[["m"]][["1920"]][[1]], ticks = ticks)


# ADL5
graphics.off()
# females
dev.new(height=10,width=5)
par(mfrow=c(3,1))
SurfMap(Results[["f"]][["1910"]][[2]], ticks = ticks)
SurfMap(Results[["f"]][["1915"]][[2]], ticks = ticks)
SurfMap(Results[["f"]][["1920"]][[2]], ticks = ticks)
# males
dev.new(height=10,width=5)
par(mfrow=c(3,1))
SurfMap(Results[["m"]][["1910"]][[2]], ticks = ticks)
SurfMap(Results[["m"]][["1915"]][[2]], ticks = ticks)
SurfMap(Results[["m"]][["1920"]][[2]], ticks = ticks)

save(Results,file="Data/ResultsIADL_ADL.Rdata")


