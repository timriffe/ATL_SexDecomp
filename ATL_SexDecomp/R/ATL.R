
# Author: tim
###############################################################################


# make new ATL surface plot.

Dat <- local(get(load("/home/tim/git/ATL_SexDecomp/ATL_SexDecomp/Data/SmArrays.Rdata")))
names(Dat)
library(RColorBrewer)
ramp <-
names()
X <- Dat[["SRHpoor"]]$Females$Surf[,,3]


