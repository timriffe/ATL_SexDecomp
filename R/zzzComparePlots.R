# Author: tim
###############################################################################

if (system("hostname",intern=TRUE) %in% c("triffe-N80Vm","tim-ThinkPad-L440")){
	# if I'm on the laptop
	setwd("/home/tim/git/ATL_SexDecomp/ATL_SexDecomp")
} else {
	# in that case I'm on Berkeley system, and other people in the dept can run this too
	setwd(paste0("/data/commons/",system("whoami",intern=TRUE),"/git/ATL_SexDecomp/ATL_SexDecomp"))
}

SurfaceList <- local(get(load("Data/SurfaceList.Rdata")))
vars <- names(SurfaceList)
varname <- vars[1]
pdf("Figures/SurfCompare1915.pdf",width=10,height=4)
for (varname in vars){
	compareSexes(varname, "1915",SurfaceList)
}
dev.off()

