#SurfA <- SurfaceList[[1]]$Male
#Surf <-SurfA[,,4]

# Surf <- A
SurfMap <- function (Surf, 
		thano = as.integer(rownames(Surf)), 
		chrono = as.integer(colnames(Surf)), 
		colramp = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "OrRd")), space = "Lab"), 
		napprox = 10, 
		xlab = "Years lived, a", 
		ylab = "Years left, y",
		contour=TRUE,
		ticks = NULL,
		legnd = TRUE,
		outline = TRUE,
		mai = c(.5,.5,.5,1.5),
		bg=FALSE) 
{
	if (is.null(ticks)){
		ticks <- pretty(Surf,n=napprox)
	}
	zlim  <- range(ticks)
	n     <- length(ticks)-1
	col   <- rev(colramp(n))
	
	SurfCol <- as.character(cut(Surf,breaks=ticks,labels=col))
	dim(SurfCol) <- dim(Surf)
	
	x <- col(Surf) - 1 + min(chrono)
	y <- row(Surf) - 1 + min(thano)
#	x <- chrono
#	y <- thano
	par(xaxs="i",yaxs="i",xpd=TRUE,mai=mai)
	plot(NULL, type = "n", xlim = c(70,101),ylim=c(0,16),xlab="",ylab="",axes=FALSE,asp=1)
	
	if (bg){
		# draw background in light gray
		rect(70,0,101,15,border=NA,col=gray(.9))
		segments(70,seq(0,15,by=5),101,seq(0,15,by=5),col="white")
		segments(seq(70,100,by=5),0,seq(70,100,by=5),15,col="white")
		segments(c(0,0,0,5,10,15,20,25,30)+70,c(5,10,15,15,15,15,15,15,15),
				c(5,10,15,20,25,30,31,31,31)+70,c(0,0,0,0,0,0,4,9,14),col="white")
	}
	
	# draw cells
	rect(x,y,x+1,y+1,border=NA,col = SurfCol)
	
	
	# outline area
	if (outline){
		ages <- as.integer(colnames(Surf))
		thana <- as.integer(rownames(Surf))
		MinL <- ages[colSums(Surf,na.rm=TRUE) > 0][1]
		MaxL <- rev(ages[colSums(Surf,na.rm=TRUE) > 0])[1]
		Th <- rev(thana[rowSums(Surf,na.rm=TRUE) > 0])[1]
		segments(MinL,0,MinL,Th+1)
		segments(MinL,0,MaxL+1,0)
		segments(MinL,Th+1,MaxL - Th + 1,Th+1)
		segments((MaxL - Th + 1):MaxL,Th:1,(MaxL - Th + 2):(MaxL+1),Th:1)
		segments((MaxL - Th + 1):(MaxL+1),(Th + 1):1,(MaxL-Th+1):(MaxL+1),Th:0)
	}
	
	# axes
	segments(70,seq(0,15,by=5),69.6,seq(0,15,by=5))
	segments(seq(70,100,by=5),0,seq(70,100,by=5),-.4)
	text(69.6,seq(0,15,by=5),seq(0,15,by=5),pos=2,cex=.8)
	text(seq(70,100,by=5),-.4,seq(70,100,by=5),pos=1,cex=.8)
	
	if (legnd){
		# legend
		ticksat <- (ticks - min(ticks))/diff(range(ticks)) * 12 + 2.5
		rect(102,ticksat[-length(ticks)],103.5,ticksat[-1],col=col, border = gray(.6))
		text(103.5, ticksat,ticks,cex=.8,pos=4)
		
		# labels
		text(85,-4,xlab)
		text(68,17,ylab,pos=4)
	}
	if (contour){
		contour(chrono+.5, thano+.5, t(Surf), add = TRUE, col = gray(.2), levels = ticks, labcex = .8)
	}
}


#SurfMap(SurfaceList[[2]]$Male[,,"1915"], bg = TRUE, contour = FALSE)

compareSexes <- function(varname,cohort="1915",SurfaceList){
	Fem <- SurfaceList[[varname]]$Female[,,cohort]
	Mal <- SurfaceList[[varname]]$Male[,,cohort]
	ticks <- pretty(c(Fem,Mal),n=napprox)
	par(mfrow = c(1, 2))
	SurfMap(Mal,contour = FALSE,ticks=ticks)
	text(80,19,paste("Males",varname,cohort,"Cohort" ))
	SurfMap(Fem,contour = FALSE,ticks=ticks)
	text(80,19,paste("Females",varname,cohort,"Cohort" ))
}
