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
		ticks,
		legnd = TRUE,
		outline = TRUE,
		mai = c(.5, .5, .5, 1.5),
		bg = FALSE,
		xlim = range(chrono) + c(0,1),
		ylim = range(thano) + c(0,1)) 
{
	if (missing(ticks)){
		ticks <- pretty(Surf, n = napprox)
	}
	zlim         <- range(ticks)
	n            <- length(ticks) - 1
	col          <- rev(colramp(n))
	
	SurfCol      <- as.character(cut(Surf, breaks = ticks, labels = col))
	dim(SurfCol) <- dim(Surf)
	
	x            <- col(Surf) - 1 + min(chrono)
	y            <- row(Surf) - 1 + min(thano)

	par(xaxs = "i", yaxs = "i", xpd = TRUE, mai = mai)
	plot(NULL, type = "n", xlim = xlim, ylim = ylim, xlab = "", ylab = "", axes=FALSE,asp=1)
	
	xlines <- chrono[chrono %% 5 == 0]
	ylines <- thano[thano %% 5 == 0]
	if (bg){
		# draw background in light gray
		
		left   <- min(chrono)
		right  <- max(chrono) + 1
		bottom <- min(thano)
		top    <- max(thano) + 1
		rect(left, bottom, right, top, border = NA, col = gray(.9))
		
		
		segments(left, ylines, right, ylines, col = "white")
		segments(xlines, bottom, xlines, top, col = "white")
		
		# diagonals are always a pain in the ass to figure out.
		# do it in steps:
		
		x2      <- seq(min(xlines), max(xlines) + 5 * length(ylines), 5) + bottom
        x1      <- x2 - (top - bottom)
		
		y2      <- rep(bottom, length(x2))
		y1      <- rep(top, length(x1))
		
		# trim 1
		diff1               <- x1[x1 < left] - left
		y1[1:length(diff1)] <- y1[1:length(diff1)] + diff1
		x1[x1 < left]       <- left
		
		# trim 2
		diff2                                           <- x2[x2 > right] - right
		y2[(length(y2) - length(diff2) + 1):length(y2)] <- diff2
		x2[x2 > right]                                  <- right
		segments(x1, y1, x2, y2, col = "white")
	}
	
	# draw cells
	rect(x, y, x + 1, y + 1, border = NA, col = SurfCol)
	
	
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
	segments(left, ylines, left - .4, ylines)
	segments(xlines, bottom, xlines, bottom - .4)
	text(left - .4, ylines, ylines, pos = 2, cex = .8)
	text(xlines, bottom - .4, xlines, pos = 1, cex = .8)
	
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

