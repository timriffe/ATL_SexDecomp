
# Getting period g'(a) from a fixed 1915 cohort g(y). 
# Calculating period HLE in 1995 and 2010 assuming stationary period mortality
# Demonstrating that decompositions will show a role for morbidity
# When all that is changing is the mortality scheduled implied by g'(a,y)

#-----------------------------------------------------------------
# First, determining the g(y) from the 1915 birth cohorts
#-----------------------------------------------------------------

# data

ADL3 <- local(get(load("C:\\Dropbox\\sex decomp Tim\\ADL3.Rdata"))) # modify as necessary
# can't have negatives
imp0 <- function(x){
  x[x < 0]<- 0
  x
}
Males   <- imp0(ADL3$Male$Surf)
Females <- imp0(ADL3$Female$Surf)
# these are arrays. the third level is the 1915-1919 cohort.
# that's the one we'll work with.

# quick peek to figure out how to proceed.
# image(t(Females[,,3]))

# note that we're missing a piece of the upper triangle (bottom of matrix)
# too few data points there to get anything out of it.
# let's just impute it with something reasonable-seeming
# so we can move forward.

#image(t(Males[,,3]))
#a <- 70:100
#plot(a,Males["12",,3],type='l',col="green",ylim=c(0,.06))
#lines(a,Males["11",,3],col="royalblue")
#lines(a,Males["10",,3],col="blue")   # (looks similar for females)

# let's just apply the same proportional change that we see btwn
# the last two rows to each extrapolated row. If it drops
# below zero, we impute zero. Easy, cheap.
extendOmatic <- function(M){
  # some columns may have all NAs, be sure to keep it that way
  # assume data are in this very particular shape. OK, because
  # this is a once-off hack (we hope)
  start <- M[nrow(M), ]
  ratio <- start /  M[nrow(M) - 1, ]
  
  Col <- ncol(M)
  Row <- nrow(M)
  
  # the part we need to fill in with a triangle of extrapolated data
  NewBlock <- matrix(nrow = Col - Row, ncol = Col)
  
  for (i in 1:nrow(NewBlock)){
    NewBlock[i, ] <- start * ratio ^ i # cumprod
  }
  
  # append matrix
  Mout <- rbind(M, NewBlock)
  
  # need to respect NA diagonal.
  Nums <- row(Mout) + col(Mout)
  
  # need a marker on the right side in the top row
  Ind  <- max((1:ncol(Mout))[!is.na(Mout[1,])])
  Mout[Nums > Nums[1,Ind]] <- NA
  
  # clean extrap, just in case
  Mout[Mout < 0] <- 0
  Mout[Mout > 1] <- 1
  
  # include rownames
  rownames(Mout) <- 1:nrow(Mout) - 1
  Mout
}
# cool, now we can do stuff!
M1915e <- extendOmatic(Males[,,3])     
F1915e <- extendOmatic(Females[,,3])


# make sure NA placement matches, you never know...
all(colSums(!is.na(M1915e)) == colSums(!is.na(F1915e)))
all(rowSums(!is.na(M1915e)) == rowSums(!is.na(F1915e)))

# slice off unneeded columns (i.e. ages not observed for this cohort)
M1915e <- M1915e[, colSums(!is.na(M1915e)) > 0] 
F1915e <- F1915e[, colSums(!is.na(F1915e)) > 0] 
M1915e <- M1915e[rowSums(!is.na(M1915e)) > 0, ] 
F1915e <- F1915e[rowSums(!is.na(F1915e)) > 0, ] 

# chrono ages + thano ages, should be same lengths
(ca <- as.integer(colnames(M1915e)))
(ta <- as.integer(rownames(M1915e)))

#image(t(M1915e-F1915e))
#
par(mfrow=c(1,2))
image(t(M1915e), zlim = c(0,.6))
image(t(F1915e), zlim = c(0,.6))

############################################################################
# note: the people alive in x are the people that die in x+

# we learn that there is another direction in which we ought to extend
# the data. Since we only have up to age 95, there are people alive in l(x)
# that we don't have in this window of d(x)...

# That last open age group of d(x)
# is therefore present in each preceding age of l(x), and if we want to see the
# chrono consequences of these morbidity patterns then we need to be able to 
# compose an l(x) out of it's various d(x) components, taking into account
# the thanatological trajectories of each d(x). So we need values for those 
# missing d(x), for the ages > 95. I'm inclined to just repeat the experience of those
# that died at age 95, shifting it right. Then the matrix will get bigger.
# another problem will be that we need to extend cMx until the cohort
# is mostly extinct. In that case, how about doing the HMD extendo-O-matic
# thing and borrowing from the preceding cohorts.

# in short, we need to 'close out' both the mortality data and the morbidity data. 
# In the end, this will have little leverage on results. If mortality were
# a lot lower, it would have leverage. Realistically, but the time mortality
# ever does get that low, we'll actually have the data and not need to extrapolate!!

############################################################
# OK, step 1, for fudging a working example:
# extend our morbidity surfaces to include people out to omega.
# let's say omega is 110. Because HMD.

openm       <- rev(diag(M1915e[nrow(M1915e):1,]))
openf       <- rev(diag(F1915e[nrow(F1915e):1,]))
# again, for the sake of closing out. More thought could go into this.
anew        <- 96:110
RightSide   <- matrix(nrow=nrow(M1915e),ncol=length(anew),dimnames=list(rownames(M1915e),anew))
M1915ee     <- cbind(M1915e, RightSide)
F1915ee     <- cbind(F1915e, RightSide)
N           <- ncol(M1915e)
for (i in 1:ncol(RightSide)){
  for (z in 1:length(openm)){
    M1915ee[z,N+1+i-z] <- openm[z]
    F1915ee[z,N+1+i-z] <- openf[z]
  }
}

M1915ee <- rbind(M1915ee,M1915ee[1:15,] * 0)
F1915ee <- rbind(F1915ee,F1915ee[1:15,] * 0)

M1915ee[upper.tri(M1915ee)[nrow(M1915ee):1,]] <- NA
F1915ee[upper.tri(F1915ee)[nrow(F1915ee):1,]] <- NA
rownames(F1915ee) <- rownames(M1915ee) <- 1:nrow(M1915ee) - 1

Mgy <- apply(M1915ee,1,mean,na.rm="T")
Fgy <- apply(F1915ee,1,mean,na.rm="T")


#--------------------------------------------------------------------
#-- next: get period mortality data from HMD
#--------------------------------------------------------------------
USAFLT <- read.table("C:\\hmd_statistics\\lt_female\\fltper_1x1\\USA.fltper_1x1.txt",
                     head=T,skip=2,na.strings=".")
USAMLT <-  read.table("C:\\hmd_statistics\\lt_male\\mltper_1x1\\USA.mltper_1x1.txt",
                      head=T,skip=2,na.strings=".")
USAFLT[,2] <- rep(0:110,length.out=dim(USAFLT)[1])
USAMLT[,2] <- rep(0:110,length.out=dim(USAMLT)[1])

fmxp1 <- subset(USAFLT,Year==1995 & Age>=74 & Age<=110)$mx
fmxp2 <- subset(USAFLT,Year==2010 & Age>=74 & Age<=110)$mx

mmxp1 <- subset(USAMLT,Year==1995 & Age>=74 & Age<=110)$mx
mmxp2 <- subset(USAMLT,Year==2010 & Age>=74 & Age<=110)$mx


#---------------------------------------------------------------
# Calculating the resulting g'(a) from the stationary mortality 
# distributions of 1995 and 2010
#---------------------------------------------------------------

# little functions
  mx2lx <- function(mx){
    N <- length(mx)
    mx <- c(mx[-c(N-1,N)],sum(mx[(N-1):N]))
    c(1, exp(-cumsum(mx)))
  }
  
  lx2dx <- function(lx){
    -diff(c(lx,0))
  }
  
  lx2e0 <- function(lx){
    sum((lx + c(lx[-1],0)) / 2)
  }
  
  lx2ex <- function(lx){
    lx <- lx / lx[1]
    Lx <- lx2Lx(lx)
    sum(Lx)
  }

  lx2Lx <- function(lx){
    (lx + c(lx[-1],0)) / 2
    
  }


  expit <- function(g){
    exp(g) / (1+exp(g))
  }

  qx2lx <- function(qx){
    cumprod(1-c(0,qx))
  }

  mx2qx <- function(mx){
    mx / (1 + .5 * mx)
  }



# needed quantities
flx1 <- mx2lx(fmxp1)
mlx1 <- mx2lx(mmxp1)

flx2 <- mx2lx(fmxp2)
mlx2 <- mx2lx(mmxp2)

fdx1 <- lx2dx(flx1)
fdx2 <- lx2dx(flx2)

mdx1 <- lx2dx(mlx1)
mdx2 <- lx2dx(mlx2)


# Putting gy on a matrix
n <- length(fdx2)
FThx <- MThx <- matrix(0,n,n)
# columns = crono ; rows = lifespans. despues ponderamos con d(x)
for (i in 1:n){
  FThx[i,1:(n+1-i)] <- Fgy[(n+1-i):1]
  MThx[i,1:(n+1-i)] <- Mgy[(n+1-i):1]
 }

#image(FThx)
head(FThx)
# longest lifespans at top, going right over rows.
# that's why we do rev dx
plot(fdx1)
FM1 <- FThx*rev(fdx1)
MM1 <- MThx*rev(mdx1)
FM2 <- FThx*rev(fdx2)
MM2 <- MThx*rev(mdx2)


plot(Fgap1)



sum(FThx*rev(fdx1)) # esperanza de morbididad 1
sum(FThx*rev(fdx2)) # esperanza de morbididad 2

sum(MThx*rev(mdx1)) # esperanza de morbididad 1
sum(MThx*rev(mdx2)) # esperanza de morbididad 2


lx2e0(flx1) - sum(FThx*rev(fdx1)) # esperanza de salud 1
lx2e0(flx2) - sum(FThx*rev(fdx2)) # esperanza de salud 2 # caso 2 es una mejora muy fuerte!

lx2e0(mlx1) - sum(MThx*rev(mdx1)) # esperanza de salud 1
lx2e0(mlx2) - sum(MThx*rev(mdx2)) # esperanza de salud 2 # caso 2 es una mejora muy fuerte!

fga1 <- colSums(FThx*rev(fdx1))/flx1
fga2 <- colSums(FThx*rev(fdx2))/flx2
mga1 <- colSums(MThx*rev(mdx1))/mlx1
mga2 <- colSums(MThx*rev(mdx2))/mlx2

#---------------------plotting
Age <- 74:105
APadlf3 <- local(get(load("C:\\Dropbox\\sex decomp Tim\\APadlf3.Rdata"))) 
APadlm3 <- local(get(load("C:\\Dropbox\\sex decomp Tim\\APadlm3.Rdata"))) 
fmorbp1 <- APadlf3[10:41,4]
fmorbp2 <- APadlf3[10:41,19]
mmorbp1 <- APadlm3[10:41,4]
mmorbp2 <- APadlm3[10:41,19]

x11(height=6,width=5.5)

plot(Age,fga1[1:32],lwd=2,t="l",ylab="Disability prevalence g'(a)",
     main="Observed period mortality, \n Cohort 1915 g(y)",col="darkorange",ylim=c(0,0.45))
lines(Age,fga2[1:32],lwd=2,col="red")
lines(Age,mga1[1:32],lwd=2,col="cornflowerblue")
lines(Age,mga2[1:32],lwd=2,col="navyblue")
legend("topleft",legend=c("Female 1995", "Female 2010","Male 1995", "Male 2010"),
       lty=1,lwd=2,col=c("darkorange","red","cornflowerblue","navyblue"))


plot(Age,fmorbp1,lwd=2,t="l",ylab="Disability prevalence g'(a)",
     main="Observed period morbidity",col="darkorange",ylim=c(0,0.9))
lines(Age,fmorbp2,lwd=2,col="red")
lines(Age,mmorbp1,lwd=2,col="cornflowerblue")
lines(Age,mmorbp2,lwd=2,col="navyblue")


# -------------------------------------------Decomposing














# unhealthy life years 
geteUx <- function(mx,morb){
  lx <- mx %>% mx2qx %>% qx2lx 
  Lx <- lx2Lx(lx)
  Lx  <- Lx / lx[1]
  if (length(Lx) > length(morb)){
    morb <- c(morb, morb[length(morb)])
  }  
  sum(Lx*morb)
}

# healthy life years expectancy
geteHx <- function(mx, morb){
  lx  <- qx2lx(mx2qx(mx[-length(mx)]))
  lx2ex(lx) - geteUx(mx, morb)
}


#-----------------------------------------------------------
library(DecompHoriuchi)
library(magrittr)
mxgy2eH <- function(rates){
  mx   <- rates[1:41]
  gy   <- rates[-c(1:41)]
  Morb <- matrix(gy, length(gy), length(gy))
  geteHx(mx,Morb) 
}

fnogy <- DecompContinuousOrig(mxgy2eH, rates1 = c(fmxp1,Fgy), 
                              rates2 = c(fmxp2,Fgy),N=20)
fmxcontrib <- fnogy[1:41]
fgycontrib <- fnogy[-c(1:41)]




















#--#-------------------------
# some helper functions:
  
  
  mx2ex <- function(mx){
    mx %>% mx2qx %>% qx2lx %>% lx2ex
  }

# remaining life expectancies at age 65

fe74p1 <- mx2ex(fmxp1)
me74p1 <- mx2ex(mmxp1)

fe74p2 <- mx2ex(fmxp2)
me74p2 <- mx2ex(mmxp2)

# dx for the periods

fdxp1 <- lx2dx(qx2lx(mx2qx(fmxp1)))
mdxp1 <- lx2dx(qx2lx(mx2qx(mmxp1)))

fdxp2 <- lx2dx(qx2lx(mx2qx(fmxp2)))
mdxp2 <- lx2dx(qx2lx(mx2qx(mmxp2)))


# putting dx on a matrix
dxMp1      <- mdxp1[row(M1915ee) + col(M1915ee) - 1]
dim(dxMp1) <- dim(M1915ee)
dxFp1      <- fdxp1[row(F1915ee) + col(F1915ee) - 1]
dim(dxFp1) <- dim(F1915ee)

dxMp2      <- mdxp2[row(M1915ee) + col(M1915ee) - 1]
dim(dxMp2) <- dim(M1915ee)
dxFp2      <- fdxp2[row(F1915ee) + col(F1915ee) - 1]
dim(dxFp2) <- dim(F1915ee)


# verify that the column sums match lx:

flxp1 <- qx2lx(mx2qx(fmxp1))
mlxp1 <- qx2lx(mx2qx(mmxp1))
flxp2 <- qx2lx(mx2qx(fmxp2))
mlxp2 <- qx2lx(mx2qx(mmxp2))

# marginal difference (4th decimal) in the first entry
all(colSums(dxMp1,na.rm=TRUE) == mlxp1[-length(mlxp1)])
all(colSums(dxFp1,na.rm=TRUE) == flxp1[-length(flxp1)])
all(colSums(dxMp2,na.rm=TRUE) == mlxp2[-length(mlxp2)])
all(colSums(dxFp2,na.rm=TRUE) == flxp2[-length(flxp2)])

# so, we can multiply in neat ways:
Mtp1 <- dxMp1 * M1915ee
Ftp1 <- dxFp1 * F1915ee

Mtp2 <- dxMp2 * M1915ee
Ftp2 <- dxFp2 * F1915ee

# the column sums of this morbidity prevalence matrix give the chronological
# age pattern of prevalence within the cohort (Assuming this radix)

Mc <- colSums(Mtc, na.rm=TRUE) / lxm[-length(lxm)]
Fc <- colSums(Ftc, na.rm=TRUE) / lxf[-length(lxf)]

plot(74:110, Mc, type = 'l', ylim = c(0, .6), main = "marginal chronological age pattern (artifactual in this case)")
lines(74:110, Fc, col = "red")




# Part of the difference between males and females is due to the underlying
# morbidity pattern, which translates to the difference in the distribution
# of lifelines. 

# Further (what I tried to demonstrate in the presentation),
# if one holds the morbidity pattern constant, but changes
# the underlying mortality, we will draw different conclusions about the expected
# change in prevalence than if we simply assume the marginal chronological age pattern
# and then perturb mortality. 
# leverage gained by this added insight into morbidity variation could help
# better predict convergence/divergence between males and females, for instance.


#---------------------------------------------------------------------------------------
# code copied from my Colmex class code. Coments left in place.
# I might translate them sometime, who knows. Just need some
# quick plots.
# this script uses the HMDresults object to search for common patterns to the various defined measures.

# metodo barato de covertir m(x) en l(x)

######################################################
# ahora, definimos dos regimenes de mortalidad. Uno alto, el otro mas bajo.

# x = la edad
x   <- 0:100
# hacemos m(x)
mx1 <- mxgomp(a=.0005,b=.075,x)
mx2 <- mxgomp(a=.0005,b=.065,x) # el azul es mas bajo

# mirar las dos m(x)
plot(mx1,log='y', type = 'l')
lines(mx2,col = "blue")

# hacemos l(x) desde m(x)
lx1 <- mx2lx(mx1)
lx2 <- mx2lx(mx2)
# cuales son los e(0) implicados?
lx2e0(lx1);lx2e0(lx2) # diferencia de unos 7 años

# miramos l(x) nada raro aqui

pdf("Figures/LabPres/Z1PopAPopB.pdf",width=4.5,height=4.5)
par(mai=c(.9,.9,.1,.1))
plot(x,lx1,type='l', xlab = "chrono age", ylab = "l(x)")
lines(x,lx2,col = "blue")
legend("topright",col=c("black","blue"),lty=1,
       legend=c(paste0("l(x) A | e(0 )= ",round(lx2e0(lx1),1)),
                paste0("l(x) B | e(0) = ",round(lx2e0(lx2),1))),bty="n")
dev.off()
# hacemos d(x) para los dos casos
dx1 <- lx2dx(lx1)
dx2 <- lx2dx(lx2)

############################################################
# gran paso: imaginamos una pauta de alguna caracteristica
# que varie solo por edad cronologica y otra funcion de edad
# que es tanatologica.
# chrono and thano patterns:

# g(x) es alguna caracteristica cronologica (dolor de espalda?)

gx <- expit(seq(-6,4,length.out=101))
# t(y) es alguna caracteristica tanatologica (algun ADL?)
ty <- expit(seq(1,-80,length.out = 101))

# miramos g(x)

pdf("Figures/LabPres/Z2gx.pdf",width=4.5,height=4.5)
par(mai=c(.9,.9,.1,.1))
plot(x,gx,type='l', xlab = "chrono age", ylab = "prop with condition g")
dev.off()
pdf("Figures/LabPres/Z3ty.pdf",width=4.5,height=4.5)
par(mai=c(.9,.9,.1,.1))
plot(x,ty,type='l', xlab = "thano age", ylab = "prop with condition t")
dev.off()
# como cambia la prevalencia en la poblacion estacionaria
# si solo cambiamos el regimen de mortalidad. especificamente,
# si bajamos la mortalidad?

pdf("Figures/LabPres/Z4GxAB.pdf")
par(mai=c(.9,.9,.1,.1))
plot(x,lx1,type='l', xlab = "chrono age", ylab = "l(x) and prevalence, G(x)")
lines(x,gx,col="red")
polygon(c(x,rev(x)),c(rep(0,101),rev(gx*lx1)),col="#FF000050")
lines(x,lx2,col = "blue")
polygon(c(x,rev(x)),c(rep(0,101),rev(gx*lx2)),col="#0000FF50")
legend(73,.8,col=c("black","blue","red"),
       lty=c(1,1,1),
       legend=c("l(x) A","l(x) B","g(x)"),bty="n")
legend(73,.65,col=c("#FF000050","#0000FF50"),
       fill=c("#FF000050","#0000FF50"),
       legend=c("G(x) under A","G(x) under B"),bty="n")
dev.off()
# la respuesta: sube la prevalencia de condicion 'g' si es que la condicion es
# cronologica.

##########################################
# unos indices para el caso cronologico: #
##########################################
# Crono: aumenta de esperanza de vida aumenta la prevalencia de la condicion g
sum(gx * lx1) # esperanza de morbididad 1
sum(gx * lx2) # esperanza de morbididad 2 (mas grande en caso de variacion cronologica)
lx2e0(lx1) - sum(gx*lx1) # e0 saludabre 1
lx2e0(lx2) - sum(gx*lx2) # e0 saludabre 2 (ha subido, pero la morbididad ha subido mas!)

############################################################################
# ahora miramos que pasa si es que la caracteristica es funcion de la edad
# tanatologica SOLO.
############################################################################

# miramos t(y): solo sube al final de la vida: recuerda que el tiempo mueve a la izquierda!
par(mai=c(.9,.9,.1,.1))
plot(x, ty, type = 'l', xlab = "edad tanatologica",ylab =  "proporcion con condicion T")

#####################################
# para calcular, mejor tener t(y) en una matriz:
Thx <- matrix(0,101,101)
# columnas = crono ; filas = duracion de vida. despues ponderamos con d(x)
for (i in 1:101){
  Thx[i,1:(102-i)] <- ty[(102-i):1]
  #	if ( i > 1){
  #	Thx[i,(102-i+1):101] <- NA
}
#}
image(Thx)
# dibuja como es la pauta en la matriz para que se vea



plot(x,dx1,type='l')
lines(x,dx2,col="blue")

######################################
# como cambia la prevalencia de T bajo los dos regimenes de mortalidad?
pdf("Figures/LabPres/Z6TxATxB.pdf")
par(mai=c(.9,.9,.1,.1))
plot(x,lx1,type='l', xlab = "edad cronologica", ylab = "l(x), T(x), T(x)/l(x)")
polygon(c(x,rev(x)),c(rep(0,101),rowSums(Thx*rev(dx1))),
        col="#FF000050")

lines(x,lx2,col="blue")
polygon(c(x,rev(x)),c(rep(0,101),rowSums(Thx*rev(dx2))),col="#0000FF50")
lines(x,colSums(Thx*rev(dx1))/lx1, col = "black",lty=2)
lines(x, colSums(Thx * rev(dx2))/lx2, col = "blue",lty=2)

legend(65,1,col=c("black","blue","black","blue"),
       lty=c(1,1,2,2),
       legend=c("l(x) under A","l(x) under B","t'(x) under A","t'(x) under B"),bty="n")
legend(65,.83,col=c("#FF000050","#0000FF50","#00000050"),
       fill=c("#FF000050","#0000FF50","#00000050"),
       legend=c("T(x) under A","T(x) under B", "T'(x) under B|t'(x)A"),bty="n")
#polygon(c(x,rev(x)),c(rep(0,101),rev((colSums(Thx*rev(dx1))/lx1)*lx2)),col="#000000")
txproy <- colSums(Thx*rev(dx1))/lx1 * lx2
polygon(c(x,rev(x)),c(rep(0,101),rev(txproy)),col="#00000050")
dev.off()



# la prevalencia baja cuando la mortalidad baja!
####################################################
pdf("Figures/LabPres/Z5txAtxB.pdf",width=4.5,height=4.5)
par(mai=c(.9,.9,.1,.1))
plot(x,colSums(Thx*rev(dx1))/lx1, type = 'l',xlab = "chrono age", ylab = "prop with condition T")
lines(x, colSums(Thx * rev(dx2))/lx2, col = "blue")
legend("topleft",1,col=c("black","blue",NA,"#FF000050","#0000FF50"),
       lty=c(1,1),
       legend=c("t'(x) under A","t'(x) under B"),bty="n")
dev.off()
####################################################
# indices de resumen:
# cambio trivial en la esperanza de morbididad
sum(Thx*rev(dx1)) # esperanza de morbididad 1
sum(Thx*rev(dx2)) # esperanza de morbididad 2
lx2e0(lx1) - sum(Thx*rev(dx1)) # esperanza de salud 1
lx2e0(lx2) - sum(Thx*rev(dx2)) # esperanza de salud 2 # caso 2 es una mejora muy fuerte!
sum(txproy)
lx2e0(lx2)-sum(txproy)

100 * (1 - sum(Thx*rev(dx1)) / lx2e0(lx1)) # esperanza % saludabre
100 * (1 - sum(Thx*rev(dx2)) / lx2e0(lx2)) # algo mejor ---
# tampoco puedes esperar mucho mas aqui, estamos acercando a 100...



