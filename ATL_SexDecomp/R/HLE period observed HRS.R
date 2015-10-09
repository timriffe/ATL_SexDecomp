# US data decompo
library(magrittr)

#--data
# These data are produced as follows: 
# Individuals are asked if they have problems with ADLs (either a list of 3 or 5).
# The proportions above are the total number of problems divided by 3 or 5 respectively


# females
APadlf3 <- local(get(load("C:\\Dropbox\\sex decomp Tim\\APadlf3.Rdata"))) 
#APadlf5 <- local(get(load("C:\\Dropbox\\sex decomp Tim\\APadlf5.Rdata"))) 

# males
APadlm3 <- local(get(load("C:\\Dropbox\\sex decomp Tim\\APadlm3.Rdata"))) 
#APadlm5 <- local(get(load("C:\\Dropbox\\sex decomp Tim\\APadlm5.Rdata"))) 


#---let's use ADL 3 to start with

# Visualizing the trends for different age groups - 69,79,89,99
# pretty flat for most age groups, maybe increasing for oldest women,
# and up and downs for 89 yr olds
Year <- 1992:2013

plot(Year,APadlf3[5,],ylim=c(0,0.6),t="o",lwd=2)
lines(Year,APadlf3[15,],pch=2,t="o",lwd=2)
lines(Year,APadlf3[25,],pch=4,t="o",lwd=2)
lines(Year,APadlf3[35,],pch=5,t="o",lwd=2)

plot(Year,APadlm3[5,],ylim=c(0,0.6),t="o",lwd=2)
lines(Year,APadlm3[15,],pch=2,t="o",lwd=2)
lines(Year,APadlm3[25,],pch=4,t="o",lwd=2)
lines(Year,APadlm3[35,],pch=5,t="o",lwd=2)

# First we want to see using the 'traditional method', how HLE has changed
# and decompose this difference in morbidity and mortality

# Period 1: 1995.  Period 2: 2010

# HMD LT mx for the two periods

USAFLT <- read.table("C:\\hmd_statistics\\lt_female\\fltper_1x1\\USA.fltper_1x1.txt",
                     head=T,skip=2,na.strings=".")
USAMLT <-  read.table("C:\\hmd_statistics\\lt_male\\mltper_1x1\\USA.mltper_1x1.txt",
                      head=T,skip=2,na.strings=".")
USAFLT[,2] <- rep(0:110,length.out=dim(USAFLT)[1])
USAMLT[,2] <- rep(0:110,length.out=dim(USAMLT)[1])

fmxp1 <- subset(USAFLT,Year==1995 & Age>=74 & Age<=105)$mx
fmxp2 <- subset(USAFLT,Year==2010 & Age>=74 & Age<=105)$mx

mmxp1 <- subset(USAMLT,Year==1995 & Age>=74 & Age<=105)$mx
mmxp2 <- subset(USAMLT,Year==2010 & Age>=74 & Age<=105)$mx


#-------------------------
# some helper functions:
  mx2qx <- function(mx){
    mx / (1 + .5 * mx)
  }
  
  
  qx2lx <- function(qx){
    cumprod(1-c(0,qx))
  }
  
  lx2dx <- function(lx){
    -diff(c(lx,0))
  }
  
  lx2Lx <- function(lx){
    (lx + c(lx[-1],0)) / 2
    
  }

  lx2ex <- function(lx){
    lx <- lx / lx[1]
    Lx <- lx2Lx(lx)
    sum(Lx)
  }

  mx2ex <- function(mx){
    mx %>% mx2qx %>% qx2lx %>% lx2ex
  }

# remaining life expectancies at age 74

fe74p1 <- mx2ex(fmxp1)
me74p1 <- mx2ex(mmxp1)

fe74p2 <- mx2ex(fmxp2)
me74p2 <- mx2ex(mmxp2)


# morbidity prevalence at ages 74 to 105, 1995 and 2010
fmorbp1 <- APadlf3[10:41,4]
fmorbp2 <- APadlf3[10:41,19]

mmorbp1 <- APadlm3[10:41,4]
mmorbp2 <- APadlm3[10:41,19]

fLx1 <- lx2Lx(qx2lx(mx2qx(fmx1)))


# unhealthy life years at age 74
  geteUx <- function(mx,morb){
    lx <- mx %>% mx2qx %>% qx2lx 
    Lx <- lx2Lx(lx)
    Lx  <- Lx / lx[1]
    if (length(Lx) > length(morb)){
      morb <- c(morb, morb[length(morb)])
    }  
    sum(Lx*morb)
  }



# healthy expectancy
# likewise we build from mx
  geteHx <- function(mx, morb){
    lx  <- qx2lx(mx2qx(mx[-length(mx)]))
    lx2ex(lx) - geteUx(mx, morb)
  }



# proportion of life healthy
  getPropH <- function(mx, morb){
    geteHx(mx,morb) / lx2ex(qx2lx(mx2qx(mx)))
  }


geteUx(fmxp1,fmorbp1)
geteHx(fmxp1,fmorbp1)
fe74p1
getPropH(fmxp1,fmorbp1)

geteUx(fmxp2,fmorbp2)
geteHx(fmxp2,fmorbp2)
fe74p2
getPropH(fmxp2,fmorbp2)

geteUx(mmxp1,mmorbp1)
geteHx(mmxp1,mmorbp1)
me74p1
getPropH(mmxp1,mmorbp1)

geteUx(mmxp2,mmorbp2)
geteHx(mmxp2,mmorbp2)
me74p2
getPropH(mmxp2,mmorbp2)

#------------------------------------------------------------------------