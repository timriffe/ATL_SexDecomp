# Assuming stationary populations, examining change in morbidity and mortality components 
# given fixed gy, all HMD period 1x1 life tables, 10 years apart.

library(plyr)
library(dplyr)
options(scipen=3,stringsAsFactors=FALSE)
setwd("C:\\Dropbox\\sex decomp Tim\\Data\\")
load("gyall.Rdata")
load("SmArrays.Rdata")
PrevTTD$TTD <- as.numeric(PrevTTD$TTD)


# averaging over cohorts
meanprev <- ddply(PrevTTD,.(Morbidity,Sex,TTD), summarize, Morb=mean(Prev,na.rm=T))
head(DatForAlyson)

# adding data from DatForAlyson which has calculated proportion >= 1,2,3 ADLs and IADLs
SRHM <- rowMeans(DatForAlyson$SRHpoor$Males$Surf,na.rm=T)
SRHF <- rowMeans(DatForAlyson$SRHpoor$Females$Surf,na.rm=T)
ADL3M <- rowMeans(DatForAlyson$ADL3$Males$Surf,na.rm=T)
ADL3F <- rowMeans(DatForAlyson$ADL3$Females$Surf,na.rm=T)
IADL1M <- rowMeans(DatForAlyson$IADL1$Males$Surf,na.rm=T)
IADL1F <- rowMeans(DatForAlyson$IADL1$Females$Surf,na.rm=T)

ADL2M <- rowMeans(DatForAlyson$ADL2$Males$Surf,na.rm=T)
ADL2F <- rowMeans(DatForAlyson$ADL2$Females$Surf,na.rm=T)
ADL1M <- rowMeans(DatForAlyson$ADL1$Males$Surf,na.rm=T)
ADL1F <- rowMeans(DatForAlyson$ADL1$Females$Surf,na.rm=T)
IADL2M <- rowMeans(DatForAlyson$IADL2$Males$Surf,na.rm=T)
IADL2F <- rowMeans(DatForAlyson$IADL2$Females$Surf,na.rm=T)
IADL3M <- rowMeans(DatForAlyson$IADL3$Males$Surf,na.rm=T)
IADL3F <- rowMeans(DatForAlyson$IADL3$Females$Surf,na.rm=T)


# meanprev <- rbind(meanprev,
#                   data.frame(Morbidity="SRHpoor",Sex='m',TTD=0:12,Morb=SRHM),
#                   data.frame(Morbidity="SRHpoor",Sex='f',TTD=0:12,Morb=SRHF),
#                   data.frame(Morbidity="ADL3",Sex='m',TTD=0:12,Morb=ADL3M),
#                   data.frame(Morbidity="ADL3",Sex='f',TTD=0:12,Morb=ADL3F),
#                   data.frame(Morbidity="IADL1",Sex='m',TTD=0:12,Morb=IADL1M),
#                   data.frame(Morbidity="IADL1",Sex='f',TTD=0:12,Morb=IADL1F),
#                   data.frame(Morbidity="ADL2",Sex='m',TTD=0:12,Morb=ADL2M),
#                   data.frame(Morbidity="ADL2",Sex='f',TTD=0:12,Morb=ADL2F),
#                   data.frame(Morbidity="ADL1",Sex='m',TTD=0:12,Morb=ADL1M),
#                   data.frame(Morbidity="ADL1",Sex='f',TTD=0:12,Morb=ADL1F),
#                   data.frame(Morbidity="IADL2",Sex='m',TTD=0:12,Morb=IADL2M),
#                   data.frame(Morbidity="IADL2",Sex='f',TTD=0:12,Morb=IADL2F),
#                   data.frame(Morbidity="IADL3",Sex='m',TTD=0:12,Morb=IADL3M),
#                   data.frame(Morbidity="IADL3",Sex='f',TTD=0:12,Morb=IADL3F))
# 

meanprev <- rbind(meanprev,
                  data.frame(Morbidity="srh",Sex='m',TTD=0:12,Morb=SRHM),
                  data.frame(Morbidity="srh",Sex='f',TTD=0:12,Morb=SRHF),
                  data.frame(Morbidity="ADL3",Sex='m',TTD=0:12,Morb=ADL3M),
                  data.frame(Morbidity="ADL3",Sex='f',TTD=0:12,Morb=ADL3F),
                  data.frame(Morbidity="IADL1",Sex='m',TTD=0:12,Morb=IADL1M),
                  data.frame(Morbidity="IADL1",Sex='f',TTD=0:12,Morb=IADL1F),
                  data.frame(Morbidity="ADL2",Sex='m',TTD=0:12,Morb=ADL2M),
                  data.frame(Morbidity="ADL2",Sex='f',TTD=0:12,Morb=ADL2F),
                  data.frame(Morbidity="ADL1",Sex='m',TTD=0:12,Morb=ADL1M),
                  data.frame(Morbidity="ADL1",Sex='f',TTD=0:12,Morb=ADL1F),
                  data.frame(Morbidity="IADL2",Sex='m',TTD=0:12,Morb=IADL2M),
                  data.frame(Morbidity="IADL2",Sex='f',TTD=0:12,Morb=IADL2F),
                  data.frame(Morbidity="IADL3",Sex='m',TTD=0:12,Morb=IADL3M),
                  data.frame(Morbidity="IADL3",Sex='f',TTD=0:12,Morb=IADL3F))


# combining with data from decomp...

# WCres1 <- read.csv("C:\\Dropbox\\sex decomp Tim\\Data\\Within_country_decomp2.csv",head=T)
# BCres1 <- read.csv("C:\\Dropbox\\sex decomp Tim\\Data\\between_country_decomp2.csv",head=T)
# WCres2 <- read.csv("C:\\Dropbox\\sex decomp Tim\\Data\\Within_country_decomp3.csv",head=T)
# BCres2 <- read.csv("C:\\Dropbox\\sex decomp Tim\\Data\\between_country_decomp3.csv",head=T)
# 
# WCres <- rbind(WCres1,WCres2)
# BCres <- rbind(BCres1,BCres2)
# 
# # maximum prevalence by disability type
# maxpr <- ddply(meanprev,.(Sex,Morbidity),summarize,maxprev=max(Morb))
# colnames(maxpr)=c("Sex","Morbtype","maxprev")
# 
# # joining to WC and BC res
# WCres <- join(WCres,maxpr)
# BCres <- join(BCres,maxpr)
# 
# write.csv(WCres,file="C:\\Dropbox\\sex decomp Tim\\Data\\WCres.csv",row.names=F)
# write.csv(BCres,file="C:\\Dropbox\\sex decomp Tim\\Data\\BCres.csv",row.names=F)


plot(meanprev$TTD[meanprev$Morbidity=="nh_now" & meanprev$Sex=="f"][1:12],cex.axis=1.1,cex.lab=1.1,
     meanprev$Morb[meanprev$Morbidity=="nh_now" & meanprev$Sex=="f"][1:12],t="l",lwd=2,
     ylim=c(0,0.7),ylab="Disability prevalance",xlab="Time to Death (Years)")
lines(meanprev$TTD[meanprev$Morbidity=="srh" & meanprev$Sex=="f"][1:12],lwd=2,
      meanprev$Morb[meanprev$Morbidity=="srh" & meanprev$Sex=="f"][1:12],t="l",col="red")
lines(meanprev$TTD[meanprev$Morbidity=="name_mo" & meanprev$Sex=="f"][1:12],lwd=2,
      meanprev$Morb[meanprev$Morbidity=="name_mo" & meanprev$Sex=="f"][1:12],t="l",col="blue")
grid()
legend("topright",legend=c("In nursing home","Poor self-rated health","Unable to name month"),
       bty="n",col=c("black","red","blue"),
       lty=1,ncol=1,lwd=2,cex=0.8)





cairo_pdf("C:\\Dropbox\\sex decomp Tim\\figures\\DisbyTTD_pres.pdf",width=6,height=4)
par(mfrow=c(1,2))

plot(meanprev$TTD[meanprev$Morbidity=="ADL1" & meanprev$Sex=="f"],cex.axis=1.1,cex.lab=1.1,
     meanprev$Morb[meanprev$Morbidity=="ADL1" & meanprev$Sex=="f"],t="l",lwd=2,
     ylim=c(0,0.7),ylab="Disability prevalance",xlab="Time to Death (Years)")
lines(meanprev$TTD[meanprev$Morbidity=="ADL2" & meanprev$Sex=="f"],lwd=2,
      meanprev$Morb[meanprev$Morbidity=="ADL2" & meanprev$Sex=="f"],t="l",col="red")
lines(meanprev$TTD[meanprev$Morbidity=="ADL3" & meanprev$Sex=="f"],lwd=2,
      meanprev$Morb[meanprev$Morbidity=="ADL3" & meanprev$Sex=="f"],t="l",col="blue")
grid()
legend("topright",legend=c(paste0("ADL ","\u2265 1"," of 5" ),
                           paste0("ADL ","\u2265 2"," of 5" ),
                           paste0("ADL ","\u2265 3"," of 5")),
       bty="n",col=c("black","red","blue"),
       lty=1,ncol=1,lwd=2,cex=0.8)

plot(meanprev$TTD[meanprev$Morbidity=="nh_now" & meanprev$Sex=="f"][1:12],cex.axis=1.1,cex.lab=1.1,
     meanprev$Morb[meanprev$Morbidity=="nh_now" & meanprev$Sex=="f"][1:12],t="l",lwd=2,
     ylim=c(0,0.7),ylab="Disability prevalance",xlab="Time to Death (Years)")
lines(meanprev$TTD[meanprev$Morbidity=="SRHpoor" & meanprev$Sex=="f"][1:12],lwd=2,
      meanprev$Morb[meanprev$Morbidity=="SRHpoor" & meanprev$Sex=="f"][1:12],t="l",col="red")
lines(meanprev$TTD[meanprev$Morbidity=="name_mo" & meanprev$Sex=="f"][1:12],lwd=2,
      meanprev$Morb[meanprev$Morbidity=="name_mo" & meanprev$Sex=="f"][1:12],t="l",col="blue")
grid()
legend("topright",legend=c("In nursing home","Poor self-rated health","Unable to name month"),
       bty="n",col=c("black","red","blue"),
       lty=1,ncol=1,lwd=2,cex=0.8)

dev.off()



cairo_pdf("C:\\Dropbox\\sex decomp Tim\\figures\\DisbyTTD.pdf",width=8,height=4)
par(mfrow=c(1,3))

plot(meanprev$TTD[meanprev$Morbidity=="IADL1" & meanprev$Sex=="f"],cex.axis=1.5,cex.lab=1.5,
     meanprev$Morb[meanprev$Morbidity=="IADL1" & meanprev$Sex=="f"],t="l",lwd=2,
     ylim=c(0,0.7),ylab="Disability prevalance",xlab="Time to Death (Years)")
lines(meanprev$TTD[meanprev$Morbidity=="IADL2" & meanprev$Sex=="f"],lwd=2,
      meanprev$Morb[meanprev$Morbidity=="IADL2" & meanprev$Sex=="f"],t="l",col="red")
lines(meanprev$TTD[meanprev$Morbidity=="IADL3" & meanprev$Sex=="f"],lwd=2,
      meanprev$Morb[meanprev$Morbidity=="IADL3" & meanprev$Sex=="f"],t="l",col="blue")
grid()
legend("topright",legend=c(paste0("IADL ","\u2265 1"," of 5" ),paste0("IADL ","\u2265 2"," of 5" ),
                           paste0("IADL ","\u2265 3"," of 5" )),
       bty="n",col=c("black","red","blue"),
       lty=1,ncol=1,lwd=2,cex=1.2)



plot(meanprev$TTD[meanprev$Morbidity=="ADL1" & meanprev$Sex=="f"],cex.axis=1.5,cex.lab=1.5,
     meanprev$Morb[meanprev$Morbidity=="ADL1" & meanprev$Sex=="f"],t="l",lwd=2,
     ylim=c(0,0.7),ylab="Disability prevalance",xlab="Time to Death (Years)")
lines(meanprev$TTD[meanprev$Morbidity=="ADL2" & meanprev$Sex=="f"],lwd=2,
     meanprev$Morb[meanprev$Morbidity=="ADL2" & meanprev$Sex=="f"],t="l",col="red")
lines(meanprev$TTD[meanprev$Morbidity=="ADL3" & meanprev$Sex=="f"],lwd=2,
     meanprev$Morb[meanprev$Morbidity=="ADL3" & meanprev$Sex=="f"],t="l",col="blue")
grid()
legend("topright",legend=c(paste0("ADL ","\u2265 1"," of 5" ),
                           paste0("ADL ","\u2265 2"," of 5" ),
                           paste0("ADL ","\u2265 3"," of 5")),
       bty="n",col=c("black","red","blue"),
       lty=1,ncol=1,lwd=2,cex=1.2)

plot(meanprev$TTD[meanprev$Morbidity=="nh_now" & meanprev$Sex=="f"][1:12],cex.axis=1.5,cex.lab=1.5,
     meanprev$Morb[meanprev$Morbidity=="nh_now" & meanprev$Sex=="f"][1:12],t="l",lwd=2,
     ylim=c(0,0.7),ylab="Disability prevalance",xlab="Time to Death (Years)")
lines(meanprev$TTD[meanprev$Morbidity=="SRHpoor" & meanprev$Sex=="f"][1:12],lwd=2,
      meanprev$Morb[meanprev$Morbidity=="SRHpoor" & meanprev$Sex=="f"][1:12],t="l",col="red")
lines(meanprev$TTD[meanprev$Morbidity=="name_mo" & meanprev$Sex=="f"][1:12],lwd=2,
      meanprev$Morb[meanprev$Morbidity=="name_mo" & meanprev$Sex=="f"][1:12],t="l",col="blue")
grid()
legend("topright",legend=c("In nursing home","Poor self-rated health","Unable to name month"),
       bty="n",col=c("black","red","blue"),
       lty=1,ncol=1,lwd=2,cex=1.2)

dev.off()




cairo_pdf("C:\\Dropbox\\sex decomp Tim\\figures\\DisbyTTD2.pdf",width=6,height=6)

plot(meanprev$TTD[meanprev$Morbidity=="IADL1" & meanprev$Sex=="f"],
     meanprev$Morb[meanprev$Morbidity=="IADL1" & meanprev$Sex=="f"],t="l",lwd=2,col="navyblue",
     ylim=c(0,0.7),ylab="Disability prevalance",xlab="Time to death")
lines(meanprev$TTD[meanprev$Morbidity=="IADL2" & meanprev$Sex=="f"],lwd=2,
      meanprev$Morb[meanprev$Morbidity=="IADL2" & meanprev$Sex=="f"],t="l",col="blue")
lines(meanprev$TTD[meanprev$Morbidity=="IADL3" & meanprev$Sex=="f"],lwd=2,
      meanprev$Morb[meanprev$Morbidity=="IADL3" & meanprev$Sex=="f"],t="l",col="deepskyblue")

lines(meanprev$TTD[meanprev$Morbidity=="ADL1" & meanprev$Sex=="f"],col="red4",
     meanprev$Morb[meanprev$Morbidity=="ADL1" & meanprev$Sex=="f"],t="l",lwd=2)
lines(meanprev$TTD[meanprev$Morbidity=="ADL2" & meanprev$Sex=="f"],lwd=2,
      meanprev$Morb[meanprev$Morbidity=="ADL2" & meanprev$Sex=="f"],t="l",col="red3")
lines(meanprev$TTD[meanprev$Morbidity=="ADL3" & meanprev$Sex=="f"],lwd=2,
      meanprev$Morb[meanprev$Morbidity=="ADL3" & meanprev$Sex=="f"],t="l",col="red")


lines(meanprev$TTD[meanprev$Morbidity=="nh_now" & meanprev$Sex=="f"][1:13],
     meanprev$Morb[meanprev$Morbidity=="nh_now" & meanprev$Sex=="f"][1:13],t="l",lwd=2)     
lines(meanprev$TTD[meanprev$Morbidity=="SRHpoor" & meanprev$Sex=="f"][1:13],lwd=2,
      meanprev$Morb[meanprev$Morbidity=="SRHpoor" & meanprev$Sex=="f"][1:13],t="l",col="green3")
lines(meanprev$TTD[meanprev$Morbidity=="name_mo" & meanprev$Sex=="f"][1:13],lwd=2,
      meanprev$Morb[meanprev$Morbidity=="name_mo" & meanprev$Sex=="f"][1:13],t="l",col="grey80")


legend("topright",legend=c(paste0("IADL ","\u2265 1"," of 5" ),paste0("IADL ","\u2265 2"," of 5" ),
                           paste0("IADL ","\u2265 3"," of 5" ),paste0("ADL ","\u2265 1"," of 5" ),
                           paste0("ADL ","\u2265 2"," of 5" ),paste0("ADL ","\u2265 3"," of 5" ),
                           "In nursing home","Poor self-rated health","Unable to name month"),
       bty="n",col=c("navyblue","blue","deepskyblue","red4","red3","red","black","green3","grey80"),
       lty=1,ncol=1,lwd=2)

dev.off()



#---------------------------------------------


Decomp.DLY <- function(Ctry,Year1,Year2,sex,Morbtype){
  # Ctry is HMD short form, sex is 'm' or 'f', Morbtypes can be seen unique(meanprev$Morbidity)
  Morbdat <- filter(meanprev,Morbidity==Morbtype,Sex==sex)

  LTF <- read.table(paste("C:\\hmd_statistics\\lt_female\\fltper_1x1\\",Ctry,
                        ".fltper_1x1.txt",sep=""),head=T,skip=2,na.strings=".")
  LTM <- read.table(paste("C:\\hmd_statistics\\lt_male\\mltper_1x1\\",Ctry,
                        ".mltper_1x1.txt",sep=""),head=T,skip=2,na.strings=".")
             
  if(sex=="m")
      LT = LTM  else LT = LTF
  
  LT$Age <- rep(0:110,length(unique(LT$Year)))
  
  LT1 <- filter(LT,Year==Year1,Age>=60)
  LT2 <- filter(LT,Year==Year2,Age>=60)
  
  # vector of remaining lifespans  
  ages <- 60:110
  rls1 <- list()   
  rls2 <- list()
    for (i in 1:length(ages)){
      rls1[[i]] <- (LT1$dx[LT1$Age>=ages[i]])/sum(LT1$dx[LT1$Age>=ages[i]])  
      rls2[[i]] <- (LT2$dx[LT2$Age>=ages[i]])/sum(LT2$dx[LT2$Age>=ages[i]])  
    }
  
  # adding NAs to make each vector the same length
  rlswNA1 <- lapply(rls1, function(x) {c(x, rep(NA, length(ages) - length(x)))})
  rlswNA2 <- lapply(rls2, function(x) {c(x, rep(NA, length(ages) - length(x)))})
  
  # putting into matrix format
  rlsmat1 <- do.call(rbind, rlswNA1)
  rlsmat2 <- do.call(rbind, rlswNA2)
  
  # disability prevalence by age
  DP1 <- colSums(t(rlsmat1[,1:12])*Morbdat$Morb[1:12],na.rm=T)
  DP2 <- colSums(t(rlsmat2[,1:12])*Morbdat$Morb[1:12],na.rm=T)
  
  # at old ages where there are zero deaths DP is zero. Instead let's replace that with the 
  # maxixum disability, which is generally at the highest age with death counts
#  DP1[DP1==0] <- max(DP1)
#  DP2[DP2==0] <- max(DP2)
#   oldages <- 61:111
#   plot(oldages,DP1,t="l",lwd=2)
#   lines(oldages,DP2,col="red",lwd=2)
  
  # Person-years lived from life tables
  Lx1 <- LT1$Lx/100000
  Lx2 <- LT2$Lx/100000
  
  # person-years with disability at age x
  PYD1 <- Lx1*DP1
  PYD2 <- Lx2*DP2
  
  # change in person-years with disability over time
  dPYD <- (DP1+DP2)/2*(Lx2-Lx1)+(DP2-DP1)*(Lx1+Lx2)/2
  
  # decompose difference in person-years with disability between time 1 and time 2 
  # into mortality and disability prevalence component # from Nusselder and Looman (2004)
  MOR.C <- (DP1+DP2)/2*(Lx2-Lx1)
  DIS.C <- (Lx1+Lx2)/2*(DP2-DP1)
  
  res <- data.frame(Ctry=Ctry,Year1=Year1,Year2=Year2,total.diff=sum(dPYD),mortality.comp=sum(MOR.C),disability.comp=sum(DIS.C),
           e60.T1=LT1$ex[1],e60.T2=LT2$ex[1],e60.change=LT2$ex[1]-LT1$ex[1],
           DLY.T1=sum(PYD1),DLY.T2=sum(PYD2),Sex=sex,Morbtype=Morbtype)
  return(res)
  }


#--------------running on all ctries

Ctries <- read.csv("country_codes_HMD.csv",head=T)
Ctry <- Ctries$Ctry


# within country comparison
Ctrydf <- list()
for(c in 1:length(Ctry)){
  Ctrydf[[c]] <- data.frame(Ctry=Ctries$Ctry[c],
                            Year1=seq(Ctries$fyear10[c],Ctries$lyear10[c]-10,10),
                            Year2=seq(Ctries$fyear10[c]+10,Ctries$lyear10[c],10))
}
Ctrydf <- rbind.fill(Ctrydf)
expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))
varchoose <- expand.grid.df(Ctrydf,data.frame(Sex=c("f","m")),
                            data.frame(morbtype=c("ADL3","IADL1","nh_now","name_mo","SRHpoor")))
varchoose <- filter(varchoose,Year1>=1950)


Bigres <- list()
for (i in 1:dim(varchoose)[1]){
  Bigres[[i]] <- Decomp.DLY(Ctry=varchoose[i,1],Year1=varchoose[i,2],Year2=varchoose[i,3],
                            sex=varchoose[i,4],Morbtype=varchoose[i,5])
}
Bigres <- rbind.fill(Bigres)

write.csv(Bigres,file="Within_country_decomp.csv2",row.names=F)

#-----------------Between-country comparison
            #------ need to alter function here----------

Decomp.DLY.BC <- function(CtryA,CtryB,Year1,sex,Morbtype){
  # Ctry is HMD short form, sex is 'm' or 'f', Morbtypes can be seen unique(meanprev$Morbidity)
  Morbdat <- filter(meanprev,Morbidity==Morbtype,Sex==sex)
  
  LTAF <- read.table(paste("C:\\hmd_statistics\\lt_female\\fltper_1x1\\",CtryA,
                          ".fltper_1x1.txt",sep=""),head=T,skip=2,na.strings=".")
  LTAM <- read.table(paste("C:\\hmd_statistics\\lt_male\\mltper_1x1\\",CtryA,
                          ".mltper_1x1.txt",sep=""),head=T,skip=2,na.strings=".")
  
  if(sex=="m")
    LTA = LTM  else LTA = LTF
  
  LTA$Age <- rep(0:110,length(unique(LTA$Year)))
  
  
  LTBF <- read.table(paste("C:\\hmd_statistics\\lt_female\\fltper_1x1\\",CtryB,
                           ".fltper_1x1.txt",sep=""),head=T,skip=2,na.strings=".")
  LTBM <- read.table(paste("C:\\hmd_statistics\\lt_male\\mltper_1x1\\",CtryB,
                           ".mltper_1x1.txt",sep=""),head=T,skip=2,na.strings=".")
  
  if(sex=="m")
    LTB = LTM  else LTB = LTF
  
  LTB$Age <- rep(0:110,length(unique(LTB$Year)))
  
  LT1 <- filter(LTA,Year==Year1,Age>=60)
  LT2 <- filter(LTB,Year==Year1,Age>=60)
  
  # vector of remaining lifespans  
  ages <- 60:110
  rls1 <- list()   
  rls2 <- list()
  for (i in 1:length(ages)){
    rls1[[i]] <- (LT1$dx[LT1$Age>=ages[i]])/sum(LT1$dx[LT1$Age>=ages[i]])  
    rls2[[i]] <- (LT2$dx[LT2$Age>=ages[i]])/sum(LT2$dx[LT2$Age>=ages[i]])  
  }
  
  # adding NAs to make each vector the same length
  rlswNA1 <- lapply(rls1, function(x) {c(x, rep(NA, length(ages) - length(x)))})
  rlswNA2 <- lapply(rls2, function(x) {c(x, rep(NA, length(ages) - length(x)))})
  
  # putting into matrix format
  rlsmat1 <- do.call(rbind, rlswNA1)
  rlsmat2 <- do.call(rbind, rlswNA2)
  
  # disability prevalence by age
  DP1 <- colSums(t(rlsmat1[,1:12])*Morbdat$Morb[1:12],na.rm=T)
  DP2 <- colSums(t(rlsmat2[,1:12])*Morbdat$Morb[1:12],na.rm=T)
  
  # at old ages where there are zero deaths DP is zero. Instead let's replace that with the 
  # maxixum disability, which is generally at the highest age with death counts
  #  DP1[DP1==0] <- max(DP1)
  #  DP2[DP2==0] <- max(DP2)
  #   oldages <- 61:111
  #   plot(oldages,DP1,t="l",lwd=2)
  #   lines(oldages,DP2,col="red",lwd=2)
  
  # Person-years lived from life tables
  Lx1 <- LT1$Lx/100000
  Lx2 <- LT2$Lx/100000
  
  # person-years with disability at age x
  PYD1 <- Lx1*DP1
  PYD2 <- Lx2*DP2
  
  # change in person-years with disability over time
  dPYD <- (DP1+DP2)/2*(Lx2-Lx1)+(DP2-DP1)*(Lx1+Lx2)/2
  
  # decompose difference in person-years with disability between time 1 and time 2 
  # into mortality and disability prevalence component # from Nusselder and Looman (2004)
  MOR.C <- (DP1+DP2)/2*(Lx2-Lx1)
  DIS.C <- (Lx1+Lx2)/2*(DP2-DP1)
  
  res <- data.frame(CtryA=CtryA,CtryB=CtryB,Year=Year1,total.diff=sum(dPYD),mortality.comp=sum(MOR.C),disability.comp=sum(DIS.C),
                    e60.T1=LT1$ex[1],e60.T2=LT2$ex[1],e60.change=LT2$ex[1]-LT1$ex[1],
                    DLY.T1=sum(PYD1),DLY.T2=sum(PYD2),Sex=sex,Morbtype=Morbtype)
  return(res)
}


ccomp <- as.data.frame(t(combn(Ctry,2)))
names(ccomp) <- c("CtryA","CtryB")
varchoose2 <- expand.grid.df(ccomp,data.frame(Sex=c('m','f')),
                          data.frame(Year=seq(1980,2000,10)),
                          data.frame(morbtype=c("ADL3","IADL1","nh_now","name_mo","SRHpoor")))

Bigres2 <- list()
for (i in 1:dim(varchoose2)[1]){
  Bigres2[[i]] <- Decomp.DLY.BC(CtryA=varchoose2[i,1],CtryB=varchoose2[i,2],
                             Year1=varchoose2[i,4],sex=varchoose2[i,3],
                             Morbtype=varchoose2[i,5])
}
Bigres2 <- rbind.fill(Bigres2)

write.csv(Bigres2,file="between_country_decomp2.csv",row.names=F)



#-----------plotting some associations [see file 'plot decomps HMD']

#---IADL5
fresWC <- filter(Bigres,Morbtype=="iadl5_",Sex=="f")
mresWC <- filter(Bigres,Morbtype=="iadl5_",Sex=="m")

fresBC <- filter(Bigres2,Morbtype=="iadl5_",Sex=="f")
mresBC <- filter(Bigres2,Morbtype=="iadl5_",Sex=="m")

# x11()
# par(mfrow=c(3,2))
pdf(file="IADL_WC.pdf")
plot(fresWC$e60.change,fresWC$disability.comp,main="IADL 5 \n Within-country, 10 yrs apart",
     xlab=bquote(paste("Increase in ",e[60])),ylab="Change in disability component",
     xlim=range(c(fresWC$e60.change,mresWC$e60.change)),
     ylim=range(c(fresWC$disability.comp,mresWC$disability.comp)),col="darkred",
     cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
points(mresWC$e60.change,mresWC$disability.comp,col="cornflowerblue",pch=2)
abline(h=0,lty="dotted")
abline(v=0,lty="dotted")
legend("topright",c("Female","Male"),col=c("darkred","cornflowerblue"),pch=1:2,cex=1.5,bty="n")
dev.off()

# summary(lm(fresWC$disability.comp ~ fresWC$e60.change))
# summary(lm(mresWC$disability.comp ~ mresWC$e60.change))

pdf(file="IADL_WCnBC_f.pdf")
plot(fresWC$e60.change,fresWC$disability.comp,main="IADL 5 \n Females, WC and BC comparison",
     xlab=bquote(paste("Increase in ",e[60])),ylab="Change in disability component",
     xlim=range(c(fresWC$e60.change,fresBC$e60.change)),
     ylim=range(c(fresWC$disability.comp,fresBC$disability.comp)),col="darkred",
     cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
points(fresBC$e60.change,fresBC$disability.comp,col="darkorange",pch=2)
abline(h=0,lty="dotted")
abline(v=0,lty="dotted")
legend("topright",c("Within-country","Between-country"),col=c("darkred","darkorange"),pch=1:2,cex=1.5,bty="n")
dev.off()




plot(res1950$e60.change,res1950$disability.comp,main="IADL 5, since 1950",
     xlab=bquote(paste("Increase in ",e[60])),ylab="Change in disability component",
     xlim=range(res$e60.change),ylim=range(res$disability.comp))
abline(h=0,lty="dotted")
abline(v=0,lty="dotted")
text(2,0.1,bquote(paste(R^2,"=0.978")))
summary(lm(res1950$disability.comp ~ res1950$e60.change))


plot(res$e60.change,res$mortality.comp,main="IADL 5, all points",
     xlab=bquote(paste("Increase in ",e[60])),ylab="Change in mortality component",
     xlim=range(res$e60.change),ylim=range(res$mortality.comp))
abline(h=0,lty="dotted")
abline(v=0,lty="dotted")
text(-1,0.5,bquote(paste(R^2,"=0.838")))
summary(lm(res$mortality.comp ~ res$e60.change))



plot(res1950$e60.change,res1950$mortality.comp,main="IADL 5, since 1950",
     xlab=bquote(paste("Increase in ",e[60])),ylab="Change in mortality component",
     xlim=range(res$e60.change),ylim=range(res$mortality.comp))
abline(v=0,lty="dotted")
text(-1,0.3,bquote(paste(R^2,"=0.969")))
summary(lm(res1950$mortality.comp ~ res1950$e60.change))




plot(res$mortality.comp[res$e60.change<0],res$disability.comp[res$e60.change<0],
     col="darkorange",xlab="change in DLY owing to mortality",
     ylab="change in DLY owing to morbidity",
     xlim=range(res$mortality.comp),ylim=range(res$disability.comp),
     main="All points")
points(res$mortality.comp[res$e60.change>=0 & res$e60.change<1],
       res$disability.comp[res$e60.change>=0 & res$e60.change<1],
       col="red",pch=2)
points(res$mortality.comp[res$e60.change>=1 & res$e60.change<4],
       res$disability.comp[res$e60.change>=1 & res$e60.change<4],
       col="navyblue",pch=5)
abline(h=0,lty="dotted")
abline(v=0,lty="dotted")
# legend("topright",legend=c("change in e60 < 0","change in e60 > 0 < 1",
#                            "change in e60 > 1"), pch=c(1,2,5),
#           col=c("darkorange","red","navyblue"),bty="n")
# 


plot(res1950$mortality.comp[res1950$e60.change<0],res1950$disability.comp[res1950$e60.change<0],
    col="darkorange",xlab="change in DLY owing to mortality",
    ylab="change in DLY owing to morbidity",
    xlim=range(res$mortality.comp),ylim=range(res$disability.comp),
    main="Points since 1950")
points(res1950$mortality.comp[res1950$e60.change>=0 & res1950$e60.change<1],
       res1950$disability.comp[res1950$e60.change>=0 & res1950$e60.change<1],
       col="red",pch=2)
points(res1950$mortality.comp[res1950$e60.change>=1 & res1950$e60.change<4],
       res1950$disability.comp[res1950$e60.change>=1 & res1950$e60.change<4],
       col="navyblue",pch=5)
abline(h=0,lty="dotted")
abline(v=0,lty="dotted")
legend("topright",legend=c("change in e60 < 0","change in e60 > 0 < 1",
                           "change in e60 > 1"), pch=c(1,2,5),
       col=c("darkorange","red","navyblue"),bty="n")
