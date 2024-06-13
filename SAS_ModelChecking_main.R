rm(list=ls()) 
  
#######################################################################################
#definition of some quantities  

#ns       simulation times
#alpha    significant level
#N        full data size
#n        subsample data size
#n0       sample size for pilot study
#p        dimension of covarites
#delta    the Signal strength under alternative
#dm       the structural dimension under null(model dimension)
#d        the structural dimension of residuals
#h        the bandwidth for KDE f()
#hE       the bandwidth for estimation of E(\vareplsion|t(B)X)
#cE       parameter in estimating bandwidth
#c        parameter in testing bandwidth
#fx       sampling distribution at each sample point
#Lc       control limit(critical value)
#Model    Model type under H0
#struc    covariate structure
#eps      the distribution of error term
#replace  sampling with or without replacement
#wb       wild bootstrap time for Guo(2016) and SZ(2002)
#ker      kernel function in nonparametric estimation
#sdr      sufficient dimension reduction method 


#######################################################################################
#R package
library(MASS)
library(dr)
library(mgcv)
library(flexclust)
library(MAVE)
library(Matrix)
library(KernSmooth)
library(locfit)
#library(tvReg)
library(ggplot2)
#library(np)
# library(GrassmannOptim)
# library(ldr)
#library(ICtest)

#######################################################################################
######oracle, SAS and US##############################################################

struc <- "IID"
eps <- "norm"
sdr <- "MAVE"
ker <- "Epanechnikov"
c <- 0.1
cE <- 0.5
delta <- seq(0,1,0.125)

Res_temp <- rep(0,length(delta))
TypeI <- list(Res_US = Res_temp, Res_SAS_n0200 = Res_temp,
              Res_SAS_n0300 = Res_temp,Res_SAS_n0400 = Res_temp,
              Res_SAS_n0500 = Res_temp, Res_SAS_n0600 = Res_temp,
              Res_SAS_n01000 = Res_temp, Res_OS = Res_temp,
              Res_US0 = Res_temp)


tic <- proc.time()
#sm="US",n+n0 random sampling
func0 <- function(x){Process_SAS(delta[x],struc,eps,n0=300,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="US",Lc,ns)$TypeI}
TypeI$Res_US <- sapply(1:length(delta),func0)

#sm="SAS", n0=200
func200 <- function(x){Process_SAS(delta[x],struc,eps,n0=200,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="SAS",Lc,ns)$TypeI}
TypeI$Res_SAS_n0200 <- sapply(1:length(delta),func200)
#sm="SAS", n0=300
func300 <- function(x){Process_SAS(delta[x],struc,eps,n0=300,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="SAS",Lc,ns)$TypeI}
TypeI$Res_SAS_n0300 <- sapply(1:length(delta),func300)
#sm="SAS", n0=400
func400 <- function(x){Process_SAS(delta[x],struc,eps,n0=400,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="SAS",Lc,ns)$TypeI}
TypeI$Res_SAS_n0400 <- sapply(1:length(delta),func400)
#sm="SAS", n0=500
func500 <- function(x){Process_SAS(delta[x],struc,eps,n0=500,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="SAS",Lc,ns)$TypeI}
TypeI$Res_SAS_n0500 <- sapply(1:length(delta),func500)
#sm="SAS", n0=600
func600 <- function(x){Process_SAS(delta[x],struc,eps,n0=600,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="SAS",Lc,ns)$TypeI}
TypeI$Res_SAS_n0600 <- sapply(1:length(delta),func600)

#sm="OS", oracle sampling
func2 <- function(x){Process_SAS(delta[x],struc,eps,n0=300,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="OS",Lc,ns)$TypeI}
TypeI$Res_OS <- sapply(1:length(delta),func2)

#sm="US0",n random sampling
func00 <- function(x){Process_SAS(delta[x],struc,eps,n0=300,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="US0",Lc,ns)$TypeI}
TypeI$Res_US0 <- sapply(1:length(delta),func00)

tic2 <- proc.time()-tic
print(paste("Sam",Model,"-n",n,"-c",c,"-time:",tic2[3],sep=""))

TypeI
TypeI_thetatrue <- TypeI
write.csv(TypeI,file='TypeI.csv')
write.csv(TypeI_thetatrue,file='TypeI_thetatrue.csv')


#######################################################################################
#####testing bandwidth for SAS ########################################################

struc <- "IID"
eps <- "norm"
sdr <- "MAVE"
ker <- "Epanechnikov"

delta <- seq(0,0.5,0.25)
c <- c(0,0.1,0.2,0.3,0.4)
c <- c(0.05,0.1,0.15,0.2,0.25)
cE <- rep(0.5,5)

Res_temp <- matrix(0,length(delta),length(c))
TypeI <- list(Res_SAS_n0200 = Res_temp,
              Res_SAS_n0300 = Res_temp,
              Res_SAS_n0400 = Res_temp,
              Res_SAS_n0500 = Res_temp,
              Res_SAS_n0600 = Res_temp)

for (i in 1:length(c)){
  tic <- proc.time()

  #sm="SAS", n0=200
  func200 <- function(x){Process_SAS(delta[x],struc,eps,n0=200,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE[i],c[i],sm="SAS",Lc,ns)$TypeI}
  TypeI$Res_SAS_n0200[,i] <- sapply(1:length(delta),func200)

  #sm="SAS", n0=300
  func300 <- function(x){Process_SAS(delta[x],struc,eps,n0=300,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE[i],c[i],sm="SAS",Lc,ns)$TypeI}
  TypeI$Res_SAS_n0300[,i] <- sapply(1:length(delta),func300)

  #sm="SAS", n0=400
  func400 <- function(x){Process_SAS(delta[x],struc,eps,n0=400,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE[i],c[i],sm="SAS",Lc,ns)$TypeI}
  TypeI$Res_SAS_n0400[,i] <- sapply(1:length(delta),func400)

  #sm="SAS", n0=500
  func500 <- function(x){Process_SAS(delta[x],struc,eps,n0=500,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE[i],c[i],sm="SAS",Lc,ns)$TypeI}
  TypeI$Res_SAS_n0500[,i] <- sapply(1:length(delta),func500)

  #sm="SAS", n0=600
  func600 <- function(x){Process_SAS(delta[x],struc,eps,n0=600,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE[i],c[i],sm="SAS",Lc,ns)$TypeI}
  TypeI$Res_SAS_n0600[,i] <- sapply(1:length(delta),func600)

  tic2 <- proc.time()-tic
  print(paste("Sam",Model,"-n",n,"-c",c[i],"-time:",tic2[3],sep=""))

}

TypeI
TypeI_bandwidth <- TypeI
write.csv(TypeI,file='TypeI.csv')
write.csv(TypeI_bandwidth,file='TypeI_bandwidth.csv')


#######################################################################################
##### power comparison between SAS and US #############################################
#under IID and COR

eps <- "norm"       #norm  chisq  t3
sdr <- "MAVE"
ker <- "Epanechnikov"

n0 <-300
c <- 0.1
cE <- 0.5

delta <- seq(0,1,0.125)
Res_temp <- matrix(0,length(delta),length(c))
TypeI <- list(Res_US_IID = Res_temp, Res_SAS_IID = Res_temp,
              Res_US_COR = Res_temp, Res_SAS_COR = Res_temp,
              Res_US0_IID = Res_temp, Res_US0_COR = Res_temp)

tic <- proc.time()

#sm="US", IID
func0 <- function(x){Process_SAS(delta[x],struc="IID",eps,n0,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="US",Lc,ns)$TypeI}
TypeI$Res_US_IID <- sapply(1:length(delta),func0)

#sm="SAS", IID
func1 <- function(x){Process_SAS(delta[x],struc="IID",eps,n0,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="SAS",Lc,ns)$TypeI}
TypeI$Res_SAS_IID <- sapply(1:length(delta),func1)

#sm="US", COR
func00 <- function(x){Process_SAS(delta[x],struc="COR",eps,n0,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="US",Lc,ns)$TypeI}
TypeI$Res_US_COR <- sapply(1:length(delta),func00)

#sm="SAS", COR
func11 <- function(x){Process_SAS(delta[x],struc="COR",eps,n0,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="SAS",Lc,ns)$TypeI}
TypeI$Res_SAS_COR <- sapply(1:length(delta),func11)

#sm="US0", IID
func000 <- function(x){Process_SAS(delta[x],struc="IID",eps,n0,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="US0",Lc,ns)$TypeI}
TypeI$Res_US0_IID <- sapply(1:length(delta),func000)

#sm="US0", COR
func111 <- function(x){Process_SAS(delta[x],struc="COR",eps,n0,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="US0",Lc,ns)$TypeI}
TypeI$Res_US0_COR <- sapply(1:length(delta),func111)

tic2 <- proc.time()-tic
print(paste("Sam",Model,"-n",n,"-c",c,"-time:",tic2[3],sep=""))


TypeI
TypeI_struc <- TypeI
write.csv(TypeI,file='TypeI.csv')
write.csv(TypeI_struc,file='TypeI_struc.csv')

#######################################################################################
##### power of SAS for different SDR methods###########################################

eps <- "norm"
ker <- "Epanechnikov"
n0 <- 400
c <- 0.1
cE <- 0.5

delta <- seq(0,0.5,0.25)
Res_temp <- matrix(0,length(delta),length(c))
TypeI <- list(Res_SAS_IID_MAVE = Res_temp, Res_SAS_COR_MAVE = Res_temp,
              Res_SAS_IID_SAVE = Res_temp, Res_SAS_COR_SAVE = Res_temp,
              Res_SAS_IID_DR = Res_temp, Res_SAS_COR_DR = Res_temp)

tic <- proc.time()
#sm="SAS",struc="IID",sdr="MAVE"
func1_IIDMAVE <- function(x){Process_SAS(delta[x],struc="IID",eps,n0,n,N,p,d,sdr="MAVE",ker,Model,dm,replace,alpha,cE,c,sm="SAS",Lc,ns)$TypeI}
TypeI$Res_SAS_IID_MAVE <- sapply(1:length(delta),func1_IIDMAVE)

#sm="SAS",struc="COR",sdr="MAVE"
func1_CORMAVE <- function(x){Process_SAS(delta[x],struc="COR",eps,n0,n,N,p,d,sdr="MAVE",ker,Model,dm,replace,alpha,cE,c,sm="SAS",Lc,ns)$TypeI}
TypeI$Res_SAS_COR_MAVE <- sapply(1:length(delta),func1_CORMAVE)

#sm="SAS",struc="IID",sdr="SAVE"
func1_IIDSAVE <- function(x){Process_SAS(delta[x],struc="IID",eps,n0,n,N,p,d,sdr="SAVE",ker,Model,dm,replace,alpha,cE,c,sm="SAS",Lc,ns)$TypeI}
TypeI$Res_SAS_IID_SAVE <- sapply(1:length(delta),func1_IIDSAVE)

#sm="SAS",struc="COR",sdr="SAVE"
func1_CORSAVE <- function(x){Process_SAS(delta[x],struc="COR",eps,n0,n,N,p,d,sdr="SAVE",ker,Model,dm,replace,alpha,cE,c,sm="SAS",Lc,ns)$TypeI}
TypeI$Res_SAS_COR_SAVE <- sapply(1:length(delta),func1_CORSAVE)

#sm="SAS",struc="IID",sdr="DR"
func1_IIDDR <- function(x){Process_SAS(delta[x],struc="IID",eps,n0,n,N,p,d,sdr="DR",ker,Model,dm,replace,alpha,cE,c,sm="SAS",Lc,ns)$TypeI}
TypeI$Res_SAS_IID_DR <- sapply(1:length(delta),func1_IIDDR)

#sm="SAS",struc="COR",sdr="DR"
func1_CORDR <- function(x){Process_SAS(delta[x],struc="COR",eps,n0,n,N,p,d,sdr="DR",ker,Model,dm,replace,alpha,cE,c,sm="SAS",Lc,ns)$TypeI}
TypeI$Res_SAS_COR_DR <- sapply(1:length(delta),func1_CORDR)

tic2 <- proc.time()-tic
print(paste("Sam",Model,"-n",n,"-c",c,"-time:",tic2[3],sep=""))


TypeI
TypeI_sdr <- TypeI
write.csv(TypeI,file='TypeI.csv')
write.csv(TypeI_sdr,file='TypeI_sdr.csv')

#######################################################################################
##### power of SAS for different kernel function ######################################

eps <- "norm"
sdr <- "MAVE"
n0 <- 300
c <- 0.1
cE <- 0.5

delta <- seq(0,0.5,0.25)
Res_temp <- matrix(0,length(delta),length(c))
TypeI <- list(Res_SAS_IID_Epan = Res_temp, Res_SAS_COR_Epan = Res_temp,
              Res_SAS_IID_Quartic = Res_temp, Res_SAS_COR_Quartic = Res_temp,
              Res_SAS_IID_Triweight = Res_temp, Res_SAS_COR_Triweight = Res_temp,
              Res_SAS_IID_Triangular = Res_temp, Res_SAS_COR_Triangular = Res_temp)

tic <- proc.time()
#sm="SAS",struc="IID",ker = "Epanechnikov"
func1_IID_Epan <- function(x){Process_SAS(delta[x],struc="IID",eps,n0,n,N,p,d,sdr,ker="Epanechnikov",Model,dm,replace,alpha,cE,c,sm="SAS",Lc,ns)$TypeI}
TypeI$Res_SAS_IID_Epan <- sapply(1:length(delta),func1_IID_Epan)

#sm="SAS",struc="COR",ker = "Epanechnikov"
func1_COR_Epan <- function(x){Process_SAS(delta[x],struc="COR",eps,n0,n,N,p,d,sdr,ker="Epanechnikov",Model,dm,replace,alpha,cE,c,sm="SAS",Lc,ns)$TypeI}
TypeI$Res_SAS_COR_Epan <- sapply(1:length(delta),func1_COR_Epan)

#sm="SAS",struc="IID",ker = "Quartic"
func1_IID_Quartic <- function(x){Process_SAS(delta[x],struc="IID",eps,n0,n,N,p,d,sdr,ker="Quartic",Model,dm,replace,alpha,cE,c,sm="SAS",Lc,ns)$TypeI}
TypeI$Res_SAS_IID_Quartic <- sapply(1:length(delta),func1_IID_Quartic)

#sm="SAS",struc="COR",ker = "Quartic'
func1_COR_Quartic <- function(x){Process_SAS(delta[x],struc="COR",eps,n0,n,N,p,d,sdr,ker="Quartic",Model,dm,replace,alpha,cE,c,sm="SAS",Lc,ns)$TypeI}
TypeI$Res_SAS_COR_Quartic <- sapply(1:length(delta),func1_COR_Quartic)

#sm="SAS",struc="IID",ker = "Triweight"
func1_IID_Triweight <- function(x){Process_SAS(delta[x],struc="IID",eps,n0,n,N,p,d,sdr,ker="Triweight",Model,dm,replace,alpha,cE,c,sm="SAS",Lc,ns)$TypeI}
TypeI$Res_SAS_IID_Triweight <- sapply(1:length(delta),func1_IID_Triweight)

#sm="SAS",struc="COR",ker="Triweight"
func1_COR_Triweight <- function(x){Process_SAS(delta[x],struc="COR",eps,n0,n,N,p,d,sdr,ker="Triweight",Model,dm,replace,alpha,cE,c,sm="SAS",Lc,ns)$TypeI}
TypeI$Res_SAS_COR_Triweight <- sapply(1:length(delta),func1_COR_Triweight)

#sm="SAS",struc="IID",ker = "Triangular"
func1_IID_Triangular <- function(x){Process_SAS(delta[x],struc="IID",eps,n0,n,N,p,d,sdr,ker="Triangular",Model,dm,replace,alpha,cE,c,sm="SAS",Lc,ns)$TypeI}
TypeI$Res_SAS_IID_Triangular <- sapply(1:length(delta),func1_IID_Triangular)

#sm="SAS",struc="COR",ker="Triangular"
func1_COR_Triangular<- function(x){Process_SAS(delta[x],struc="COR",eps,n0,n,N,p,d,sdr,ker="Triangular",Model,dm,replace,alpha,cE,c,sm="SAS",Lc,ns)$TypeI}
TypeI$Res_SAS_COR_Triangular <- sapply(1:length(delta),func1_COR_Triangular)
tic2 <- proc.time()-tic

print(paste("Sam",Model,"-n",n,"-c",c,"-time:",tic2[3],sep=""))

TypeI
TypeI_ker <- TypeI
write.csv(TypeI,file='TypeI.csv')
write.csv(TypeI_ker,file='TypeI_ker.csv')



###########################################################################################
##### power comparison between SAS and US for different error terms ###################
#Under varying coefficient model

struc <- "IID"
Model <- "Varying"; p <- 4; dm <- 4

sdr <- "MAVE"
ker <- "Epanechnikov"

n <- 1000
n0 <- 300
c <- 0.1
cE <- 0.5

delta <- seq(0,1,0.125)
Res_temp <- matrix(0,length(delta),length(c))
TypeI <- list(Res_US_norm = Res_temp, Res_SAS_norm = Res_temp,
              Res_US_chisq = Res_temp, Res_SAS_chisq = Res_temp,
              Res_US_t3 = Res_temp, Res_SAS_t3 = Res_temp,
              Res_US0_norm = Res_temp,Res_US0_chisq = Res_temp,
              Res_US0_t3 = Res_temp)

tic <- proc.time()
#sm="US", eps="norm"
func0_norm <- function(x){Process_SAS(delta[x],struc,eps="norm",n0,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="US",Lc,ns)$TypeI}
TypeI$Res_US_norm <- sapply(1:length(delta),func0_norm)

#sm="SAS",eps="norm"
func1_norm <- function(x){Process_SAS(delta[x],struc,eps="norm",n0,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="SAS",Lc,ns)$TypeI}
TypeI$Res_SAS_norm <- sapply(1:length(delta),func1_norm)

#sm="US", eps="chisq"
func0_chisq <- function(x){Process_SAS(delta[x],struc,eps="chisq",n0,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="US",Lc,ns)$TypeI}
TypeI$Res_US_chisq <- sapply(1:length(delta),func0_chisq)

#sm="SAS",eps="chisq"
func1_chisq <- function(x){Process_SAS(delta[x],struc,eps="chisq",n0,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="SAS",Lc,ns)$TypeI}
TypeI$Res_SAS_chisq <- sapply(1:length(delta),func1_chisq)

#sm="US", eps="t3"
func0_t3 <- function(x){Process_SAS(delta[x],struc,eps="t3",n0,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="US",Lc,ns)$TypeI}
TypeI$Res_US_t3 <- sapply(1:length(delta),func0_t3)

#sm="SAS",eps="t3"
func1_t3 <- function(x){Process_SAS(delta[x],struc,eps="t3",n0,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="SAS",Lc,ns)$TypeI}
TypeI$Res_SAS_t3 <- sapply(1:length(delta),func1_t3)

#sm="US0", eps="norm"
func00_norm <- function(x){Process_SAS(delta[x],struc,eps="norm",n0,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="US0",Lc,ns)$TypeI}
TypeI$Res_US0_norm <- sapply(1:length(delta),func00_norm)

#sm="US0", eps="chisq"
func00_chisq <- function(x){Process_SAS(delta[x],struc,eps="chisq",n0,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="US0",Lc,ns)$TypeI}
TypeI$Res_US0_chisq <- sapply(1:length(delta),func00_chisq)

#sm="US", eps="t3"
func00_t3 <- function(x){Process_SAS(delta[x],struc,eps="t3",n0,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="US0",Lc,ns)$TypeI}
TypeI$Res_US0_t3 <- sapply(1:length(delta),func00_t3)

tic2 <- proc.time()-tic
print(paste("Sam",Model,"-n",n,"-c",c,"-time:",tic2[3],sep=""))


TypeI
TypeI_eps <- TypeI
write.csv(TypeI,file='TypeI.csv')
write.csv(TypeI_eps,file='TypeI_eps.csv')



#######################################################################################
#time comparison under varying coefficient model#######################################
#eps = norm, fix n

##simulation one time
struc <- "IID"
Model <- "Varying"; p <- 4; dm <- 4

eps <- "norm"
sdr <- "MAVE"
ker <- "Epanechnikov"

n <- 1000
n0 <- 300
c <- 0
cE <- 0.5
N <- c(10^3,10^3,5*(10^3),10^4,5*(10^4),10^5,5*(10^5),10^6,5*(10^6),10^7)
delta <- 0
seedNo <- 1
result_temp <- matrix(0,length(N),6)
rownames(result_temp) <- N
colnames(result_temp) <- c("tt","Tn","Lc","time_pilot","time_PI","time_stat")
result <- list(result_US0 = result_temp, result_US = result_temp,
               result_SAS = result_temp, result_FULL = result_temp)
time_temp <- matrix(0,length(N),length(n))
rownames(time_temp) <- N
time <- list(time_US0 = time_temp,time_US = time_temp,
             time_SAS = time_temp, time_FULL = time_temp)


tic <- proc.time()
for (i in 1:length(N)){

  set.seed(1234567)

  Sampleall <- Gendata(N[i], p, delta, struc, Model, eps)

  #sm="US0"
  result$result_US0[i,] <- SAS(Sampleall,n0,n,N[i],p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="US0",Lc,seedNo)
  time$time_US0[i,] <- sum((result$result_US0[i,])[4:6])

  #sm="US"
  result$result_US[i,] <- SAS(Sampleall,n0,n,N[i],p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="US",Lc,seedNo)
  time$time_US[i,] <- sum((result$result_US[i,])[4:6])

  #sm="SAS"
  result$result_SAS[i,] <- SAS(Sampleall,n0,n,N[i],p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="SAS",Lc,seedNo)
  time$time_SAS[i,] <- sum((result$result_SAS[i,])[4:6])

  #FULL, sm ="US"
  result$result_FULL[i,] <- SAS(Sampleall,n0,n=N[i],N[i],p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm="US",Lc,seedNo)
  time$time_FULL[i,] <- sum((result$result_FULL[i,])[4:6])

  tic2 <- proc.time()-tic
  print(paste("Sam",Model,"-N",N[i],"-c",c,"-time:",tic2[3],sep=""))

} 
time
result

write.csv(time,file='time.csv')
write.csv(result,file='result.csv')
