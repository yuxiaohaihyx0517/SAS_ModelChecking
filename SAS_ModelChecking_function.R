#######################################################################################
#generate null and alternative model
Gendata <- function(n, p, delta, struc, Model, eps){
  
  #two covariance structure
  if (struc == "IID") {Sigma <- diag(rep(1,p))}
  if (struc == "COR") {
    Sigma <- matrix(rep(0,p^2),p,p)
    for (i in 1:p){
      for(j in 1:p){
        Sigma[i,j] <- 0.5^(abs(i-j))
      }
    }
  }
  
  #generate covarites
  X <- mvrnorm(n,rep(0,p),Sigma)
  
  #gernerate error term
  if (eps =="norm")  {ep <- rnorm(n,0,1)/sqrt(1)}
  if (eps =="chisq") {ep <- (rchisq(n,1,0)-1)/sqrt(2)}
  if (eps =="t3")    {ep <- rt(n,3,0)/sqrt(3)}
  
  
  #model: Y=G(X,beta,g)+ep
  if (Model == "LM"){
    Y <- X%*%(c(rep(1,p/2),rep(-1,p/2))/sqrt(p)) + delta*0.4*abs(0+X%*%(c(1,-1,rep(0,p-2))/sqrt(2)))^3 + ep
    Y0 <- X%*%(c(rep(1,p/2),rep(-1,p/2))/sqrt(p)) + ep
  }
  if (Model == "SIM"){
    Y <- 2*exp(-1*((X[,1]-X[,2])/sqrt(2))) + delta*0.4*(0+(X[,3]-X[,4])/sqrt(2))^2 + ep  
    Y0 <- 2*exp(-1*((X[,1]-X[,2])/sqrt(2))) + ep 
  }
  if (Model == "MIM"){
    Y <- (X%*%c(1,0,0,0)/sqrt(1))^2 + (X%*%c(0,1,0,0)/sqrt(1))^2 + delta*0.8*exp(-0.4*(X%*%(c(0,0,1,0)/sqrt(1)))) + 1*ep  
    Y0 <- (X%*%c(1,0,0,0)/sqrt(1))^2 + (X%*%c(0,1,0,0)/sqrt(1))^2 + 1*ep  
  }
  if (Model == "Varying"){
    Y <- ((X[,2])/sqrt(1))*(sin((X[,1])/sqrt(1))+cos((X[,1])/sqrt(1)))+2*(X[,1]*(1-X[,1]))*X[,3]+ delta*0.3*((X[,2]+X[,3]+X[,4])/sqrt(3))^3 + ep
    Y0 <- ((X[,2])/sqrt(1))*(sin((X[,1])/sqrt(1))+cos((X[,1])/sqrt(1)))+2*(X[,1]*(1-X[,1]))*X[,3] + ep
  }  
  
  return(list(X = X, Y = Y, Y0 = Y0, ep = ep, Sigma = Sigma)) 
}
 
#######################################################################################
#estimate the dimension reduction subspace B######

MAVE <- function(d, X, Y){
  A <- NULL
  SDR_MAVE <- mave(Y~X,method ='meanOPG',max.dim = 10) #CSOPG,meanOPG,CSMAVE,meanMAVE
  A <- coef(SDR_MAVE,d)
  return(A)
}
 
SAVE <- function(d, X, Y){
  A <- NULL
  SDR_SAVE <- dr(Y~X,data = as.data.frame(cbind(X,Y)),method ='save',nslice = 4)   
  A <- SDR_SAVE$evectors[,1:d]                 
  return(A)
}


matpower <- function(a,al){
  a <- round((a+t(a))/2,7)
  tem <- eigen(a)
  return(tem$vectors%*%diag((tem$values)^al)%*%t(tem$vectors))
}
discretize <- function(Y, h){
  n <- length(Y)
  m <- floor(n/h)
  Y <- Y+.00001*mean(Y)*rnorm(n)
  yord <- Y[order(Y)]
  divpt <- numeric();for(i in 1:(h-1)) divpt <- c(divpt,yord[i*m+1])
  Y1 <- rep(0,n)
  Y1[Y<divpt[1]] <- 1
  Y1[Y>=divpt[h-1]] <- h
  for(i in 2:(h-1)) Y1[(Y>=divpt[i-1])&(Y<divpt[i])] <- i
  return(Y1)
}
DR <- function(d, X, Y){
  p <- ncol(X)
  n <- nrow(X)
  h <- 10
  signrt <- matpower(var(X),-1/2)
  xc <- t(t(X)-apply(X,2,mean))
  xst <- xc%*%signrt
  #if(ytype=="continuous") 
  ydis <- discretize(Y,h)
  # if(ytype=="categorical") ydis=y
  yless <- ydis
  ylabel <- numeric()
  for(i in 1:n) {
    if(var(yless)!=0) {ylabel <- c(ylabel,yless[1])
    yless <- yless[yless!=yless[1]]}}
  ylabel <- c(ylabel,yless[1])
  prob <- numeric()
  for(i in 1:h) prob <- c(prob,length(ydis[ydis==ylabel[i]])/n)
  vxy <-  array(0,c(p,p,h))
  exy <- numeric()
  for(i in 1:h) {
    vxy[,,i] <- var(xst[ydis==ylabel[i],])
    exy <- rbind(exy,apply(xst[ydis==ylabel[i],],2,mean))
  }
  mat1 <-  matrix(0,p,p)
  mat2 <- matrix(0,p,p)
  for(i in 1:h){
    mat1 <-  mat1+prob[i]*(vxy[,,i]+exy[i,]%*%t(exy[i,]))%*%
      (vxy[,,i]+exy[i,]%*%t(exy[i,]))
    mat2 <- mat2+prob[i]*exy[i,]%*%t(exy[i,])}
  out <- 2*mat1+2*mat2%*%mat2+2*sum(diag(mat2))*mat2-2*diag(p)
  return(signrt%*%eigen(out)$vectors[,1:d])
}

SDR <- function(d, X, Y, sdr){
  
  if (sdr == "SIR")  {B <- SIR(d, X, Y)}
  if (sdr == "SAVE") {B <- SAVE(d, X, Y)}
  if (sdr == "MAVE") {B <- MAVE(d,X,Y)}
  if (sdr == "DR")   {B <- DR(d, X, Y)}
  
  return(B)
}


#######################################################################################
#estimate semi-parametric model to get fitted residual

Resd <- function(X, Y, sdr, Model, dm){
 
  if (Model == "LM"){
    beta <- solve(t(X)%*%X)%*%t(X)%*%Y
    Yhat <- X%*%beta
    eps_hat <- Y - Yhat
  }
  
  if (Model == "SIM"){
    beta <- SDR(dm, X, Y, "MAVE")
    data_train <- data.frame(Y = Y, dir1 = X%*%beta)
    gam <- gam(Y~s(dir1), data = data_train)
    eps_hat <- gam$residuals
  }
  
  if (Model =="MIM"){
    beta <- SDR(dm, X, Y,"MAVE")
    data_train <- data.frame(dir1 = (X%*%beta)[,1],dir2 =(X%*%beta)[,2] ,Y = Y)
    gam <- gam(Y~s(dir1)+s(dir2),data = data_train)
    eps_hat <- gam$residuals
  }
  
  if (Model == "Varying"){
    beta <- diag(1,dm)
    data_train <- data.frame(Y = Y, u = X[,1], Z1 = X[,2], Z2 = X[,3])
    gam <- gam(Y~s(u,by=Z1)+s(u,by=Z2),data = data_train)
    eps_hat <- gam$residuals
  }
  
  return(list(eps_hat = eps_hat, beta = beta))
}

#######################################################################################
#distance matrix for all paired reduced sample with reduced subspace 

Dist <- function(x, y, B){
  dist <- dist2(x%*%B, y%*%B, method ="euclidean")
  return(dist)
}

#######################################################################################
#kernel function

kernel <- function(x, y, B, id1, id2, h, d, ker){
 
  dist <- dist2(x[id1,]%*%B,y[id2,]%*%B, method="euclidean") #nrow(X)*nrow(Y)
  dmat <- dist
  
  if (ker == "Epanechnikov")   {KI <- 0.75*(1-(dmat/h)^2)/(h^d)}    #Epanechnikov
  if (ker == "Quartic")        {KI <- 15/16*(1-(dmat/h)^2)^2/(h^d)}      #Quartic
  if (ker == "Triweight")      {KI <- 35/32*(1-(dmat/h)^2)^3/(h^d)}    #Triweight
  if (ker == "Triangular")     {KI <- (1-abs(dmat/h))/(h^d)}            #triangular
 
   KI[KI<0] <- 0                 #remove the case when u>1
  
   return(KI)
  
}

#######################################################################################
#normalized test statistic Tn

Stat <- function(X, B, N, n, p, d, h, idx, fx, resd , ker){
  
  K <- kernel(X, X, B, idx, idx, h, d, ker)     #Extract Submatrix by idx
  diag(K) <- 0               #remove the diag
  pie <- 1/(fx[idx]*N)     #the probability of each data point BX
 
  #compute the test statistics
  temp <- resd*sqrt(pie)
  Vn <- t(temp)%*%K%*%temp/(n*(n-1))
  Vn <- n*h^(d/2)*Vn*sqrt((n-1)/n)
  
  resdH0 <- resd - mean(resd)
  tempV <- resdH0*sqrt(pie)
  var <- t(tempV^2)%*%(K^2)%*%(tempV^2)*(h^d)*2/(n*(n-1))  #covarinace Vn
  sd <- sqrt(var)
  
  Tn <- Vn/sd
 
  return(as.numeric(Tn))
}

#######################################################################################
####optimal sampling distribution
fx_optimal <- function(X, B, id1, y1, id2, hE, d, xi,ker){
  #compute the kernel distance as a weight:length(x)*length(x1)
  weight <- kernel(X, X, B, id1, id2, hE, d, ker)  #n*N
  
  #compute the sampling density at each reduced sample location of t(B)%*%x1
  Eomega <- sapply(id2,function(a){sum(weight[,a]*y1)/sum(weight[,a])}) #length(x) is the sample size of x
  Eomega[is.na(Eomega)] <- 0    #set 0 for missing 
  
  gs <- Eomega^2
  gs <- sapply(id2,function(a){max(gs[a],xi)})  #N*1
  fx <- gs/sum(gs)
  
  return(fx)
} 

#true sample density for original data point
fx_True <- function(N, resd, fx0, xi){
  
  gs <- (resd)^2
  gs <- sapply(1:N,function(a) {max(gs[a],xi)})
  fxTrue <- gs/sum(gs)
 
  if(sum(resd)^2==0) {fxTrue <- fx0}
 
 return(fxTrue)
}

#uniform sampling
fx_0 <- function(N){
  fx0 <- rep(1/N,N)
  
  return(fx0)
} 

#######################################################################################
#optimal sample test

SAS <- function(Sampleall,n0,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm,Lc,seedNo){
  
  set.seed(123456+10000*seedNo)
  
  #get the original data
  Xall <- Sampleall$X
  Yall <- Sampleall$Y
  YallH0 <- Sampleall$Y0
  epall <- Sampleall$ep
  
  ###pilot study###############################################
  #get the initial sampling index, subsample in stage 1 from Unifom sampling
  tic_pilot <- proc.time()
  idx0 <- sample(1:N, size = n0, replace = replace)
  X0 <- Xall[idx0,]
  Y0 <- Yall[idx0]        #exact subsample,use to estimate projection direction and sampling density
   
 
  #compute the residual in pilot study,estimate beta and get residual
  resd0 <- Resd(X0, Y0, sdr, Model, dm)$eps_hat
   
  #estimate the dimension reduction subsapce
  B0 <- SDR(d, X0, resd0, sdr)                 #estimate dimension reduction subspace, (X,residual)
  
  #bandwidth
  hE0 <- cE*n0^(-1/(4+d))*sum(apply(X0%*%B0,2,sd)) #bandwidth for estimation in pilot
  h0 <- hE0^(2+c)  
  time_pilot <- (proc.time()-tic_pilot)[3]
  
  
  #threshold
  xi <- 0.01/(N*h0^d)   
  #xi <- 0.01/N
  
  ##test stage############################################
  tt <- 0   #reject-1, accept-0
  Tn <- 0    #statistics in optimal sampling test
  
  tic_PI <- proc.time()
  #select a sampling scheme, sampling method
  if (sm == "US0"){
    fx <- fx_0(N)
    idx <- sample(1:N, size = n, replace = replace, prob=fx)
  } 
  if (sm == "US"){
    fx <- fx_0(N)
    idx <- sample(1:N, size = n, replace = replace, prob=fx)
    idx <- c(idx0,idx)
    n <- length(idx)
  }
  if (sm == "SAS"){
    if(n < N){
      fx <- fx_optimal(Xall, B0, idx0, resd0, 1:N, hE0, d, xi,ker)
      idx <- sample(1:N,size = n,replace = replace, prob = fx)
      #idx=c(idx0,idx)
      #n=length(idx)
    }
    if (n == N){
      fx <- fx_0(N)
      idx <- 1:N
    }
  }  
  #sampling according to the estimated l(t(B)%*%X) density
  if (sm == "OS"){
    if (Model == "LM")      {B0 <- c(1,-1,0,0)/sqrt(2) }              #true dimension
    if (Model == "SIM")     {B0 <- c(0,0,1,-1,0,0)/sqrt(2)}
    if (Model == "MIM")     {B0 <- c(0,0,1,0)/sqrt(1)}
    if (Model == "Varying") {B0 <- c(0,1,1,1)/sqrt(3)}  #no dimension structure in varying coefficent model
    
    fx0 <- fx_0(N)
    fxTrue <- fx_True(N,Yall-YallH0,fx0,xi)   
    fx <- fxTrue
    idx <- sample(1:N,size = n,replace = replace, prob = fx)
  }     
  time_PI <- (proc.time()-tic_PI)[3]
  
  #sampling
  tic_stat <- proc.time()
  X <- Xall[idx,]
  Y <- Yall[idx]
  YH0 <- YallH0[idx]
  
  resd <- Resd(X,Y,sdr,Model,dm)$eps_hat
  B <- B0     #consider the estimate of B in pilot study is credible
  
  hE <- cE*n^(-1/(4+d))*sum(apply(X%*%B,2,sd))
  h <- hE^(2+c)
  
  #compute the test statistic
  Tn <- Stat(Xall, B, N, n, p, d ,h , idx, fx, resd,ker)
  
  time_stat <- (proc.time()-tic_stat)[3]
  
  if(Tn>=Lc){tt=1}              #one-side
   return(c(tt,Tn,Lc,time_pilot,time_PI,time_stat))
}

Process_SAS <- function(delta,struc,eps,n0,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm,Lc,ns){
  
  source("SAS_ModelChecking_function.R",local=TRUE)
  
  set.seed(1234567)
  
  Sampleall <- Gendata(N, p, delta, struc, Model, eps)
  Xall <- Sampleall$X
  Yall <- Sampleall$Y
  YallH0 <- Sampleall$Y0
  epall <- Sampleall$ep
  
  #model checking test based on optimal sampling
  result <- foreach(seedNo=1:ns,.combine='rbind',.packages=c("dr","flexclust","MASS","MAVE","mgcv"))%dopar%
    SAS(Sampleall,n0,n,N,p,d,sdr,ker,Model,dm,replace,alpha,cE,c,sm,Lc,seedNo)
  
  TypeI <- mean(result[,1])
  Tn <- result[,2]
  Lc <- result[,3]
  return(list(TypeI=TypeI,result=result))
  
} 
