
init_mlHMM <- function(simdat,M, forwardAlgorithm=TRUE, particleFilter=FALSE){
#browser()
#simdat <- round(simdat)
nVars <- dim(simdat)[3]
n <- dim(simdat)[2]
nTime <- dim(simdat)[1]
init_b.prec <- array(0, c(nVars^2,nVars^2,M))
resCovInit <- array(0, c(n,nVars, nVars, M))
bRInits <- array(0, c(n, nVars^2,M))
b1Inits <- array(0, c(n, nVars,M))
betaMat <- array(0, c(n, nVars^2))
mm <- array(NA, c(nTime,n))
mlVarMod <- list(list())
mlVarDat <- array(simdat, c(n*nTime,nVars))
mlVarDat <- as.data.frame(mlVarDat)
mlVarDat$ID <- rep(1:n, each=nTime)
colnames(mlVarDat) <- c(paste("VAR", 1:nVars, sep=""), "ID")
km.res <- list()

#browser()
km.all <- kmeans(mlVarDat[,1:nVars],centers = M,nstart = nTime,iter.max = 50)

for(i in 1:n){

  km.res[[i]] <- try(kmeans(mlVarDat[mlVarDat$"ID"==i,1:nVars],centers = km.all$centers,nstart = nTime,
                   iter.max = 50),silent = TRUE)
  if(inherits(km.res[[i]], "try-error")){
    km.res[[i]] <- km.all
  }
  # if(any(km.res[[i]]$size < round(nTime/100))){
  #   km.res[[i]] <- km.all
  # }

#  mm[,i] <- as.vector(km.res[[i]]$cluster)

  for(state in 1:M){
    #betaMat <- array(NA, c(n, nVars^2))

     newDat <- mlVarDat[mlVarDat$"ID"==i&km.res[[i]]$cluster==state,]
    #
    # mlVarMod<- try(graphicalVAR(data = newDat[,1:nVars],verbose = FALSE,
    #                         lambda_kappa = 0,lambda_beta = 0,
    #                         centerWithin = TRUE,scale=FALSE),silent = TRUE)
      if(nrow(newDat)<(nVars+1)){
        resCovInit[i,,,state] <- cov(newDat[,1:nVars])
        bRInits[i,,state] <- 0

        b1Inits[i,,state] <- km.res[[i]]$centers[state,]
      }else{

      resCovInit[i,,,state] <- cov(newDat[,1:nVars])

      b1Inits[i,,state] <- km.res[[i]]$centers[state,]
      }
  }
}


  a <- array(NA, c(n,2, nTime))
  vInit <- array(0, c(n,M,M))

  #km.res[[i]]$cluster

  for(i in 1:n){
    #browser()
    if(km.res[[i]]$tot.withinss==km.all$tot.withinss){
      clusters <- km.res[[i]]$cluster[(n*nTime-nTime+1):(n*nTime)]
      a[i,1,] <- c(NA, clusters[-nTime])
      a[i,2,] <- clusters
    }else{
      a[i,1,] <- c(NA, km.res[[i]]$cluster[-nTime])
      a[i,2,] <- km.res[[i]]$cluster
    }

    for(t in 2:(nTime)){
      vInit[rbind(c(i,a[i,,t]))] <- vInit[rbind(c(i,a[i,,t]))] + 1
    }

    if(sum(vInit[i,,]<(rowSums(vInit[i,,])/100))>0){
   #   browser()
      vInit[i,,][vInit[i,,]<(rowSums(vInit[i,,])/100)] <- matrix(rowSums(vInit[i,,]),M,M)[vInit[i,,]<(rowSums(vInit[i,,])/100)]
    }

    vInit[i,,] <-(vInit[i,,]/rowSums(vInit[i,,]))
    #diag(vInit[i,,]) <- 1/M
    #vInit[i,,] <- vInit[i,,]+(1/M)

    vInit[i,,] <- log(vInit[i,,]/diag(vInit[i,,]))

    diag(vInit[i,,]) <- NA
    vInit[i,,] <- t(vInit[i,,])
  }


  vMean <- apply(vInit, 2:3, mean)
  vCov <- apply(vInit, 2:3, var)

  b1.precInits <- apply(b1Inits, 2:3, sd)^2
  b1.precInits[b1.precInits==0] <- var(simdat)/nVars

  scaleResCov <- apply(resCovInit, 2:4, mean)
  for(state in 1:M){
    scaleResCov[,,state] <- scaleResCov[,,state]*((km.all$size/(n*nTime))*(n+nVars)*M)[state]
  }



  inits <- list(
    probs1 = km.all$size/(n*nTime),
    b1.means = t(km.all$centers),
   # bRaneff=bRInits,
    b1=b1Inits,
    b1.prec =b1.precInits,
  #  b.prec =init_b.prec,
   # b.means = t(apply(bRInits,2,colMeans)),
    v.means = vMean[!is.na(vMean)],
    v.prec = as.matrix(nearPD(diag(vCov[!is.na(vCov)]))$mat),
    scaleResCov = scaleResCov,
    dfResCov=ifelse(((km.all$size/(n*nTime))*(n+nVars)*M)>(nVars+4),
                    ((km.all$size/(n*nTime))*(n+nVars)*M),nVars+4),
    resCov = resCovInit,
   raneffs=array(vInit[!is.na(vInit)],c(n,M*M-M)))

  if(!forwardAlgorithm&!particleFilter){
    inits$mm <- mm
  }
  return(inits)
}



