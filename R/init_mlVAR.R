
init_mlVAR <- function(simdat){

nVars <- dim(simdat)[3]
n <- dim(simdat)[2]
nTime <- dim(simdat)[1]
init_b.prec <- array(0, c(nVars^2,nVars^2))
resCovInit <- array(0, c(n,nVars, nVars))
bRInits <- array(0, c(n, nVars*nVars))
b1Inits <- array(0, c(n, nVars))
betaMat <- array(0, c(n, nVars^2))
mm <- array(NA, c(nTime,n))
mlVarMod <- list(list())
mlVarDat <- array(simdat, c(n*nTime,nVars))
mlVarDat <- as.data.frame(mlVarDat)
mlVarDat$ID <- rep(1:n, each=nTime)
colnames(mlVarDat) <- c(paste("VAR", 1:nVars, sep=""), "ID")
km.res <- list()

#browser()

mlVarMod<- try(mlVAR(mlVarDat, idvar="ID", vars = colnames(mlVarDat)[-(nVars+1)], scale=FALSE, scaleWithin = FALSE),silent = TRUE)
for(i in 1:n){



  #mm[,i] <- as.vector(km.res[[i]]$cluster)





    if(inherits(mlVarMod, "try-error")){


      resCovInit[i,,] <- cov(simdat[,i,])

      arMat <- matrix(0,nVars, nVars)
      for(var in 1:nVars){
        arvar <- ar(newDat[,var], order=1)$ar
        if(length(arvar>0)){
          arMat[var,var] <- arvar
        }else{
          arMat[var,var] <- 0
        }

      }

      bRInits[i,] <- as.vector(arMat)


      b1Inits[i,] <- colMeans(simdat[,i,])

    }else{


    resCovInit[i,,] <- mlVarMod$results$Theta$cov$subject[[i]]

    bRInits[i,] <- as.vector(mlVarMod$results$Beta$subject[[i]])

    b1Inits[i,] <- colMeans(simdat[,i,])

    }

}
 # for(state in 1:M){


    init_b.prec <- diag(nVars*nVars)
    #diag(init_b.prec) <- as.vector((mlVarMod$results$Beta$SD)^2)
    #init_b.prec <- as.matrix(nearPD(init_b.prec)$mat)
  #}


  b1.precInits <- apply(b1Inits,2,var)
  b1.precInits[b1.precInits==0] <- var(simdat)/nVars

  scaleResCov <- apply(resCovInit,2:3,mean)*n
  #for(state in 1:M){
    #scaleResCov <- scaleResCov*((km.all$size/(n*nTime))*(n+nVars)*M)
  #}



  inits <- list(
   # probs1 = km.all$size/(n*nTime),
    b1.means = apply(simdat,3, mean),
    bRaneff=bRInits,
    b1=b1Inits,
    b1.prec =b1.precInits,
    b.prec =init_b.prec,
    b.means = colMeans(bRInits),
  #  v.means = vMean[!is.na(vMean)],
   # v.prec = as.matrix(nearPD(diag(vCov[!is.na(vCov)]))$mat),
    scaleResCov = scaleResCov,
    dfResCov=n,
    resCov = resCovInit)
  return(inits)
}



