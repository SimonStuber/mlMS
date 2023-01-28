init_randomMLAR <- function(y, xOutcome,nTime, constants){
  resVar <- c()#
  ar <- c()
  int <- c()
  for(i in 1:ncol(y)){
    mod <- ar(y[,i],order.max = 1, aic=FALSE, method="ols")

    if(!length(mod$ar)==0){
      ar[i] <- mod$ar
      resVar[i] <- var(mod$resid, na.rm=TRUE)
      int[i] <- mean(y[,i], na.rm=TRUE)
    }else{
      ar[i] <- 0
      resVar[i] <- var(y[,i], na.rm=TRUE)
      int[i] <- mean(y[,i], na.rm=TRUE)
    }
  }

  # yLag <- rbind(NA,y[-nTime,])
  # for(i in 1:ncol(y)){
  #   yLag[,i] <- yLag[,i] - mean(yLag[,i], na.rm=TRUE)
  #
  # }
  # lmerShapeDat <- cbind(as.vector(y),as.vector(yLag), rep(1:ncol(y), each=nrow(y)))
  # colnames(lmerShapeDat) <- c("y", "yLag","id")
  #
  # LmerMod <- lmer(y ~ yLag + (yLag|id), data=as.data.frame(lmerShapeDat),
  #                 REML = FALSE,
  #                 control = lmerControl(optimizer ="bobyqa"))
  # #lm(dat$xAsOutcome ~ ranef(LmerMod)$id[,1] + ranef(LmerMod)$id[,2])
  # intVarMlEst <- var(int)
  # arVarMlEst <- var(ar)
  # resVarMlEst <- var(resVar)
##
  lm <- log((mean(resVar)^2)/sqrt((mean(resVar)^2)+(var(resVar))))
  lv <- log(1+((var(resVar))/(mean(resVar)^2)))

  b0 <- int
  b1 <- ar
  res <- resVar
  sds <- sqrt(c(var(int),var(ar),lv))

  eff <- cbind(b0, b1, log(res))
  covMat <- cov(eff)
  diag(covMat) <- sds^2
  U <- chol(covMat)

  inits <- list(effMeans=c(mean(int),mean(ar),lm),
                U=U,
                Ustar=U/sds,
                b0=b0,
                b1=b1,
                res=res,
                sds=sds,
                eff=eff)




  if(constants$predAr){
    inits$bArPred <- 0
  }
  if(constants$predResVar){
    inits$bResVarPred <- 0
  }
  if(constants$predMean){
    inits$bMeanPred <- 0
  }
  if(constants$predX){
    bb <- c()
    bb[1] <- mean(xOutcome)
    bb[2] <- 0
    bb[3] <- 0
    bb[4] <- 0
    inits$bb <- bb
    inits$xOutResVar <- 1/var(xOutcome)
  }
  return(inits)
}
