init_randomMLAR <- function(y, xOutcome,nTime, constants){
  resVar <- c()#
  ar <- c()
  int <- c()
  for(i in 1:ncol(y)){
    mod <- ar(y[,i],order.max = 1)
    resVar[i] <- var(mod$resid, na.rm=TRUE)
    if(!length(mod$ar)==0){
      ar[i] <- mod$ar
    }else{
      ar[i] <- 0
    }

    int[i] <- mod$x.mean
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
#
  lm <- log((mean(resVar)^2)/sqrt((mean(resVar)^2)+(var(resVar))))
  lv <- log(1+((var(resVar))/(mean(resVar)^2)))
  U <- chol(diag(c(var(int),var(ar),lv)))
  inits <- list(effMeans=c(mean(int),mean(ar),lm),
                U=U,
                Ustar=U/sqrt(c(var(int),var(ar),lv)))

  inits$b0 <- int
  inits$b1 <- ar
  inits$res <- resVar
  inits$sds <- sqrt(c(var(int),var(ar),lv))

  eff <- cbind(inits$b0, inits$b1, exp(inits$res))

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
    inits$bb[4] <- 0
    inits$bb[2] <- 0
    inits$bb[3] <- 0
    inits$bb[1] <- mean(xOutcome)
    inits$xOutResVar <- var(xOutcome)
  }
  return(inits)
}
