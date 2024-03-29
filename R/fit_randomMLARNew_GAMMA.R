library(nimble)
library(lme4)
library(BBmisc)




fit_randomMLAR <- function(y, niter=30000, nburnin=20000,
                             thin=1, nchains=1, inits=NULL,
                             summary=TRUE,
                             arOnX=NULL,
                             resVarOnX=NULL,
                             meanOnX=NULL,
                             xOutcome=NULL,
                             xOn=c("mean","ar","resVar"),
                             randomRes=TRUE,
                             estimateAffectDynamics=TRUE,
                             buildOnly=FALSE,
                             xOnFullAD=FALSE){

  #browser()
  # uppertri_mult_diag <- nimbleFunction(
  #   run = function(mat = double(2), vec = double(1)) {
  #     returnType(double(2))
  #     p <- length(vec)
  #     out <- matrix(nrow = p, ncol = p, init = FALSE)
  #     for(i in 1:p)
  #       out[ , i] <- mat[ , i] * vec[i]
  #     return(out)
  #   })

  nTime <- dim(y)[1]
  N <- dim(y)[2]
  constants <- list(N=N,
                    nTime=nTime)
  if(xOnFullAD){
    constants$fullAD <- 1
  }else{
    constants$fullAD <- 0
  }

  if(estimateAffectDynamics){
    constants$estimateAffectDynamics <- 1
  }else{
    constants$estimateAffectDynamics <- 0
  }

  if(!is.null(arOnX)){
    constants$predAr <- 1
  }else{
    constants$predAr <- 0
  }
  if(!is.null(resVarOnX)){
    constants$predResVar <- 1
  }else{
    constants$predResVar <- 0
  }
  if(!is.null(meanOnX)){
    constants$predMean <- 1
  }else{
    constants$predMean <- 0
  }

  if(!is.null(xOutcome)){
    constants$predX <- 1
  }else{
    constants$predX <- 0
  }

  constants$randomRes <- 1*randomRes


  lastT <- apply(y, 2,function(x)which.last(!is.na(x)))
  firstT <- apply(y, 2,function(x)which.first(!is.na(x)))

  constants$firstT <- firstT
  constants$lastT <- lastT

  if(!constants$randomRes){
    print("No random residual variances specified.
        Possible predictors for random residual variances are beeing ignored.
        Affect dynamics are calculated based on a fixed residual variance.")
  }

  modelBaseline <- nimbleCode({
    if(randomRes){
      for(i in 1:N){
        yhat[firstT[i],i] <- b0[i]
        resStart[i] <- res[i]/(1-(b1[i]^2))
        y[firstT[i], i] ~ dnorm(yhat[firstT[i],i], var=resStart[i])
        for(t in (firstT[i]+1):lastT[i]){
          yhat[t,i] <- b0[i] + b1[i]*(y[t-1,i] - b0[i])
          y[t, i] ~ dnorm(yhat[t,i], var=res[i])
        }
      }

    }else{
      for(i in 1:N){
        yhat[firstT[i],i] <- b0Start[i]
        resStart[i] <- res/(1-(b1[i]^2))
        y[firstT[i], i] ~ dnorm(yhat[firstT[i],i], var=resStart[i])
        for(t in (firstT[i]+1):lastT[i]){
          yhat[t,i] <- b0[i] + b1[i]*(y[t-1,i] - b0[i])
          y[t, i] ~ dnorm(yhat[t,i], var=res)
        }
      }
    }


    for(i in 1:N){
      if(predAr){
        b1[i] <- eff[i, 2] + bArPred*arOnX_c[i]
      }else{
        b1[i] <- eff[i, 2]
      }

      if(predMean){
        b0[i] <- eff[i, 1] + bMeanPred*meanOnX_c[i]
      }else{
        b0[i] <- eff[i, 1]
      }


      if(randomRes){
        if(predResVar){
          # res[i] <- exp( eff[i,3]+ bResVarPred*resVarOnX_c[i])
        }else{
          #    res[i] <- exp(eff[i,3])
        }
        eff[i,1:2] ~ dmnorm(effMeans[1:2], cholesky = U[1:2, 1:2], prec_param = 0)
      }else{
        eff[i,1:2] ~ dmnorm(effMeans[1:2], effPrec[1:2,1:2])
      }
    }

    if(predX){
      if(fullAD){
        for(i in 1:5){
          bb[i] ~ dnorm(0,.001)
        }
      }else{
        for(i in 1:4){
          bb[i] ~ dnorm(0,.001)
        }
      }
      xOutResVar ~ dgamma(.5,.5)

      if(fullAD){
        for(i in 1:N){
          xOutHat[i] <-  bb[1] + bb[2]*b0[i] + bb[3]*b1[i] + bb[4]*mssd[i] +bb[5]*var[i]
          xOutcome[i] ~ dnorm(xOutHat[i], xOutResVar)
        }

      }else{
        for(i in 1:N){
          xOutHat[i] <-  bb[1] + bb[2]*(b0[i]-mean(b0[1:N])) + bb[3]*(b1[i]-mean(b1[1:N])) + bb[4]*(res[i]-mean(res[1:N]))
          xOutcome[i] ~ dnorm(xOutHat[i], xOutResVar)
        }
      }

    }




    if(!randomRes){
      res <- exp(fixRes)
    }

    if(predAr){
      bArPred ~ dnorm(0,0.001)
      arOnX_c[1:N] <- (arOnX[1:N] -
                         mean(arOnX[1:N]))

    }

    if(predMean){
      bMeanPred ~ dnorm(0,0.001)
      meanOnX_c[1:N] <- (meanOnX[1:N] -
                           mean(meanOnX[1:N]))

    }

    if(predResVar){
      if(randomRes){
        bResVarPred ~ dnorm(0,0.001)
        resVarOnX_c[1:N] <- (resVarOnX[1:N] -
                               mean(resVarOnX[1:N]))
      }
    }



    effMeans[1] ~ dnorm(0, 0.001)
    effMeans[2] ~ dunif(-1,1)

    if(randomRes){
      # effMeans[3] ~ dnorm(0, 0.001)
      #   effPrec[1:3,1:3] ~ dwish(effPrecPriorMat[1:3, 1:3], 3)
      # Ustar[1:2,1:2] ~ dlkj_corr_cholesky(1, 2)
      # for(nsd in 1:3){
     # sds[1] ~ dunif(0,100)
    #  sds[2] ~ dunif(0,2)
      #sds[3] ~ dunif(0,100)
      #}

      if(randomRes){
      #  effMeans[3] ~ dnorm(0, 0.001)
        #   effPrec[1:3,1:3] ~ dwish(effPrecPriorMat[1:3, 1:3], 3)
        Ustar[1:3,1:3] ~ dlkj_corr_cholesky(1.5, 3)
        # for(nsd in 1:3){
        sds[1] ~ dunif(0,100)
        sds[2] ~ dunif(0,2)
        res.var ~ dunif(0,100)
        sds[3] <- sqrt(res.var)
        #}
        U[1:3,1:3] <- uppertri_mult_diag(Ustar[1:3, 1:3], sds[1:3])
        effVar[1:3,1:3] <- t(U[1:3,1:3])%*%(U[1:3,1:3])
}
      # effC ~ dnorm(0,.001)

      # effVar[1,1] <- sds[1]^2
      # effVar[2, 2] <- sds[2]^2
      # effVar[1, 2] <- effC
      # effVar[2, 1] <- effC

      phi <- (res.mean^2)/(res.var)
      mu <- res.mean
      shape <- phi+2
      scale <- mu*(1+phi)

      res.mean ~ dunif(0,100)
   #   res.var ~ dunif(0,100)
      for(i in 1:N){
        res[i] ~ dinvgamma(shape=shape, scale=scale)
        #  res[i] <- 1/resPrec[i]
        #var[i] <- (1/res[i])/(1-(b1[i]^2))
        #mssd[i] <- 2*var[i]*(1-b1[i])

      }
      #effVar[1:2,1:2] <- t(U[1:2,1:2])%*%(U[1:2,1:2])
    }else{
      fixRes ~ dnorm(0, 0.0000001)
      effPrec[1:2,1:2] ~ dwish(effPrecPriorMat[1:2, 1:2], 2)
      effVar[1:2,1:2]  <- inverse(effPrec[1:2,1:2])
    }


    if(estimateAffectDynamics){
      for(i in 1:N){
        if(randomRes){
          var[i] <- (res[i])/(1-(b1[i]^2))
          mssd[i] <- 2*var[i]*(1-b1[i])
        }else{
          var[i] <- (res)/(1-(b1[i]^2))
          mssd[i] <- 2*var[i]*(1-b1[i])
        }
      }
    }

  })

  yLag <- rbind(NA,y[-nTime,])

  for(i in 1:ncol(y)){
    yLag[,i] <- yLag[,i] - mean(yLag[,i], na.rm=TRUE)

  }


  dataList <- list("y"=y)

  if(constants$predAr){
    dataList$arOnX <- arOnX
  }

  if(constants$predMean){
    dataList$meanOnX <- meanOnX
  }

  if(constants$predResVar&constants$randomRes){
    dataList$resVarOnX <- resVarOnX
  }
  #
  if(constants$predX){
    dataList$xOutcome <- xOutcome
  }
  if(buildOnly){
    if(constants$randomRes){
      #  dataList$effPrecPriorMat <- diag(3)
      #  initEffMeans <- c(fixef(LmerMod)[1],fixef(LmerMod)[2], log(var(resid(LmerMod))))
    }else{
      dataList$effPrecPriorMat <- diag(2)
      #   initEffMeans <- c(fixef(LmerMod)[1],fixef(LmerMod)[2])
    }
  }else{
    if(constants$randomRes){
      #  dataList$effPrecPriorMat <- diag(3)
      #    initEffMeans <- c(fixef(LmerMod)[1],fixef(LmerMod)[2], log(var(resid(LmerMod))))
    }else{
      dataList$effPrecPriorMat <- diag(2)
      #     initEffMeans <- c(fixef(LmerMod)[1],fixef(LmerMod)[2])
    }
  }



  if(is.null(inits)){
    inits <- init_randomMLAR_G(y, xOutcome, nTime, constants)
  }
  #browser()
  if(xOnFullAD){
    inits$bb[5] <- 0
  }

  monitorPars <- c("effMeans", "b1", "b0", "res", "eff", "effVar", "res.mean", "res.var")

  if(constants$predAr){
    monitorPars <- c(monitorPars, "bArPred")
  }
  if(constants$predResVar&constants$randomRes){
    monitorPars <- c(monitorPars, "bResVarPred")
  }
  if(constants$predMean){
    monitorPars <- c(monitorPars, "bMeanPred")
  }
  if(constants$predX){
    monitorPars <- c(monitorPars, "bb", "xOutResVar")
  }
  if(constants$estimateAffectDynamics){
    monitorPars <- c(monitorPars, "mssd", "var")
  }
  #browser()

  # compileNimble(uppertri_mult_diag)
  buildMod <- nimbleModel(modelBaseline,
                          data = dataList,
                          constants=constants,
                          inits=inits)

  compMod <- compileNimble(buildMod)
  # browser()

  mcmcConfig <- configureMCMC(buildMod,print = TRUE, monitors = monitorPars)

  # mcmcConfig$removeSampler(c("eff"))
  # mcmcConfig$addSampler(type = 'AF_slice',
  #                       target=c("eff")),
  #                       control=list("sliceAdaptFactorInterval"=500))
  #   mcmcConfig$removeSampler(c("effMeans"))
  #
  #   mcmcConfig$removeSampler(c("bb"))
  #    mcmcConfig$removeSampler(c("eff"))
  #    mcmcConfig$removeSampler(c("Ustar"))
  #
  #
  #    mcmcConfig$addSampler(type = 'RW_block_lkj_corr_cholesky',
  #                          target=c("Ustar"),
  #                          control=list("adaptFactorExponent"=.7,
  #                                       "adaptInterval"=4000,
  #                                       "tries"=1))
  #
  #
  # # mcmcConfig$addSampler(type = 'RW_wishart',
  # #                       target=c("Ustar"),
  # #                       control=list("adaptFactorExponent"=.01,
  # #                                    "adaptInterval"=200,
  # #                                    "tries"=1))
  #
  #    mcmcConfig$addSampler(type = 'RW',
  #                          target=c("Ustar[3 ,3]"),
  #                          control=list("adaptFactorExponent"=.6,
  #                                       "adaptInterval"=300,
  #                                       "tries"=1))
  #    mcmcConfig$addSampler(type = 'RW',
  #                          target=c("Ustar[3 ,3]"),
  #                          control=list("adaptFactorExponent"=.1,
  #                                       "adaptInterval"=100,
  #                                       "tries"=1))
  #
  #
  #
 # mcmcConfig$removeSampler(c("effMeans",  "bb"))
  mcmcConfig$addSampler(type = 'RW_block',
                        target=c("bb"))
  mcmcConfig$removeSampler(c("res.mean"))
  mcmcConfig$addSampler(type = 'RW',
                        target=c("res.mean"),
                        control=list(reflective=FALSE,
                                     adaptInterval=1000))


  #
  #
  # mcmcConfig$addSampler(type="RW",
  #                       target="effMeans[3]",
  #                       control=list("tries"=2,
  #                                  "scale"=inits$effMeans[3]))
  #

  # eff <- c()
  # res <- c()
  #
  # mcmcConfig$removeSampler(c("eff"))
  # for(i in 1:N){
  #    eff[i] <- paste("eff[",i," ,", "1:2]", sep="")
  #    res[i] <- paste("res[",i,"]", sep="")
  #
  # mcmcConfig$addSampler(type = 'RW_block',
  #                       target=c(eff[i], res[i]),
  #                       control=list("adaptFactorExponent"=.8,
  #                                    "adaptInterval"=200,
  #                                    "propCov"=diag(c(inits$b0[i], inits$b1[i], inits$res[i]))))

  # mcmcConfig$addSampler(type="RW",
  #                       target=eff3[i],
  #                       control=list(scale=inits$eff[i, 3],
  #                                    tries=2,
  #                                    scale=inits$eff[i,3]))
  #}


  #   #
  # mcmcConfig$removeSampler(c("eff", "res"))
  # mcmcConfig$addSampler(type = 'RW_block',
  #                       target=c("eff", "res"),
  #                       control=list(tries=3,
  #                                    adaptInterval=100,
  #                                    adaptFactorExponent=.2))

  # if(constants$predX){
  mcmcConfig$removeSampler(c("sds", "res.var", "xOutResVar"))
  # mcmcConfig$addSampler(type = 'RW_block',
  #                       target=c("sds", "effC"))
  mcmcConfig$addSampler(type = 'RW',
                        target=("res.var"),
                        control=list(reflective=FALSE,
                                     scale=inits$res.var/2,
                                     adaptFactorExponent=.8,
                                     adaptInterval=200,
                                     tries=2))
  mcmcConfig$removeSampler(c("Ustar"))
  mcmcConfig$addSampler(type = "RW_block_lkj_corr_cholesky",
                        target=c("Ustar"),
                        control=list("adaptFactorExponent"=.9,
                        "adaptInterval"=1000,
                        "tries"=1))

  mcmcConfig$addSampler(type = 'RW',
                        target=("xOutResVar"),
                        control=list(reflective=FALSE,
                                     scale=inits$xOutResVar))
  mcmcConfig$addSampler(type = 'RW',
                        target=("sds[1]"),
                        control=list(reflective=FALSE,
                                     scale=inits$sds[1]))


  mcmcConfig$addSampler(type = 'RW',
                        target=("sds[2]"),
                        control=list(reflective=FALSE,
                                     scale=inits$sds[2]))

  # mcmcConfig$addSampler(type = 'RW',
  #                       target=("effC"),
  #                       control=list(scale=abs(inits$effC)))

  mcmcConfig$removeSampler(c("res"))
  for(i in 1:N){

    mcmcConfig$addSampler(type = 'RW',
                          target=paste("res[",i,"]", sep=""),
                          control=list(reflective=TRUE,
                                       scale=inits$res[i],
                                       adaptInterval=200))
  }

  #
  # mcmcConfig$addSampler(type = 'AF_slice',
  #                       target=c("res"))
  # }



  if(constants$predMean){
    mcmcConfig$removeSampler(c("bMeanPred"))
    mcmcConfig$addSampler(type = 'slice',
                          target=c("bMeanPred"),
                          control=list(
                            "adaptInterval"=500))

  }
  if(constants$predAr){
    mcmcConfig$removeSampler(c("bArPred"))
    mcmcConfig$addSampler(type = 'slice',
                          target=c("bArPred"),
                          control=list(
                            "adaptInterval"=500))

  }
  if(constants$predResVar){
    mcmcConfig$removeSampler(c("bResVarPred"))
    mcmcConfig$addSampler(type = 'slice',
                          target=c("bResVarPred"),
                          control=list(
                            "adaptInterval"=500))

  }





  if(buildOnly){
    mcmcMod <- buildMCMC(mcmcConfig, monitors = buildMod$getNodeNames(stochOnly = TRUE,
                                                                      includeData = FALSE))
  }else{

    mcmcMod <- buildMCMC(mcmcConfig, monitors=monitorPars)
  }


  cMcmcMod <- compileNimble(mcmcMod, project = buildMod)

  if(!buildOnly){
    start <- Sys.time()
    samplesList <- runMCMC(cMcmcMod, niter = niter,nburnin = nburnin,
                           nchains = nchains, inits = inits, thin=thin, summary=summary)

    finish <- Sys.time()
    runtime <- diff.Date(c(start,finish))
    print(paste("Completed sampling in", runtime))
    return(list(samplesList))
  }else{
    return(cMcmcMod)
  }
}

