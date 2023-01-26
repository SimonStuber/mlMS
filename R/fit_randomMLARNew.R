library(nimble)
library(lme4)
library(BBmisc)

#step 1: install nimble and all requirements
#following https://r-nimble.org/html_manual/cha-installing-nimble.html

#step 2: source this R-script

#step 3: store data of of the uncentered dependent variable
#in a matrix with dimensions nTime*N

#step 4: (optional) store data of uncentered Lag-2
#predictors in separate vectors
#for the prediction of the autocorrelation (arOnX)
#or the random residual variance (resVarOnX)

#step 5: specify reasonable initial values and store them in a list.
#
#here, we start with starting initial for the
#covariance matrix of the (specified) random effects
#for instance, for a model with random intercepts,
#random slopes and random residual variances, initial values can be specified
#in the following way:
# effVar <- diag(c(0,0,0))
# effVar[1,1] <- 1 # assume a variance of 1 for the random intercept
# effVar[2,2] <- 0.1 # assume a variance of .1 for the random autocorrelations
# effVar[3,3] <- log(3) # assume a variance of 3,
#for the random residual variance
# note that the variance for the random residual variance needs
#to be log-transformed
#
#we continue with starting values for the means of the random effects
#effMeans=c(10,.1, log(3)) # for an average intercept of 10,
#an average autocorrelation of .1 and an average residual variance of 3
#
#next, we specify initial values for the effect of the Level-2 predictors
#
# bResVarPred=0
# bArPred=0
# bMeanPred=0
#
#and for the mean at the first measurement occasion
# b0Start=rep(10, N)
#
#now, store all initial values in a list
# inits <- list(effMeans=effMeans,
#            effVar=effVar,
#            bResVarPred=bResVarPred,
#            bArPred=bArPred,
#           bMeanPred=bMeanPred,
#            b0Start=b0Start)


#step 6: now, we can specify certain constants, that define the model.
#This information needs
#to be stored in a list and provided to the fit_randomMLAR function
#through the constants argument
#
# constants <- list(N=dim(y)[2],
##number of individuals
#                nTime=dim(y)[1],
## maximum number of time points
#                predAr=1,
##should the random autocorrelation be predicted by a covaraite (0=no, 1=yes)?
#                predResVar=1,
##should the random residual variance be predicted by a covaraite (0=no, 1=yes)?
#                predMean=1,
##should the random intercepts be predicted by a covaraite (0=no, 1=yes)?
#                randomRes=1,
##should random effects in the residual variance be estimated (0=no, 1=yes)?
#                lastT=lastT,
## a vector with the first measurement occasion for every individual
#                firstT=firstT
## a vector with the last measurement occasion for every inidividual
#)

# #step 7: fit the model
# res <- fit_randomMLAR(y, arOnX=x,
#                       resVarOnX=x,
#                       meanOnX=x,
#                       niter = 60000,
#                       nburnin = 30000,
#                       inits=inits,
#                       constants=constants)

#step 8: inspect results
#summary(res)



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
          res[i] <- exp( eff[i,3]+ bResVarPred*resVarOnX_c[i])
        }else{
          res[i] <- exp(eff[i,3])
        }
        eff[i,1:3] ~ dmnorm(effMeans[1:3], cholesky = U[1:3, 1:3], prec_param = 0)
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
      xOutResVar ~ dunif(0,100)
      if(fullAD){
        for(i in 1:N){
          xOutHat[i] <-  bb[1] + bb[2]*b0[i] + bb[3]*b1[i] + bb[4]*mssd[i] +bb[5]*var[i]
          xOutcome[i] ~ dnorm(xOutHat[i], sd=xOutResVar)
        }

      }else{
        for(i in 1:N){
          xOutHat[i] <-  bb[1] + bb[2]*(b0[i]-mean(b0[1:N])) + bb[3]*(b1[i]-mean(b1[1:N])) + bb[4]*(res[i]-mean(res[1:N]))
          xOutcome[i] ~ dnorm(xOutHat[i], sd=xOutResVar)
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
      effMeans[3] ~ dnorm(0, 0.001)
   #   effPrec[1:3,1:3] ~ dwish(effPrecPriorMat[1:3, 1:3], 3)
      Ustar[1:3,1:3] ~ dlkj_corr_cholesky(1, 3)
     # for(nsd in 1:3){
        sds[1] ~ dunif(0,100)
        sds[2] ~ dunif(0,2)
        sds[3] ~ dunif(0,100)
      #}

U[1:3,1:3] <- uppertri_mult_diag(Ustar[1:3, 1:3], sds[1:3])


     effVar[1:3,1:3] <-U[1:3,1:3]%*%t(U[1:3,1:3])
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
    inits <- init_randomMLAR(y, xOutcome, nTime, constants)
  }
  if(xOnFullAD){
    inits$bb[5] <- 0
  }

  monitorPars <- c("effMeans", "b1", "b0", "res", "eff", "effVar")

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

  mcmcConfig <- configureMCMC(buildMod,print = FALSE, monitors = monitorPars)

  # mcmcConfig$removeSampler(c("eff"))
  # mcmcConfig$addSampler(type = 'AF_slice',
  #                       target=c("eff")),
  #                       control=list("sliceAdaptFactorInterval"=500))
  mcmcConfig$removeSampler(c("effMeans"))
  mcmcConfig$addSampler(type = 'RW_block',
                        target=c("effMeans"),
                        control=list(tries=1,
                                     propCov=inits$U%*%t(inits$U),
                                     adaptInterval=300,
                                     adaptFactorExponent=.8))

   mcmcConfig$removeSampler(c("eff"))
  eff <- c()
  for(i in 1:N){
    eff[i] <- paste("eff[",i," ,", "1:3]", sep="")
    mcmcConfig$addSampler(type = 'RW_block',
                          target=c(eff[i]),
                          control=list(tries=2,
                                       propCov=diag(inits$eff[i,]),
                                       adaptFactorExponent=.8,
                                       adaptInterval=500))
  }

  #
  # # mcmcConfig$removeSampler(c("eff"))
  # # mcmcConfig$addSampler(type = 'RW_block',
  # #                       target=c("eff"),
  # #                       control=list(tries=3,
  # #                                    adaptInterval=100,
  # #                                    adaptFactorExponent=.2))
  # #
  if(constants$predX){
    mcmcConfig$removeSampler(c("bb"))
    mcmcConfig$addSampler(type = 'RW_block',
                          target=c("bb"),
                          control=list(tries=1,
                                       adaptInterval=500,
                                       adaptFactorExponent=.8))


   }



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
    return(list(samplesList, cMcmcMod))
  }else{
    return(cMcmcMod)
  }
}

