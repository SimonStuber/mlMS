#' fit random mlMSAR model
#'
#'
#' @param y Data
#'
#'@import nimble mlVAR
#'
#' @return mcmc object
#'
#' @examples
#' fit_randomMLAR <-
#'
#' @export



fit_mlHMM <- function(simdat, M, randomEffects="transition_probabilities",
                        niter=10000, nburnin=.5, inits=NULL, thin=1,
                        nVars, nchains=2,forwardAlgorithm=TRUE,particleFilter=FALSE,
                        ncores=2, neff=10, rhat=1.1,resetSamplerRate=50){
  #browser()
  #parallel computing stuff adapted from:
  #https://github.com/qureshlatif/QSLpersonal/blob/master/R/RunNimbleParallel.R

  ni <- niter
  par.ignore.Rht = c()
  nc <- ncores
  nb <- nburnin
  nt  <- thin

  mod.nam = "mod"

  max.samples.saved = niter

  rtrn.model = TRUE
  sav.model = TRUE

  Rht.required <- rhat
  neff.required <- neff
  print("Building and compiling model.")
nVars <- nVars
    modelBaseline <- nimbleCode({
       for(i in 1:N){

         if(forwardAlgorithm){
           y[1:nTime, i, 1:nVars] ~ dForProbsHMM(probs1 = probs1[1:M],
                                                      y.hat.err=resCov2[i,1:(nVars*nVars),1:M],
                                                #      b=bRaneff[i,1:(nVars*nVars),1:M],
                                                      b1=b1[i,1:nVars,1:M],
                                                      # b1.1=b1.1.means[1:nVars,1:M],
                                                      # resCov1=resCov1.1[1:(nVars*nVars),1:M],
                                                      probs = probs[i,1:M, 1:M])
         }else{

           for(state in 1:M){
             b2[i,1:nVars, 1:nVars,state] <- bRaneff[i,1:(nVars*nVars),state]
             resCovFirstT[i,1:nVars, 1:nVars,state] <- (b2[i,1:nVars, 1:nVars,state]%*%
                                                          resCov[i,1:nVars, 1:nVars,state]%*%
                                                          t(b2[i,1:nVars, 1:nVars,state])+
                                                          resCov[i,1:nVars, 1:nVars,state])
           }



           y.hat[1, i, 1:nVars] <- b1[i,1:nVars,mm[1,i]]
           y[1, i, 1:nVars] ~ dmnorm(y.hat[1, i, 1:nVars],resCovFirstT[i,1:nVars, 1:nVars,mm[1,i]])


           for(t in 2:nTime){
             y.hat[t, i, 1:nVars] <- b1[i,1:nVars,mm[t,i]] + b2[i,1:nVars, 1:nVars,mm[t,i]]%*%(y[t-1, i, 1:nVars]-b1[i,1:nVars,mm[t,i]])
             y[t, i, 1:nVars] ~ dmnorm(y.hat[t, i, 1:nVars], resCov[i,1:nVars, 1:nVars,mm[t,i]])
           }


             for (t in 1:nTime){

               mm[t,i] ~ dcat(ps[i,t,1:M])
             }

             for (state in 1:M){

               ps[i,1,state] <- probs1[state]


               for (t in 2:nTime){
                 ps[i,t,state] <- probs[i, mm[(t-1),i],state]
               }
             }

         }


       }

      # for(state in 1:M){
      #   allMeans[1:((nVars*nVars)+nVars),state] <- c(b1.means[1:nVars,state],b.means[1:(nVars*nVars),state])
      # }

      for(i in 1:N){
        for(state in 1:M){
         # allB[i,1:((nVars*nVars)+nVars),state] ~ dmnorm(b)
         #   bRaneff[i,1:(nVars*nVars),state] ~ dmnorm(b.means[1:(nVars*nVars),state], cov=b.prec[1:(nVars*nVars),1:(nVars*nVars),state])
           # constraint_data ~ dconstraint(max(abs(bRaneff[i,1:(nVars*nVars),state]))< 1 )
          for(v in 1:nVars){
              b1[i,v,state] ~ dnorm(b1.means[v,state], var=b1.prec[v,state])
            }
          }
      }


      #    for(i in 1:N){
      # for(state in 1:M){
      #   #  y.hat.err[i, state, 1:nVars] ~ dmnorm(resMean[i,state, 1:nVars], resCov[i, state, 1:nVars, 1:nVars])
      #   resCov1[1:nVars, 1:nVars,state] ~ dinvwish(mat2[1:nVars, 1:nVars], nVars)
      #   # resMean[i,state, 1:nVars] ~ dmnorm(zeroVector[1:nVars], mat2[1:nVars, 1:nVars])
      # }
      #   }

      #      for(i in 1:N){
      # for(state in 1:M){
      #   resCov1.1[1:(nVars*nVars),state] <-  c(resCov1[1:nVars, 1:nVars,state])
      # }
      #     }

        for(state in 1:M){
          for(item in 1:nVars){
            b1.means[item,state] ~ dnorm(0,0.0000001)
          #  b1.1.means[item,state] ~ dnorm(0,0.0000001)
          }
        }

      # for(state in 1:M){
      #   for(eff in 1:(nVars*nVars)){
      #     b.means[eff,state] ~ T(dnorm(0,.001),-.99,.99)
      #   }
      # }






      for(state in 1:M){
     #   b.prec[1:(nVars*nVars),1:(nVars*nVars),state] ~ dinvwish(mat1[1:(nVars*nVars), 1:(nVars*nVars),state],(nVars*nVars))
        for(v in 1:nVars){
          b1.prec[v,state] ~ dinvgamma(.1, .1)
        }

      }

      # for (state in 1:M){
      #   shape[state] <- pow(res.mean[state],2) / pow(res.sd[state],2)
      #   rate[state] <- res.mean[state] / pow(res.sd[state],2)
      #   res.mean[state] ~ dunif(0,100)
      #   res.sd[state] ~ dunif(0,100)
      #   for(i in 1:N){
      #     y.hat.err[i, state] ~ dgamma(shape[state], rate[state])
      #   }
      # }

      # R_Means[v,i] ~ dnorm(0,1)
      #
      # Fisher_corr_noise[v,i] ~ dnorm(R_Means[i], pow(R_SD,-2))
      #
      #Correlation matrix
      # corr_noise[pp] = (exp(2*Fisher_corr_noise[pp])-1)/ (exp(2*Fisher_corr_noise[pp])+1)
      #


#       for(state in 1:M){
#       for(v in 1:nVars){
#         R_Means[v,state] ~dunif(0,100)
#         R_var[v,state] ~ dunif(0,100)
#       }
#
#       for(i in 1:N){
#         for(v in 1:nVars){
#           Fisher_corr_noise[i,v,state] ~ dnorm(R_Means[v,state], var=R_var[v,state])
#           corr_noise[v, i,state] <-(exp(2*Fisher_corr_noise[i,v,state])-1)/ (exp(2*Fisher_corr_noise[i,v,state])+1)
#           sd_noise_log[i,v,state] ~ dnorm(0,.001)
#           sd_noise[i, v,state] <- exp(sd_noise_log[i,v,state])
#         }
#       }
#
#       for(i in 1:N){
#      #   for(state in 1:M){
#        #   for(v1 in 1:nVars){
#         #    for(v2 in 1:nVars){
#          #     if(v1==v2){
#                 resCov[i,1,1,state] <- sd_noise[i, 1,state] * sd_noise[i, 1,state]
#                 resCov[i,2,2,state] <- sd_noise[i, 2,state] * sd_noise[i, 2,state]
#                 resCov[i,3,3,state] <- sd_noise[i, 3,state] * sd_noise[i, 3,state]
#
#           #    }else{
#                # idx <- idxMat[v1,v2]
#                 resCov[i,1,2,state] <- sd_noise[i, 1,state] * corr_noise[1,i,state] * sd_noise[i,2,state]
#                 resCov[i,1,3,state] <- sd_noise[i, 1,state] * corr_noise[2,i,state] * sd_noise[i,3,state]
#                 resCov[i,2,3,state] <- sd_noise[i, 2,state] * corr_noise[3,i,state] * sd_noise[i,3,state]
#                 resCov[i,2,1,state]<-resCov[i,1,2,state]
#                 resCov[i,3,1,state]<-resCov[i,1,3,state]
#                 resCov[i,3,2,state]<-resCov[i,2,3,state]
#            #   }
#             #}
#           #}
#  #       }
#       }
#
# }

  #    for(i in 1:N){
        for(state in 1:M){
        #  y.hat.err[i, state, 1:nVars] ~ dmnorm(resMean[i,state, 1:nVars], resCov[i, state, 1:nVars, 1:nVars])

          scaleResCov[1:nVars, 1:nVars, state] ~ dinvwish(matScale[1:nVars, 1:nVars,state], nVars)
          dfResCov[state] ~ dunif(nVars+4, N*M)

         meanResCov[1:nVars, 1:nVars,state] <- scaleResCov[1:nVars, 1:nVars, state]/(dfResCov[state]-nVars-1)
         for(i in 1:nVars){
           for(j in 1:nVars){
             varResCov[i,j,state] <- ((dfResCov[state]-nVars+1)*(scaleResCov[i, j, state]^2)+(dfResCov[state]-nVars-1)*scaleResCov[i, i, state]*scaleResCov[j,j, state])/
               ((dfResCov[state]-nVars)*((dfResCov[state]-nVars-1)^2)*(dfResCov[state]-nVars-3))
           }
         }
         for(i in 1:N){
           resCov[i,1:nVars, 1:nVars,state] ~ dinvwish(scaleResCov[1:nVars, 1:nVars, state], dfResCov[state])
         }
        }
   #   }

      for(i in 1:N){
        for(state in 1:M){
          resCov2[i,1:(nVars*nVars),state] <-  c(resCov[i,1:nVars, 1:nVars,state])
        }
      }





      for(state in 1:M){
        for (prev in 1:M){
          for(i in 1:N){
            odds[i,prev,state] <- exp(v[i,prev,state])
            probs[i, prev,state] <- odds[i,prev,state]/totalodds[i,prev]
          }
        }
      }

      for(i in 1:N){
        for(state in 1:M){
          totalodds[i, state] <- sum(odds[i,state,1:M])
        }
      }



      v[1:N, 1:M, 1:M] <- makeMat(raneffs[1:N,1:(M*M-M)],N,M)


      for(i in 1:N){
        raneffs[i,1:(M*M-M)] ~ dmnorm(v.means[1:(M*M-M)], cov=v.prec[1:(M*M-M),1:(M*M-M)])
      }



      T.constant.mu <- 0
      # location parameter for the prior
      T.constant.tau <- 1/T.constant.scale.squared
      # inverse scale parameter for the prior
      T.constant.scale.squared <- T.constant.scale * T.constant.scale
      T.constant.scale <- 10
      T.constant.k <- 1
      # 1 degree of freedom; this results in a Cauchy density.

      for (par in 1:(M*M-M)){ # loop over the number of fixed logit parameters
        v.means[par] ~ dt(T.constant.mu,T.constant.tau,T.constant.k)
      }

      v.prec[1:(M*M-M),1:(M*M-M)] ~ dinvwish(v.Om[1:(M*M-M),1:(M*M-M)],(M*M-M))


      probs1[1:M] ~ ddirich(start.alphas[1:M])
    })



n <- dim(simdat)[2]
nTime <- dim(simdat)[1]

   #   mat1 <- diag(nVars^2)

   # matScale <- diag(nVars)

    if(is.null(inits)){
      inits <- init_mlMSVAR(simdat,M, forwardAlgorithm = forwardAlgorithm, particleFilter = particleFilter)
    }

#mat1 <- inits$b.prec
matScale <- inits$scaleResCov
for(state in 1:M){
 # mat1[,,state] <- diag(diag(inits$b.prec[,,state]))
  matScale[,,state] <- diag(diag(inits$scaleResCov[,,state]))
}

dataList <- list("y"=simdat,
                 "start.alphas"=rep(1,M),
                 "v.Om" = diag((M*M-M)),
                # "mat1"=mat1,
                 "matScale"=matScale)
constantsList <- list("M"=M, "nTime"=dim(simdat)[1],
                      "N"=dim(simdat)[2],
                      "nVars"=nVars,
                      "forwardAlgorithm"=forwardAlgorithm)
parameters <- c("dfResCov", "b1.means", "meanResCov", "varResCov","v.means", "b1", "resCov", "probs", "probs1")

if(ncores==1){
  print("The model will be estimated with only 1 core. This is not adviced.
        To speed up estimation and reduce memory load, set ncores > 1")
}else{
  print("The models will be estimated with parallel computing.
        To speed up estimation and reduce memory load, the model will be estimated
        in chunks. New chunks will be started until convergence is reached.
        Chunks and final results are saved to the specified working directory.")
}


cl<-makeCluster(ncores, timeout = 5184000)
clusterExport(cl, c("modelBaseline", "inits", "dataList", "constantsList","parameters",
                    "niter", "thin","ni","par.ignore.Rht","nc","nb","mod.nam",
                    "max.samples.saved","rtrn.model","sav.model","Rht.required",
                    "neff.required",
                    "dForProbsHMM", "makeMat", "makeSingleMat", "lse",
                    "init_mlMSVAR", "particleFilter"),envir=environment())

for (j in seq_along(cl)) {
  set.seed(j)
  init <<- inits
  clusterExport(cl[j], "init")
}

print("Sampling startet. This can take many hours.
      Progress bar not available with parallel computing.
      However, every niter samples, convergence statistics will be printed.")

out1 <- clusterEvalQ(cl, {
  library(nimble)
  library(coda)
  buildMod <- nimbleModel(modelBaseline,
                          constants=constantsList,
                          data = dataList,
                          inits=init)
  compMod <- compileNimble(buildMod)


  if(particleFilter){

    stateSpaceMCMCconf <- configureMCMC(buildMod,control=list("adaptInterval"=500, "adaptFactorExponent"=.5))




    stateSpaceMCMCconf$removeSampler(c("dfResCov"))#,
    # stateSpaceMCMCconf$addSampler(type = 'RW_PF_block',
    #                            target=c("dfResCov"),
    #                            control=list("adaptInterval"=500,
    #                                         "adaptFactorExponent"=.5,
    #                                         latents=c("mm"),
    #                                         timeIndex=1,
    #                                         pfType="bootstrap",
    #                                         pfNparticles=10),
    #                            pfControl=list(thresh = 0.8, saveAll = TRUE,
    #                                           smoothing = FALSE,
    #                                           timeIndex=1))
    stateSpaceMCMCconf$addSampler(type = 'AF_slice',
                                  target=c("dfResCov"),
                                  control=list("sliceAdaptFactorMaxIter"=max(15000,niter),
                                               "sliceAdaptWidthMaxIter"=max(512,niter),
                                               "sliceAdaptWidthTolerance"=.01,
                                               "sliceAdaptFactorInterval"=500,
                                               "sliceMaxSteps"=500,
                                               "maxContractions"=500))


    stateSpaceMCMCconf$removeSampler(c("resCov", "bRaneff", "b1", "scaleResCov"))#,


    b1Target <- c()
    b1MeansTarget <- c()
    b1PrecTarget <- c()




    for(i in 1:dim(simdat)[2]){
      b1Target[i] <- paste("b1[", i, ", 1:", nVars, ", 1:", M, "]", sep="")
    }

    for(i in 1:dim(simdat)[2]){
      # for(state in 1:M){
      stateSpaceMCMCconf$addSampler(type = 'RW_PF_block',
                                    target=c(b1Target[i]),
                                    control=list("adaptInterval"=500,
                                                 "adaptFactorExponent"=.5,
                                                 latents=c("mm"),
                                                 timeIndex=1,
                                                 pfType="bootstrap",
                                                 pfNparticles=10),
                                    pfControl=list(thresh = 0.8, saveAll = TRUE,
                                                   smoothing = FALSE,
                                                   timeIndex=1))
      # }
    }

    bRaneffTarget <- c()
    bMeansTarget <- c()
    bPrecTarget <- c()


    for(i in 1:dim(simdat)[2]){
      bRaneffTarget[i] <- paste("bRaneff[",i, ", 1:", nVars*nVars, ", 1:", M, "]", sep="")
    }


    for(i in 1:dim(simdat)[2]){
      # for(state in 1:M){
      stateSpaceMCMCconf$addSampler(type = 'RW_PF_block',
                                    target=c(bRaneffTarget[i]),
                                    control=list("adaptInterval"=500,
                                                 "adaptFactorExponent"=.5,
                                                 latents=c("mm"),
                                                 timeIndex=1,
                                                 pfType="bootstrap",
                                                 pfNparticles=10),
                                    pfControl=list(thresh = 0.8, saveAll = TRUE,
                                                   smoothing = FALSE,
                                                   timeIndex=1))
      #  }
    }






    target <- array(NA, c(dim(simdat)[2],M))
    for(i in 1:dim(simdat)[2]){
      for(state in 1:M){
        target[i, state] <- paste("resCov[", i, ", 1:", nVars, ", 1:", nVars, ", ", state, "]", sep="")
      }
    }
    for(i in 1:dim(simdat)[2]){
      for(state in 1:M){
        stateSpaceMCMCconf$addSampler(type = 'RW_wishart',
                                      target=target[i,state],
                                      control=list("adaptInterval"=500,
                                                   # "scale"=(1-(mean(abs(bRInits)))^2)/var(as.vector(simdat)),
                                                   "adaptFactorExponent"=.5))
      }
    }

    scaleTarget <- c()
    for(state in 1:M){
      scaleTarget[state] <- paste("scaleResCov[","1:",nVars,", 1:",nVars,", ",state, "]", sep="")

      stateSpaceMCMCconf$addSampler(type = 'RW_wishart',
                                    target=scaleTarget[state],
                                    control=list("adaptInterval"=500,
                                                 "scale"=(var(simdat)/n)*nTime,
                                                 "adaptFactorExponent"=.5))
    }



    stateSpaceMCMCconf$addMonitors(c('meanResCov', 'varResCov'))


    mcmcMod<- buildMCMC(stateSpaceMCMCconf)


    cMcmcMod <- compileNimble(mcmcMod, project = buildMod,resetFunctions = TRUE)


    start <- Sys.time()
    samplesList <- runMCMC(cMcmcMod, niter = niter,nburnin = nburnin,
                           nchains = nchains, inits = inits, thin=thin,samples = TRUE)

    finish <- Sys.time()
    runtime <- diff.Date(c(start,finish))



    # C_AFSS_mcmc <- compileNimble(AFSS_mcmc, project = buildMod)
    #
    #
    #
    # cat("Start MCMC sampling...This may even take (much) more than a while.")
    #
    #
    # st <- system.time(C_AFSS_mcmc$run(niter = niter, nburnin = nburnin, thin=1))


    # AFSS_mcmcOutput <- as.mcmc(as.matrix(C_AFSS_mcmc$mvSamples))
    #
    #
    # sum <- summary(AFSS_mcmcOutput)



    print(paste("Completed sampling in", runtime))



    return(samplesList)


  }else{

    AFSS_mcmcConfig <- configureMCMC(buildMod,control=list("adaptInterval"=500, "adaptFactorExponent"=.5))

    AFSS_mcmcConfig$removeSampler(c("dfResCov"))

    AFSS_mcmcConfig$addSampler(type = 'AF_slice',
                               target=c("dfResCov"),
                               control=list("sliceAdaptFactorMaxIter"=max(15000,niter/2),
                                            "sliceAdaptWidthMaxIter"=max(512,niter/2),
                                            "sliceAdaptWidthTolerance"=.01,
                                            "sliceAdaptFactorInterval"=500,
                                            "sliceMaxSteps"=500,
                                            "maxContractions"=500))


    AFSS_mcmcConfig$removeSampler(c("resCov", "b1", "scaleResCov"))#,

    b1Target <- matrix(NA, constantsList$N,constantsList$M)
    b1MeansTarget <- c()
    b1PrecTarget <- c()

    for(state in 1:constantsList$M){
      for(i in 1:constantsList$N){
        b1Target[i,state] <- paste("b1[", i, ", 1:", constantsList$nVars, ", ", state, "]", sep="")
      }
    }

    for(i in 1:constantsList$N){
      for(state in 1:constantsList$M){
        AFSS_mcmcConfig$addSampler(type = 'RW_block',
                                   target=c(b1Target[i,state]),
                                   control=list("adaptInterval"=500,
                                                "adaptFactorExponent"=.5,
                                                "scale"=mean(inits$b1.prec)))
      }
    }

    #bRaneffTarget <- matrix(NA, constantsList$N,constantsList$M)
    bMeansTarget <- c()
    bPrecTarget <- c()

    # for(state in 1:constantsList$M){
    #   for(i in 1:constantsList$N){
    #     bRaneffTarget[i,state] <- paste("bRaneff[",i, ", 1:", constantsList$nVars*constantsList$nVars, ", ", state, "]", sep="")
    #   }
    # }
    #
    #
    # for(i in 1:constantsList$N){
    #   for(state in 1:constantsList$M){
    #     AFSS_mcmcConfig$addSampler(type = 'RW_block',
    #                                target=c(bRaneffTarget[i,state]),
    #                                control=list("adaptInterval"=500,
    #                                             "adaptFactorExponent"=.5))
    #   }
    # }

    target <- array(NA, c(constantsList$N,constantsList$M))
    for(i in 1:constantsList$N){
      for(state in 1:constantsList$M){
        target[i, state] <- paste("resCov[", i, ", 1:", constantsList$nVars, ", 1:", constantsList$nVars, ", ", state, "]", sep="")
      }
    }
    for(i in 1:constantsList$N){
      for(state in 1:constantsList$M){
        AFSS_mcmcConfig$addSampler(type = 'RW_wishart',
                                   target=target[i,state],
                                   control=list("adaptInterval"=500,
                                                "scale"=mean(abs(inits$resCov[i,,,state])),
                                                "adaptFactorExponent"=.5))
      }
    }

    scaleTarget <- c()
    for(state in 1:constantsList$M){
      scaleTarget[state] <- paste("scaleResCov[","1:",constantsList$nVars,", 1:",constantsList$nVars,", ",state, "]", sep="")

      AFSS_mcmcConfig$addSampler(type = 'RW_wishart',
                                 target=scaleTarget[state],
                                 control=list("adaptInterval"=500,
                                              "scale"=mean(abs(inits$scaleResCov[,,state])),
                                              "adaptFactorExponent"=.5))
    }

    AFSS_mcmcConfig$addMonitors(c('meanResCov', 'varResCov',"b1", "resCov", "probs","probs1"))

    mcmcMod<- buildMCMC(AFSS_mcmcConfig)

    cMcmcMod <- compileNimble(mcmcMod, project = buildMod)


    start <- Sys.time()

    cMcmcMod$run(niter, reset = FALSE)
    finish <- Sys.time()
    runtime <- diff.Date(c(start,finish))
    print(paste("Completed sampling in", runtime))
    return(as.mcmc(as.matrix(cMcmcMod$mvSamples)))

  }

})

for(chn in 1:ncores) {
  ind.keep <- c()
  for(p in 1:length(parameters)) ind.keep <-
      c(ind.keep, which(str_detect(dimnames(out1[[chn]])[[2]], parameters[p]))) %>% unique()
  out1[[chn]] <- out1[[chn]][,ind.keep]
}

## Check convergence ##
out2 <- out1
ni.saved <- nrow(out2[[1]])
for(chn in 1:ncores) { # nc must be > 1
  out2[[chn]] <- out2[[chn]][(round(ni.saved * nburnin)+1):ni.saved,]
}
out.mcmc <- coda::as.mcmc.list(lapply(out2, coda::as.mcmc))

mod <- mcmcOutput(out.mcmc)
sumTab <- summary(mod, MCEpc = F, Rhat = T, n.eff = T, f = T, overlap0 = T, verbose = F)
sumTab <- sumTab %>%
  as_tibble() %>%
  mutate(Parameter = row.names(sumTab)) %>%
  select(Parameter, mean:f)


if(length(par.ignore.Rht) == 0) {
  mxRht <- sumTab %>% pull(Rhat) %>% max(na.rm = T)
  mn.neff <- sumTab %>% pull(n.eff) %>% min(na.rm = T)
} else {
  ind.ignore <- c()
  for(p in 1:length(par.ignore.Rht)) ind.ignore <-
      c(ind.ignore, which(str_detect(sumTab$Parameter, par.ignore.Rht[p]))) %>%
      unique()
  mxRht <- sumTab %>% slice(-ind.ignore) %>% pull(Rhat) %>% max(na.rm = T)
  mn.neff <- sumTab %>% slice(-ind.ignore) %>% pull(n.eff) %>% min(na.rm = T)
}



mod <- list(mcmcOutput = mod, summary = sumTab)
if(sav.model) R.utils::saveObject(mod, mod.nam) # If running all in one.

## If has not converged, continue sampling
if(round(mxRht, digits = 1) > Rht.required | mn.neff < neff.required) {
  n.runs <- 1

  R.utils::saveObject(out1, str_c(mod.nam, "_chunk", n.runs)) # Save samples from previous run to drive.
}
while(round(mxRht, digits = 2) > Rht.required | mn.neff < neff.required) {
  n.runs <- n.runs + 1
  print(str_c("Run = ", n.runs, ". Max Rhat = ", mxRht, " and min neff = ", mn.neff))



    if((n.runs/resetSamplerRate)%%1==0){
      out2 <- clusterEvalQ(cl, {
        cMcmcMod$run(niter, reset = TRUE, resetMV = TRUE) # Resume sampling.
        return(as.mcmc(as.matrix(cMcmcMod$mvSamples)))
        gc(verbose = F)
      })
    }else{
      out2 <- clusterEvalQ(cl, {
        cMcmcMod$run(niter, reset = FALSE, resetMV = TRUE) # Resume sampling.
        return(as.mcmc(as.matrix(cMcmcMod$mvSamples)))
        gc(verbose = F)
      })
    }


  for(chn in 1:ncores) { # ncores must be > 1
    ind.keep <- c()
    for(p in 1:length(parameters)) ind.keep <-
        c(ind.keep, which(str_detect(dimnames(out2[[chn]])[[2]], parameters[p]))) %>% unique()
    out2[[chn]] <- out2[[chn]][,ind.keep]
  }
  R.utils::saveObject(out2, str_c(mod.nam, "_chunk", n.runs)) # Save samples from previous run to drive.

  ni2 <- round(((ni / nt) * n.runs * ncores) * nb) # Anticipated number of samples to save (assuming half discarded as burn-in).
  if(ni2 > max.samples.saved) {
    nt2 <- round(1 / (max.samples.saved / ni2)) # Set additional thinning so that saved iterations don't exceed (by too much) max.samples.saved (specified by user).
  } else {
    nt2 <- 1
  }

  # Reassemble chain from chunks and apply additional thinning.
  out1 <- R.utils::loadObject(str_c(mod.nam, "_chunk", 1))
  ni.saved <- nrow(out1[[1]])
  for(chn in 1:ncores) { # ncores must be > 1
    out1[[chn]] <- out1[[chn]][seq(2, ni.saved, by = nt2),] # Starting at 2 because first iteration is NA for some reason.
  }
  for(r in 2:n.runs) {
    out.r <- R.utils::loadObject(str_c(mod.nam, "_chunk", r))
    ni.saved <- nrow(out.r[[1]])
    for(chn in 1:ncores) {
      out.r[[chn]] <- out.r[[chn]][seq(1, ni.saved, by = nt2),]
      out1[[chn]] <- rbind(out1[[chn]], out.r[[chn]])
    }
  }

  # Discard specified proportion of initial samples as burn-in
  out3 <- out1
  ni.saved <- nrow(out3[[1]])
  for(chn in 1:ncores) {
    out3[[chn]] <- out3[[chn]][(round(ni.saved/2)+1):ni.saved,]
  }
  out.mcmc.update <- coda::as.mcmc.list(lapply(out3, coda::as.mcmc))

  mod <- mcmcOutput(out.mcmc.update)
  sumTab <- summary(mod, MCEpc = F, Rhat = T, n.eff = T, f = T, overlap0 = T, verbose = F)
  sumTab <- sumTab %>%
    as_tibble() %>%
    mutate(Parameter = row.names(sumTab)) %>%
    select(Parameter, mean:f)
  if(length(par.ignore.Rht) == 0) {
    mxRht <- sumTab  %>% pull(Rhat) %>% max(na.rm = T)
    mn.neff <- sumTab %>% pull(n.eff) %>% min(na.rm = T)
  } else {
    mxRht <- sumTab %>% slice(-ind.ignore) %>% pull(Rhat) %>% max(na.rm = T)
    mn.neff <- sumTab %>% slice(-ind.ignore) %>% pull(n.eff) %>% min(na.rm = T)
  }
  gc(verbose = F)

  mod <- list(mcmcOutput = mod, summary = sumTab)
  if(sav.model) R.utils::saveObject(mod, mod.nam) # If running all in one.
}
if(exists("n.runs")) for(r in 1:n.runs) file.remove(str_c(mod.nam, "_chunk", r))
stopCluster(cl)
if(rtrn.model) return(mod)


#
# posteriorSamples <- as.matrix(compiledList$bootstrapFilter$mvEWSamples)
#
# cBootF$run(30)


}


