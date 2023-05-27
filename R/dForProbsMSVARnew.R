dForProbsMSVARnew <- nimbleFunction(
  run = function(x=double(2),
                 probs1=double(1),
                 y.hat.err=double(2),
                 b=double(2),
                 b1=double(2),
             #    b1.1=double(2),
           #      resCov1=double(2),
                 probs=double(2),
                 log = integer(0)){



    y <- x
    nVars <- dim(y)[2]
    nTime <- dim(y)[1]
    M <- length(probs1)
    # b2 <- array(0, c(nVars, nVars,M))




    #   for(state in 1:M){
    b2 <- nimArray(value=nimNumeric(length=length(b[1:(nVars*nVars),1:M]),
                                    value=b[1:(nVars*nVars),1:M]),
                   dim=c(nVars,nVars,M))



    resCov <- nimArray(value=nimNumeric(length=length(y.hat.err[1:(nVars*nVars),1:M]),
                                        value=y.hat.err[1:(nVars*nVars),1:M]),
                       dim=c(nVars,nVars,M))

    # resCov1.1 <- nimArray(value=nimNumeric(length=length(resCov1[1:(nVars*nVars),1:M]),
    #                                     value=resCov1[1:(nVars*nVars),1:M]),
    #                    dim=c(nVars,nVars,M))
    #  }
    #  }

    gamma <- probs1
    tp <- probs
    logalpha <- array(0,c(nTime,M))
    dens <- array(0,c(nTime,M))
    y.hat <- array(0, c(nTime, nVars,M))
    #sum <- 0
    #total <- 0
    for(mm in 1:M){
      y.hat[1, 1:nVars,mm] <- b1[1:nVars,mm]
      dens[1,mm] <- dmnorm_chol(x=y[1, 1:nVars], mean = y.hat[1, 1:nVars,mm], chol(b2[1:nVars, 1:nVars,mm]%*%
                                                                                     resCov[1:nVars, 1:nVars,mm]%*%
                                                                                     t(b2[1:nVars, 1:nVars,mm])+
                                                                                     resCov[1:nVars, 1:nVars,mm]), prec_param=FALSE, log=FALSE)
    }

    for(t in 2:nTime){
      for(mm in 1:M){

        y.hat[t, 1:nVars,mm] <- t(b1[1:nVars,mm]) + t(b2[1:nVars, 1:nVars,mm]%*%(y[t-1, 1:nVars]-b1[1:nVars,mm]))
        dens[t,mm] <- dmnorm_chol(x=y[t, 1:nVars], mean = y.hat[t, 1:nVars,mm], chol(resCov[1:nVars, 1:nVars,mm]), prec_param=FALSE, log=FALSE)
      }}


    #dens <- dens/max(dens)


    lalpha   <- matrix(0, nTime,M)
    alpha_prob <- matrix(0, nTime,M)
    allprobs <- dens
    delta <- probs1
#check zucchini 2016
    foo             <- delta[1:M] * allprobs[1, 1:M]
    sumfoo          <- sum(foo)
    alpha_prob[1, 1:M] <- foo/sumfoo
    lscale          <- log(sumfoo)
    lalpha[1,1:M]     <- log(alpha_prob[1,1:M]) + lscale
    for (i in 2:nTime){
      foo[1:M]              <- alpha_prob[(i - 1),1:M]%*%tp[1:M,1:M]*allprobs[i,1:M]
      sumfoo           <- sum(foo)
      alpha_prob[i,1:M]  <- foo / sumfoo
      lscale           <- lscale + log(sumfoo)
      lalpha[i,1:M]       <- log(alpha_prob[i,1:M]) + lscale
    }

    #   for(n in 1:N){

    # sumlog <- log(M)
    # for(mm in 1:M){
    #   dens[mm,t] <- dens[mm,t]-sumlog
    # }


    #}
    #
    # for(s in 1:M){
    #   logalpha[1,s] <- lse(c(dens[s,1],log(probs1[s])))
    # }
    # # sumlog <- nTime
    # #  for(mm in 1:M){
    # #    logalpha[1, mm] <-  logalpha[1, mm]-log(M)
    # #  }
    #
    #
    # for (t in 2:nTime) {
    #   for (j in 1:M) {
    #     accumulator <- array(NA, c(M))
    #     for (i in 1:M) {
    #       accumulator[i] = lse(c(logalpha[t-1, i],log(tp[i, j]),dens[j,t]))
    #     }
    #     # print(accumulator)
    #
    #     logalpha[t, j] = lse(c(accumulator))#-log(M)
    #     # logalpha[t, j] <-  logalpha[t, j]-sumlog
    #
    #
    #
    #   }
    #   #  sumlog <- lse(logalpha[t, 1:M])
    #   # for(mm in 1:M){
    #   #   logalpha[t, mm] <-  logalpha[t, mm]-sumlog
    #   # }
    #
    # }
    #
    #
    #
    #
    #
    if (log)return(lse(lalpha[nTime,1:M]))
    else return(exp(lse(lalpha[ nTime,1:M])))
    returnType(double(0))
  }
)


registerDistributions(list("dForProbsMSVARnew"=list(BUGSdist="dForProbsMSVARnew(probs1,y.hat.err,b,b1,probs)",
                                                    types=c('value=double(2)',
                                                            'probs1=double(1)',
                                                            'y.hat.err=double(2)',
                                                            'b=double(2)',
                                                            'b1=double(2)',
                                                         #   'b1.1=double(2)',
                                                         #   'resCov1=double(2)',

                                                            'probs=double(2)'))))

