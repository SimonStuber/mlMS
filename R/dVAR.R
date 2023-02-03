dVar <- nimbleFunction(
  run = function(x=double(2),
                # probs1=double(1),
                 y.hat.err=double(2),
                 b=double(1),
                 b1=double(1),
             #    b1.1=double(2),
           #      resCov1=double(2),
                 #probs=double(2),
                 log = integer(0)){



    y <- x
    nVars <- dim(y)[2]
    nTime <- dim(y)[1]
 #   M <- length(probs1)
    # b2 <- array(0, c(nVars, nVars,M))



   # b2 <- b[1:(nVars*nVars)]
    #   for(state in 1:M){
    b2 <- nimArray(value=nimNumeric(length=length(b[1:(nVars*nVars)]),
                                    value=b[1:(nVars*nVars)]),
                   dim=c(nVars,nVars))

resCov <- y.hat.err[1:nVars, 1:nVars]

    # resCov <- nimArray(value=nimNumeric(length=length(y.hat.err[1:(nVars*nVars)]),
    #                                     value=y.hat.err[1:(nVars*nVars)]),
    #                    dim=c(nVars,nVars))

    # resCov1.1 <- nimArray(value=nimNumeric(length=length(resCov1[1:(nVars*nVars),1:M]),
    #                                     value=resCov1[1:(nVars*nVars),1:M]),
    #                    dim=c(nVars,nVars,M))
    #  }
    #  }

    #gamma <- probs1
    #tp <- probs
    #logalpha <- array(0,c(nTime,M))
    dens <- array(0,c(nTime))
    y.hat <- array(0, c(nTime, nVars))
    #sum <- 0
    #total <- 0
#    for(mm in 1:M){
      y.hat[1, 1:nVars] <- b1[1:nVars]
      dens[1] <- dmnorm_chol(x=y[1, 1:nVars], mean = y.hat[1, 1:nVars], chol(b2[1:nVars, 1:nVars]%*%
                                                                                     resCov[1:nVars, 1:nVars]%*%
                                                                                     t(b2[1:nVars, 1:nVars])+
                                                                                     resCov[1:nVars, 1:nVars]), prec_param=FALSE, log=TRUE)
  #  }

    for(t in 2:nTime){
     # for(mm in 1:M){

        y.hat[t, 1:nVars] <- t(b1[1:nVars]) + t(b2[1:nVars, 1:nVars]%*%(y[t-1, 1:nVars]-b1[1:nVars]))
        dens[t] <- dmnorm_chol(x=y[t, 1:nVars], mean = y.hat[t, 1:nVars], chol(resCov[1:nVars, 1:nVars]), prec_param=FALSE, log=TRUE)
      }#}


    #dens <- dens/max(dens)


    #lalpha   <- matrix(0, nTime,M)
    #alpha_prob <- matrix(0, nTime,M)
    #allprobs <- dens
    #delta <- probs1

    # foo             <- delta[1:M] * allprobs[1, 1:M]
    # sumfoo          <- sum(foo)
    # alpha_prob[1, 1:M] <- foo/sumfoo
    # lscale          <- log(sumfoo)
    # lalpha[1,1:M]     <- log(alpha_prob[1,1:M]) + lscale
    # for (i in 2:nTime){
    #   foo[1:M]              <- alpha_prob[(i - 1),1:M]%*%tp[1:M,1:M]*allprobs[i,1:M]
    #   sumfoo           <- sum(foo)
    #   alpha_prob[i,1:M]  <- foo / sumfoo
    #   lscale           <- lscale + log(sumfoo)
    #   lalpha[i,1:M]       <- log(alpha_prob[i,1:M]) + lscale
    # }

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
    if (log)return(sum(dens[1:nTime]))
    else return(exp(sum(dens[1:nTime])))
    returnType(double(0))
  }
)


registerDistributions(list("dVar"=list(BUGSdist="dVar(probs1,y.hat.err,b,b1,probs)",
                                                    types=c('value=double(2)',
                                                       #     'probs1=double(1)',
                                                            'y.hat.err=double(2)',
                                                            'b=double(1)',
                                                            'b1=double(1)'))))
                                                         #   'b1.1=double(2)',
                                                         #   'resCov1=double(2)''))))

