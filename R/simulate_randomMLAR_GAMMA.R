simulate_randomMLAR_G <- function(n, nTime,
                                  b0.means, b1.means,b0.sd,b1.sd,res.mean,res.var,
                                  correlations){

  nTime <- nTime +100
  #browser()
  y <- array(NA, c(nTime, n))
  raneffs <- array(NA, c(n, 4))
  b0 <- c()
  b1 <- c()
  res <- c()
  x <- c()
  #
  #   if(is.null(covMat)){
  #   rate <- c()
  #   shape <-c()
  #
  #
  #   b0 <- rnorm(n = n,b0.means, b0.sd)
  #   b1 <- rnorm(n = n,b1.means, b1.sd)
  #
  #   shape <- (res.mean^2)/(res.var^2)
  #   rate <- (res.mean)/(res.var^2)
  #   res <- rgamma(n, shape,rate)
  #
  #
  #
  #
  #   for(m in 1:n){
  #
  #
  #     y[1,m] <- rnorm(1,b0[m],1/sqrt(res[m]))
  #
  #
  #       if(b1[m]<(-1)){
  #         b1[m] <- -.99999
  #       }else if(b1[m]>(1)){
  #         b1[m] <- .99999
  #       }
  #
  #
  #     for(t in 2:nTime){
  #
  #       y[t,m] <- b0[m] + b1[m]*(y[t-1,m]-b0[m]) + rnorm(1,0,1/sqrt(res[m]))
  #
  #     }
  #   }
  #
  #   res <- list(y=y,b0=b0,b=b1,res=res)
  #   return(res)
  #   }
  #   else{


  phi <- (res.mean^2)/res.var
  mu <- res.mean
  shape <- phi+2
  scale <- mu*(1+phi)
  #browser()

  i <- 1
  shrinked <- FALSE

  while(i<=n) {
    # print(i)
    #raneffs[i,] <- mvtnorm::rmvnorm(1,effMeans, effVar)

    myCop <- normalCopula(param=correlations, dim = 4, dispstr = "un")
    myMvd <- mvdc(copula=myCop, margins=c("norm", "norm", "invgamma", "norm"),
                  paramMargins=list(list(mean=b0.means, sd=b0.sd),
                                    list(mean=b1.means, sd=b1.sd),
                                    list(scale=scale,shape=shape),
                                    list(mean=5,sd=1)) )
    raneffs[i,] <- rMvdc(1,myMvd)
    b0[i] <- raneffs[i,1]
    b1[i] <- raneffs[i,2]
    res[i] <- raneffs[i,3]
    x[i] <- raneffs[i,4]


    if(abs(b1[i])>.98){
      print("Fixed and random effects have been shrinked to ensure stationarity.
              Check output for updated, shrinked parameter values")
      effVar<- effVar
      effMeans <- .9*effMeans
      shrinked <- TRUE
      i <- 1
    } else {
      shrinked <- FALSE
      i <- i+1
    }
  }





  for(m in 1:n){
   # browser()

    y[1,m] <- rnorm(1,b0[m],sqrt(res[m]))


    # if(b1[m]<(-1)){
    #   b1[m] <- -.99999
    # }else if(b1[m]>(1)){
    #   b1[m] <- .99999
    # }


    for(t in 2:nTime){

      y[t,m] <- b0[m] + b1[m]*(y[t-1,m]-b0[m]) + stats::rnorm(1,0,sqrt(res[m]))

    }
  }
  y <- y[-c(1:100),]
#browser()
  if(max(y)>100){
    shrinked <- TRUE
  }

  res <- list(y=y, xAsOutcome=x, b0=b0,b1=b1,res=res, shrinked=shrinked)
  return(res)

  #}
}

