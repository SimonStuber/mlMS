simulate_randomMLAR <- function(n, nTime,
                         effMeans, effVar,
                         resVarOnX=NULL,
                         arOnX=NULL,
                         meanOnX=NULL,
                         simXAsPredictor=FALSE,
                         simXAsOutcome=FALSE,
                         xOnMean=NULL,
                         xOnResVar=NULL,
                         xOnAr=NULL,
                         xInt=NULL, 
                         L2ResVar=1){


  nTime <- nTime + 100
  y <- array(NA, c(nTime, n))
  xAsPredictor <- c()

    raneffs <- matrix(NA, n, 3)
    b0 <- c()
    b1 <- c()
    res <- c()
    mssd <- c()
    var <- c()

    i <- 1
    shrinked <- FALSE
    while(i <= n) {
     # print(i)
      raneffs[i,] <- mvtnorm::rmvnorm(1,effMeans, effVar)
      if(simXAsPredictor){
        if(any(sapply(list(meanOnX, arOnX, resVarOnX),is.null))){
          print("Missing regression weights to predict random effects by X")
        }
        xAsPredictor[i] <- rnorm(1, 0, 1)
        b0[i] <- raneffs[i,1] + meanOnX*xAsPredictor[i]
        b1[i] <- raneffs[i,2] + arOnX*xAsPredictor[i]
        res[i] <- raneffs[i,3] + resVarOnX*xAsPredictor[i]
      }else{
        b0[i] <- raneffs[i,1]
        b1[i] <- raneffs[i,2]
        res[i] <- raneffs[i,3]
        mssd[i] <- 2*((res[i])/(1-(b1[i]^2)))*(1-b1[i])
        var[i] <- ((res[i])/(1-(b1[i]^2)))
      }

      if(abs(b1[i])>1){
        print("Fixed and random effects have been shrinked to ensure stationarity.
              Check output for updated, shrinked parameter values")
       effVar<- effVar
      effMeans <- .9*effMeans
      shrinked <- TRUE
        i <- 1
      } else {
        i <- i+1
      }
    }
    
    

    if(simXAsOutcome){
      if(any(sapply(list(xInt, xOnMean, xOnAr, xOnResVar),is.null))){
        print("Missing regression weights to predict x by the random effects")
      }
      xAsOutcome <- xInt + xOnMean*b0 + xOnAr*b1 + xOnResVar*res + rnorm(n, 0, sqrt(L2ResVar))
    }


    for(m in 1:n){
      y[1,m] <- b0[m]
      for(t in 2:nTime){
        y[t,m] <- b0[m] + b1[m]*(y[t-1,m]-b0[m]) + rnorm(1,0,sqrt(res[m]))
      }
    }

    y <- y[-c(1:100),]
    
    
    ##### HACKY SIM FIX
    if(any(b1>.99)|any(res<0)){
      shrinked <- TRUE
    }

    out <- list(y=y,
                b0=b0,
                b=b1,
                res=res,
                effMeans=effMeans,
                effVar=effVar,
                shrinked=shrinked)
    if(simXAsPredictor){
      out <-  list(y=y,
                   b0=b0,
                   b=b1,
                   res=res,
                   effMeans=effMeans,
                   effVar=effVar,
                   shrinked=shrinked,
                   resVarOnX=resVarOnX,
      arOnX=arOnX,
      meanOnX=meanOnX,
      xAsPredictor=xAsPredictor)
    }
    if(simXAsOutcome){
     out <-  list(y=y,
                   b0=b0,
                   b=b1,
                   res=res,
                   effMeans=effMeans,
                   effVar=effVar,
                   shrinked=shrinked,
                   xAsOutcome=xAsOutcome,
                   xOnMean=xOnMean,
                   xOnAr=xOnAr,
                   xOnResVar=xOnResVar,
                   xInt=xInt)
    }
    return(out)

  #}
}

