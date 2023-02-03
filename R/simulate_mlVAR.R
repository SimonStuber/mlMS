library(nimble)

simulate_mlVAR <- function(n, nTime,nItems,
                         b1.means, b.means,b1.sd,b.sd,
                         dfResCov, scaleResCov){

 # browser()

  nTime <- nTime +100

  #M <- length(probs1)
  # nVars <- dim(man.inv.var)[1]
  #  cholMat <- chol(man.inv.var)
  #v <- mvtnorm::rmvnorm(n, v.means, v.cov)
  #vmat <- makeMat(v, n, M)




 # ps <- array(NA, c(nTime,M))
  y <- array(NA, c(nTime, n, nItems))
  #s <- array(NA, nTime)
  b1 <- array(NA, c(n,nItems))
  #b <- array(NA, c(n,nItems,nItems,M))
  resCov <- array(NA, c(n,nItems,nItems))
  bRaneff <- array(NA, c(n,nItems*nItems))
  b <- array(NA, c(n, nItems, nItems))
 # scaleResCov <- scaleResCov


  #for(m in 1:M){
    for(i in 1:n){
     # scaleResCov[,,m] <- scaleResCov[,,m]*dfResCov[m]
      resCov[i,,] <- LaplacesDemon::rinvwishart(dfResCov, scaleResCov)
   # }

    b1 <- mvtnorm::rmvnorm(n = n,b1.means, b1.sd)

    i <- 1
    while(i <= n) {
      bRaneff[i,] <- mvtnorm::rmvnorm(1,b.means, b.sd)
      b[i,,] <- bRaneff[i,]
      if((any((Re(eigen(b[i,,])$values)^2 + Im(eigen(b[i,,])$values)^2)>1))){
        print("Parameters have been shrinked")
        b.sd <- .9*b.sd
        b.means <- .9*b.means
        i <- 1
      } else {
      i <- i+1
      }
    }



  }



  #tp <- vmat
  for(m in 1:n){


   # tp[m,,] <- exp(vmat[m,,])/rowSums(exp(vmat[m,,]))
  #  ps[1,] <- probs1
   # s[1] <- sample(1:M,1,FALSE,ps[1,])
    y[1,m, ] <- b1[m,] + mvtnorm::rmvnorm(1, rep(0,nItems), (resCov[m,,]))

    #  p <- matrix(NA, 1:nTime,M)
    #  p[1,] <- probs1%*%tp
    #  for(t in 2:nTime){
    #    p[t,] <- p[t-1,]%*%tp
    #  }
    # for(state in 1:M){
    #   if(b[m,state]<(-1)){
    #     b[m,state] <- -.99999
    #   }else if(b[m,state]>(1)){
    #     b[m,state] <- .99999
    #   }
    # }

    for(t in 2:nTime){
      #for(state in 1:M){
       # ps[t,state] <- tp[m,s[t-1],state]
      #}
     # s[t] <- sample(1:M,1,FALSE,ps[t,])
     # browser()

      # if(s[t-1]==s[t]){
      y[t,m, ] <- t(b1[m,] + b[m,,]%*%(y[t-1,m,]-b1[m,])) + mvtnorm::rmvnorm(1, rep(0,nItems), resCov[m,,])
      #}else{
      #  y[t,m] <- b1[s[t]] + rnorm(1,0,y.hat.err)
      #}

    }
  }

  y <- y[-c(1:100),,]

  res <- list(y=y,b1=b1,b=b,res=resCov, bRaneff=bRaneff, b.means=b.means, b.sd=b.sd)
  return(res)
}



makeMat <- nimbleFunction(
  run=function(x=double(2), N=integer(0), M=integer(0)){
    mat <- array(1,c(N,M,M))
    for(n in 1:N){
      diag(mat[n,1:M,1:M]) <- 0
      ind <- 0
      for(i in 1:M){
        for(j in 1:M){
          if(i!=j){
            ind <- ind +1
            mat[n,i,j] <- x[n,ind]
          }
        }
      }
    }
    return(mat[1:N, 1:M, 1:M])
    returnType(double(3))
  }
)

