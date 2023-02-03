library(nimble)

simulate_mlHMM <- function(probs1,  n, nTime,nItems,
                         b1.means, b1.sd,
                         v.means, v.cov, dfResCov, scaleResCov){

 #browser()

  nTime <- nTime +100

  M <- length(probs1)
  # nVars <- dim(man.inv.var)[1]
  #  cholMat <- chol(man.inv.var)
  v <- mvtnorm::rmvnorm(n, v.means, v.cov)
  vmat <- makeMat(v, n, M)




  ps <- array(NA, c(nTime,M))
  y <- array(NA, c(nTime, n, nItems))
  s <- array(NA, nTime)
  b1 <- array(NA, c(n,nItems,M))
  #b <- array(NA, c(n,nItems,nItems,M))
  resCov <- array(NA, c(n,nItems,nItems,M))
  bRaneff <- array(NA, c(n,nItems*nItems,M))
#  b <- array(NA, c(n, nItems, nItems,M))
 # scaleResCov <- scaleResCov


  for(m in 1:M){
    for(i in 1:n){
     # scaleResCov[,,m] <- scaleResCov[,,m]*dfResCov[m]
      resCov[i,,,m] <- LaplacesDemon::rinvwishart(dfResCov[m], scaleResCov[,,m])
    }

    b1[,,m] <- mvtnorm::rmvnorm(n = n,b1.means[,m], b1.sd[,,m])

   # i <- 1
#    while(i <= n) {
 #     bRaneff[,,m] <- mvtnorm::rmvnorm(1,b.means[,m], b.sd[,,m])
 #     b[i,,,m] <- bRaneff[i,,m]
#      if((any((Re(eigen(b[i,,,m])$values)^2 + Im(eigen(b[i,,,m])$values)^2)>.9))){
#        print("Parameters have been shrinked")
#        b.sd[,,m] <- .9*b.sd[,,m]
#        b.means[,m] <- .9*b.means[,m]
#        i <- 1
#      } else {
#      i <- i+1
#      }
#    }



  }



  tp <- vmat
  for(m in 1:n){


    tp[m,,] <- exp(vmat[m,,])/rowSums(exp(vmat[m,,]))
    ps[1,] <- probs1
    s[1] <- sample(1:M,1,FALSE,ps[1,])
    y[1,m, ] <- b1[m,,s[1]] + mvtnorm::rmvnorm(1, rep(0,nItems), (resCov[m,,,s[1]]))

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
      for(state in 1:M){
        ps[t,state] <- tp[m,s[t-1],state]
      }
      s[t] <- sample(1:M,1,FALSE,ps[t,])
     # browser()

      # if(s[t-1]==s[t]){
      y[t,m, ] <- t(b1[m,,s[t]]) + mvtnorm::rmvnorm(1, rep(0,nItems), resCov[m,,,s[t]])
      #}else{
      #  y[t,m] <- b1[s[t]] + rnorm(1,0,y.hat.err)
      #}

    }
  }

  y <- y[-c(1:100),,]

  res <- list(y=y,b1=b1,res=resCov, tp=tp, v=v, raneffs=v, b.sd=b.sd)
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

