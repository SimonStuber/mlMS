makeSingleMat <- nimbleFunction(
  run=function(x=double(1), M=integer(0)){
    mat <- array(1,c(M,M))

    diag(mat[1:M,1:M]) <- 0
    ind <- 0
    for(i in 1:M){
      for(j in 1:M){
        if(i!=j){
          ind <- ind +1
          mat[i,j] <- x[ind]
        }
      }
    }

    return(mat[1:M, 1:M])
    returnType(double(2))
  }
)
