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

