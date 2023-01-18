lse <- nimbleFunction(
  run = function(x=double(1)){
  c <- max(x)
  res <- c + log(sum(exp(x-c)))
  return(res)
  returnType(double(0))
  }
)
