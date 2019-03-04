ibskm <-
function(Z) 
   {
  n <-  nrow(Z)
    K <- 1 - as.matrix(dist(Z, method = "manhattan") * 0.5/max(1, ncol(Z)))
  return(K)
  }
