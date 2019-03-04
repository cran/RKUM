mdbw <-
function(X)
  {
   n <- dim(X)[1]
   ab <- X%*%t(X)
   aa <- as.matrix(diag(ab))
   Dx <-matrix(rep(aa,each=n), ncol=n, byrow=TRUE) +  matrix(rep(aa,each=n),nrow=n) - 2*ab
   Dx <- Dx-diag(diag(Dx))
   dx <- matrix(Dx, n*n,1)
   s <- sqrt(median(dx[dx!=0]))
  return (s)
 }
