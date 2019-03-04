gkm <-
function(X)
   {
     n <- dim(X)[1]
     SGX <-  mdbw(X)
     sx2  <- 2*SGX*SGX
     ab <- X%*%t(X) 
     aa <- as.matrix(diag(ab))
     D <- matrix(rep(aa,each=n), ncol=n, byrow=TRUE) 
     xx <-pmax(D + t(D) - 2*ab,  mat.or.vec(n, n))
     K <- exp(-xx/sx2)  
     return(K)
   }
