lkm <-
function(X)
   {
     K <- X%*%t(X)
     return(K)
   }
