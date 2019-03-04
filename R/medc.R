medc <-
function(A,fn=sqrt) 
 {
   e<-eigen(A)
   y<-e$vectors
   v<-e$values
   return(tcrossprod(y%*%diag(fn(v)),y))
 }
