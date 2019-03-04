rkcco <-
function(X, Y, lossfu ="Huber", kernel="rbfdot",  gamma= 0.00001) 
 {                      
   n <- nrow(X)
   RTX <- rkcm(X, lossfu, kernel)
   RTY <- rkcm(Y, lossfu, kernel)
   CKX <- RTX$rkcm
   CKY <- RTY$rkcm
   if (lossfu=="square")
     {
      W <- rep(1/n, n)
      RXX1 <- (CKX+diag(rep(gamma,n)))%*%diag(W)%*%CKX 
      RYY1 <- (CKY+diag(rep(gamma,n)))%*%diag(W)%*%CKY 
      RXY1 <- (CKX+diag(rep(gamma,n)))%*%diag(W)%*%CKY 
     }
  if (lossfu=="Hampel")
    {
    W <- rep(1/n, n) 
    error <- sqrt(abs(diag(CKX)*diag(CKY)-2*as.vector(crossprod(W,(CKX*CKY)))+as.vector(crossprod(W,  (CKX*CKY)%*%W))))
    for (k in 1:100)
    {
      Obj.old <- halofun(error)
      HV <- hadr(error)
     THV <- sum(HV)
      W <- HV/THV
      R_SecM <-  tcrossprod(CKX%*%diag(W),CKY)
      error <- sqrt(abs(diag(CKX)*diag(CKY)-2*as.vector(crossprod(W,(CKX*CKY)))+as.vector(crossprod(W,(CKX*CKY)%*%W))))
    Obj.new <-  halofun(error)
      Stop <- abs(Obj.old-Obj.new)/Obj.old
      if( Stop < 0.1^8)
        {
          break
        }
    }
     RXX1 <- (CKX+diag(rep(gamma,n)))%*%diag(W)%*%CKX 
     RYY1 <- (CKY+diag(rep(gamma,n)))%*%diag(W)%*%CKY 
     RXY1 <- (CKX+diag(rep(gamma,n)))%*%diag(W)%*%CKY 
   }
  if (lossfu=="Huber")
   {
   W <- rep(1/n, n) 
   error <- sqrt(abs(diag(CKX)*diag(CKY)-2*as.vector(crossprod(W,(CKX*CKY)))+as.vector(crossprod(W,  (CKX*CKY)%*%W))))
   for (k in 1:100)
    {
      Obj.old <- hulofun(error)
      HV <- hudr(error)
    THV <- sum(HV)
      W <- HV/THV
      R_SecM <-  tcrossprod(CKX%*%diag(W),CKY)
      error <- sqrt(abs(diag(CKX)*diag(CKY)-2*as.vector(crossprod(W,(CKX*CKY)))+as.vector(crossprod(W,(CKX*CKY)%*%W))))
    Obj.new <-  hulofun(error)
      Stop <- abs(Obj.old-Obj.new)/Obj.old
       if( Stop < 0.1^8)
        {
          break
        }
     }
     RXX1 <- (CKX+diag(rep(gamma,n)))%*%diag(W)%*%CKX 
     RYY1 <- (CKY+diag(rep(gamma,n)))%*%diag(W)%*%CKY 
     RXY1 <- (CKX+diag(rep(gamma,n)))%*%diag(W)%*%CKY 
    }
   RXX <- (RXX1 + t(RXX1))/2
   RYY <- (RYY1 + t(RYY1))/2 
   RXY <-RXY1  
  return(list(rkcmx = CKX, rkcmy = CKY, rkcooxx = RXX, rkcooyy=RYY, rkcooxy =RXY)) 
 }
