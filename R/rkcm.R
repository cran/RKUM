rkcm <-
function(X, lossfu ="Huber", kernel="rbfdot")     
   {
    n <-dim(X)[1]
  if (lossfu=="square")
    {
      if(kernel=="linear")
     {
     K <- lkm(X)
     Id <- diag(1, nrow = n)
     Id1 <- matrix(1, nrow = n, ncol = n)
     H <- Id - Id1/n
     K.robust.M <- tcrossprod(H%*%K, H)
      }
   if(kernel=="rbfdot")
     {
     K <- gkm(X)
    Id <- diag(1, nrow = n)
    Id1 <- matrix(1, nrow = n, ncol = n)
    H <- Id - Id1/n
    K.robust.M <- tcrossprod(H%*%K, H)
     }
   if(kernel=="ibs")
     {
     K <- ibskm(X)
    Id <- diag(1, nrow = n)
    Id1 <- matrix(1, nrow = n, ncol = n)
    H <- Id - Id1/n
    K.robust.M <- tcrossprod(H%*%K, H)
    }
  }
   if(lossfu=="Hampel")
   {
    if(kernel=="linear")
     {
     K <- lkm(X)
     W <- rep(1/n, n)
    aa<- sqrt(abs(diag(K)-2*as.vector(crossprod(W,K))+ as.vector(crossprod(W,K%*%W))))
    for (k in 1:100)
    {
       Obj.old<- halofun(aa)
        HV <-hadr(aa)
        THV <- sum(HV)
        W <-HV/THV
  R_mean_E <-as.vector(crossprod(W,K))
  aa<- sqrt(abs(diag(K)-2*as.vector(crossprod(W,K))+ as.vector(crossprod(W,K%*%W))))
        Obj.new<- halofun(aa)
        Stop <- abs(Obj.old-Obj.new)/Obj.old
        if( Stop < 0.1^8)
          {
            break
          }
     }
    ee<- rep(1, n)
    H <- diag(1, n)- as.vector(tcrossprod(ee,W))
    K.robust.M <-tcrossprod(H%*%K,H)
     }
    if(kernel=="rbfdot")
     {
     K <- gkm(X) 
   W <- rep(1/n, n)
    aa<- sqrt(abs(diag(K)-2*as.vector(crossprod(W,K))+ as.vector(crossprod(W,K%*%W))))
    for (k in 1:100)
    {
       Obj.old<- halofun(aa)
        HV <-hadr(aa)
        THV <- sum(HV)
        W <-HV/THV
  R_mean_E <-as.vector(crossprod(W,K))
  aa<- sqrt(abs(diag(K)-2*as.vector(crossprod(W,K))+ as.vector(crossprod(W,K%*%W))))
        Obj.new<- halofun(aa)
        Stop <- abs(Obj.old-Obj.new)/Obj.old
        if( Stop < 0.1^8)
          {
            break
          }
     }
    ee<- rep(1, n)
    H <- diag(1, n)- as.vector(tcrossprod(ee,W))
    K.robust.M <-tcrossprod(H%*%K,H)
    }
   if(kernel=="ibs")
     {
     K <- ibskm(X)
    W <- rep(1/n, n)
    aa<- sqrt(abs(diag(K)-2*as.vector(crossprod(W,K))+ as.vector(crossprod(W,K%*%W))))
    for (k in 1:100)
    {
       Obj.old<- halofun(aa)
        HV <-hadr(aa)
        THV <- sum(HV)
        W <-HV/THV
  R_mean_E <-as.vector(crossprod(W,K))
  aa<- sqrt(abs(diag(K)-2*as.vector(crossprod(W,K))+ as.vector(crossprod(W,K%*%W))))
        Obj.new<- halofun(aa)
        Stop <- abs(Obj.old-Obj.new)/Obj.old
        if( Stop < 0.1^8)
          {
            break
          }
     }
    ee<- rep(1, n)
    H <- diag(1, n)- as.vector(tcrossprod(ee,W))
    K.robust.M <-tcrossprod(H%*%K,H)
    }
  }
  if(lossfu=="Huber")
    {
    if (kernel=="linear")
     {
     K <- lkm(X)
     W <- rep(1/n, n)
    aa<- sqrt(diag(K)-2*as.vector(crossprod(W,K))+ as.vector(crossprod(W,K%*%W)))
    for (k in 1:100)
    {
        Obj.old<- hulofun(aa)
        HV <- hudr(aa)
        THV <- sum(HV)
        W <-HV/THV
      R_mean_E <-as.vector(crossprod(W,K))#as.vector(t(W)%*%K)
      aa<- sqrt(diag(K)-2*as.vector(crossprod(W,K))+ as.vector(crossprod(W,K%*%W)))
        Obj.new<- hulofun(aa)
        Stop <- abs(Obj.old-Obj.new)/Obj.old
        if( Stop < 0.1^8)
          {
            break
          }
     }
    ee<- rep(1, n)
    H <- diag(1, n)- as.vector(tcrossprod(ee,W))
    K.robust.M <-tcrossprod(H%*%K,H)
     }
    if(kernel=="rbfdot")
     {
     K <- gkm(X)
     W <- rep(1/n, n)
    aa<- sqrt(diag(K)-2*as.vector(crossprod(W,K))+ as.vector(crossprod(W,K%*%W)))
    for (k in 1:100)
    {
        Obj.old<- hulofun(aa)
        HV <- hudr(aa)
        THV <- sum(HV)
        W <-HV/THV
      R_mean_E <-as.vector(crossprod(W,K))#as.vector(t(W)%*%K)
      aa<- sqrt(diag(K)-2*as.vector(crossprod(W,K))+ as.vector(crossprod(W,K%*%W)))
        Obj.new<- hulofun(aa)
        Stop <- abs(Obj.old-Obj.new)/Obj.old
        if( Stop < 0.1^8)
          {
            break
          }
     }
    ee<- rep(1, n)
    H <- diag(1, n)- as.vector(tcrossprod(ee,W))
    K.robust.M <-tcrossprod(H%*%K,H)
    }
   if(kernel=="ibs")
     {
     K <- ibskm(X)
      W <- rep(1/n, n)
    aa<- sqrt(diag(K)-2*as.vector(crossprod(W,K))+ as.vector(crossprod(W,K%*%W)))
    for (k in 1:100)
    {
        Obj.old<- hulofun(aa)
        HV <- hudr(aa)
        THV <- sum(HV)
        W <-HV/THV
      R_mean_E <-as.vector(crossprod(W,K))#as.vector(t(W)%*%K)
      aa<- sqrt(diag(K)-2*as.vector(crossprod(W,K))+ as.vector(crossprod(W,K%*%W)))
        Obj.new<- hulofun(aa)
        Stop <- abs(Obj.old-Obj.new)/Obj.old
        if( Stop < 0.1^8)
          {
            break
          }
     }
    ee<- rep(1, n)
    H <- diag(1, n)- as.vector(tcrossprod(ee,W))
    K.robust.M <-tcrossprod(H%*%K,H)
   }
  }
   K.robust.M <- (K.robust.M + t(K.robust.M))/2
  return(list(rkcm=K.robust.M)) 
 }
