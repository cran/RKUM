ifmkcca <-
function(xx, yy, zz, kernel="rbfdot",gamma=0.00001, ncomps=1,  jth=1)
  {
  CMCOVKX <- rkcco(xx, yy, "square", kernel, gamma)
  CMCOVKZ <- rkcco(xx, zz, "square", kernel, gamma)
  CKX <- CMCOVKX$rkcmx
  CKY <- CMCOVKX$rkcmy
  CKZ <- CMCOVKZ$rkcmy
  KXX <- CMCOVKX$rkcooxx
  KYY <- CMCOVKX$rkcooyy
  KZZ <- CMCOVKZ$rkcooyy
  n <- dim(xx)[1]
  W <- rep(1/n, n)
  KXY <- tcrossprod(CKX%*%diag(W),CKY)
  KXZ <- tcrossprod(CKX%*%diag(W),CKZ)
  KYZ <- tcrossprod(CKY%*%diag(W),CKZ)
  m <- 3
  AB <- matrix(0, n * m, n * m)
  CD <- AB
  AB[0:n,0:n ] <- 0
  AB[(n + 1):(2 * n),(n + 1):(2 * n) ] <-0
  AB[(2 * n + 1):(3 * n),(2 * n + 1):(3 * n) ] <-0
  AB[0:n,(n + 1):(2 * n)] <- (KXY)
  AB[0:n,(2 * n + 1):(3 * n)] <- (KXZ)
  AB[(n+1):(2 * n),(2 * n + 1):(3 * n)] <- (KYZ)
  AB <- AB + t(AB)
  CD[1:n, 1:n] <- (KXX+ diag(rep(1e-6,n)))
  CD[(n + 1):(2 * n), (n + 1):(2 * n)] <- (KYY+ diag(rep(1e-6,n)))
  CD[(2 * n + 1):(3 * n), (2 * n + 1):(3 * n)] <- (KZZ+ diag(rep(1e-6,n)))
  CD <- (CD + t(CD))/2
  ei <-  gmedc(AB, CD)
  kcor <- as.double(ei$gvalues[1:ncomps])
  xcoef <- matrix(as.double(ei$gvectors[1:n, 1:ncomps]),   n)
  ycoef<- matrix(as.double(ei$gvectors[(n + 1):(2 * n), 1:ncomps]), n)
  zcoef<- matrix(as.double(ei$gvectors[(2 * n + 1):(3 * n), 1:ncomps]), n)
 CVX <- CKX%*%xcoef
 CVY <- CKY%*%ycoef
 CVZ <- CKZ%*%zcoef
 IFcor <- abs(2*kcor*(CVX*CVY) + 2*kcor*(CVX*CVZ)+2*kcor*(CVZ*CVY)- kcor^2*(CVX^2) - kcor^2*(CVY^2)-kcor^2*(CVZ^2))[,jth]
 return(list(ifmkcor=IFcor))
 }
