rkcca <-
function (X,  Y, lossfu ="Huber", kernel="rbfdot", gamma=0.00001, ncomps=10) 
    {
    RXY <- rkcco(X, Y, lossfu, kernel, gamma)    
    CKX <- RXY$rkcmx
    CKY <- RXY$rkcmy
    KXX <- RXY$rkcooxx
    KYY <- RXY$rkcooyy
    KXY <- RXY$rkcooxy
    n <- dim(KXX)[1]
    m <- 2
    AB <- matrix(0, n *2, 2*n)
    AB[0:n,0:n ] <- 0
    AB[0:n,(n + 1):(2 * n)] <- (KXY)
    AB[(n + 1):(2 * n),0:n ] <- t(KXY)
    AB[(n + 1):(2 * n),(n + 1):(2 * n) ] <-0
    CD <- matrix(0, n * m, n * m)
    CD[1:n, 1:n] <- (KXX+diag(rep(1e-6,n))) 
    CD[(n + 1):(2 * n), (n + 1):(2 * n)] <- (KYY+diag(rep(1e-6,n)))
    CD <- (CD + t(CD))/2
    ei <- gmedc(AB, CD)
    kcor<- abs(as.double(ei$gvalues[1:ncomps]))
    xcoef <- matrix(as.double(ei$gvectors[1:n, 1:ncomps]),n)
    ycoef<- matrix(as.double(ei$gvectors[(n + 1):(2 * n), 
        1:ncomps]), n)
    CVX <- CKX%*%xcoef
    CVY <- CKY%*%ycoef
     return(list(rkcor = kcor, rxcoef = xcoef, rycoef = ycoef,  rxcv= CVX, rycv=CVY))
 }
