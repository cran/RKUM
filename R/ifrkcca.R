ifrkcca <-
function (X, Y, lossfu ="Huber", kernel="rbfdot", gamma=0.00001, ncomps=10, jth=1) # 10^-5
 {
    n <- dim(X)[1]
    RXY <- rkcco(X, Y, lossfu, kernel, gamma)
    CKX <- RXY$rkcmx
    CKY <- RXY$rkcmx
    KXX <- RXY$rkcooxx
    KYY <- RXY$rkcooyy
    KXY <- RXY$rkcooxy
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
    cor<- as.vector(abs(round(as.double(ei$gvalues[1]))))
    xcoef <- matrix(as.double(ei$gvectors[1:n, 1:ncomps]),n)
    ycoef<- matrix(as.double(ei$gvectors[(n + 1):(2 * n),
       1:ncomps]), n)
    ss <-ei$hfinv
    CVX <- CKX%*%xcoef
    CVY <- CKY%*%ycoef
    IFcor1 <-abs(2*cor*(CVX*CVY) - cor^2*(CVX^2) - cor^2*(CVY^2))[,jth]
    IFcor <- as.vector(IFcor1/(sqrt(sum(IFcor1^2))))
    IXY <- gmi(CD) #inverse of XX and YY
    LX <- (ss[1:n,1:n]%*%gmi(ss[1:n,1:n]%*%AB[0:n,(n + 1):(2 * n)]%*%IXY[(n + 1):(2 * n),(n + 1):(2 * n)]%*%AB[0:n,(n + 1):(2 * n)]%*%ss[1:n,1:n] - diag(rep(cor[1]^2,n)))%*%ss[1:n,1:n])
    LY <- (ss[(n + 1):(2 * n),(n + 1):(2 * n)]%*%(gmi(ss[(n + 1):(2 * n),(n + 1):(2 * n)]%*%AB[(n + 1):(2 * n),(n + 1):(2 * n)]%*%IXY[1:n, 1:n]%*%AB[(n + 1):(2 * n),(n + 1):(2 * n)]%*%ss[(n + 1):(2 * n),(n + 1):(2 * n)] - diag(rep(cor[1]^2,n))))%*%ss[(n + 1):(2 * n),(n + 1):(2 * n)])
    LX <- (LX +t(LX))/2
    LY <- (LY +t(LY))/2
    P1X <-  (LX%*%CVX[,1])%*%t(cor[1]*(CVY[,1]-cor[1]*CVX[,1]))
    P2X <- (LX%*% AB[0:n,(n + 1):(2 * n)]%*%IXY[(n + 1):(2 * n),(n + 1):(2 * n)]%*%CVX[,1])%*%(CVX[,1]-cor[1]*CVY[,1])
    P3X <- 1/2*(1-CVX[,1]^2)*xcoef[,1]
    IFX <- P3X-P1X-P2X
    P1Y <-  (LY%*%CVY[,1])%*%t(cor[1]*(CVX[,1]-cor[1]*CVY[,1]))
    P2Y <- (LY%*% AB[(n + 1):(2 * n),(n + 1):(2 * n)]%*%IXY[1:n,1:n]%*%CVY[,1])%*%(CVY[,1]-cor[1]*CVX[,1])
    P3Y <- 1/2*(1-CVY[,1]^2)*ycoef[,1]
    IFY <- P3Y-P1Y-P2Y
    return(list(ifrkcor=IFcor, ifrkxcv=IFX, ifrkycv=IFY))
 }
