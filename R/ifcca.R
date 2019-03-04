ifcca <-
function (X, Y, gamma=0.00001, ncomps=2, jth=1)
    {
    Xnames <- dimnames(X)[[2]]
    Ynames <- dimnames(Y)[[2]]
    ind.names <- dimnames(X)[[1]]
    Cxx <- var(X, na.rm = TRUE, use = "pairwise") + diag(gamma,
        ncol(X))
    Cyy <- var(Y, na.rm = TRUE, use = "pairwise") + diag(gamma,
        ncol(Y))
    Cxy <- cov(X, Y, use = "pairwise")
    res <- gm3edc(Cxy, Cxx, Cyy)
    names(res) <- c("cor", "xcoef", "ycoef")
    cor = res$cor[1:ncomps]
    scores <- lcv(X, Y, res)
    CVX <- scores$xscores[,1: ncomps]
    CVY <- scores$yscores[,1: ncomps]
     IFcor1 <- abs(2*cor*(CVX*CVY) - cor^2*(CVX^2) - cor^2*(CVY^2))[,jth]
     IFcor <- as.vector(IFcor1/(sqrt(sum(IFcor1^2))))
     return(list(iflccor=IFcor))
   }
