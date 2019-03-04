snpfmridata <-
function(n=300, gamma=0.00001, ncomps=2, jth=1)
   {
      sig1 <- 0.5
sig2 <- 1
p1 <- 100
p2 <- 100
      pt <- 100
a <- rep(0, pt)
      a[1:p1] <- ranuf(p1)
b <- rep(0,pt)
b[1:p2] <- ranuf(p2) #using rndu.r file
maf  <-  runif(pt,0.2, 0.4) # (0,1)/5+0.2  # form uiform distribution
mu <- rnorm(n, 0, 1)
      mu <- sign(mu)*(abs(mu)+0.1);# normal distribution add some samll term 0.1
mu <-mu*sig1
      ind <- rep(F,pt)        #logical(pt) # logical zeors
      ind <- as.integer(ind)
ind[1:p2] <- 1 # take 1 1: p2  matrix( rnorm(N*M,mean=0,sd=1), N, M)
m1 <- mu%*%t(a)+matrix(rnorm(n*pt, 0, 1), n, pt)*sig1 #fmri
m2 <- mu%*%t(b)+matrix(rnorm(n*pt, 0, 1), n, pt)*sig2 # SNP
 for ( i in 1: pt)
 {
    p <- maf[i]
     bias <- log(1/p-1)
    if (ind[i]==TRUE)
  {
    gt <- m2[,i]
    gt <- 1/(1+exp(-gt+bias)) # for correlated variables case
        }
     else
     gt <- 1/(1+exp(bias))      #For no correlated variable #control
     temp <- runif(n,0,1)
  m2[,i] <- gt > temp
   temp <- runif(1, 0,1)
    m2[,i] <- m2[,i]+(gt > temp)
}
  E3X <- m1
  E3Y <- m2
n1 <- (0.05*n)
      ss <- sample(1:n, n1)
 m1 <-mu[ss]%*%t(a)+matrix(rnorm(n1*pt, 0, 1), n1,pt)*10 #fmri
     m2 <- mu[ss]%*%t(b)+matrix(rnorm(n1*pt, 0, 1), n1,pt)*20
       for ( i in 1: pt)
 {
    p <- maf[i]
     bias <- log(1/p-1)
    if (ind[i]==TRUE)
  {
    gt <- m2[,i]
    gt <- 1/(1+exp(-gt+bias)) # for correlated variables
        }
     else
     gt <- 1/(1+exp(bias))      #For no correlated variable
     temp <- runif(n1,0,1)
  m2[,i] <- gt > temp
   temp <- runif(1, 0,1)
    m2[,i] <- m2[,i]+(gt > temp)
}

E3Xot <- rbind(E3X[-ss,],m1)
  E3Yot <- rbind(E3Y[-ss,],m2)
      CCAID <- ifcca(E3X,E3Y, gamma, ncomps, jth)
      CCACD <- ifcca(E3Xot,E3Yot, gamma, ncomps, jth)
      KCCAID <- ifrkcca(E3X,E3Y, "square", "linear", gamma, ncomps, jth)
      KCCACD <- ifrkcca(E3Xot,E3Yot, "square", "linear", gamma, ncomps, jth)
      HARKCCAID <- ifrkcca(E3X, E3Y, "Hampel", "linear", gamma, ncomps, jth)
      HARKCCACD <- ifrkcca(E3Xot,E3Yot, "Hampel", "linear",gamma, ncomps, jth)
      HURKCCAID <- ifrkcca(E3X,E3Y, "Huber", "linear", gamma, ncomps, jth)
      HURKCCACD <- ifrkcca(E3Xot, E3Yot, "Huber", "linear",gamma, ncomps, jth)
      CCAIFID <-  CCAID$iflccor
      CCAIFCD <- CCACD$iflccor
      KCCAIFID <- KCCAID$ifrkcor
      KCCAIFCD <- KCCACD$ifrkcor
      HAKCCAIFID <- HARKCCAID$ifrkcor
      HAKCCAIFCD <- HARKCCACD$ifrkcor
      HUKCCAIFID <- HURKCCAID$ifrkcor
      HUKCCAIFCD <- HURKCCAID$ifrkcor
   {
      par(mfrow=c(4,2))
par(mar=c(5.1,5.1,4.1,2.1))
      plot(CCAIFID, pch=21:22,main=list("Ideal Data",cex =2, font=15), sub=list("Linear CCA",cex =2, font=15),ylab="",xlab=list(" ",cex =2, font=15),
     type="o",lwd=3,cex.axis=2,cex.lab=2,axes=FALSE)
     axis(2, las=0,cex.axis=2,cex.lab=2)
      plot(CCAIFCD, pch=21:22, main=list("Contaminated Data",cex =2, font=15),sub=list("Linear CCA",cex =2, font=15), ylab="",xlab=list(" ",cex =2, font=15),
     type="o",lwd=3,cex.axis=2,cex.lab=2,axes=FALSE)
     axis(2, las=0,cex.axis=2,cex.lab=2)
      plot(KCCAIFID, pch=21:22, main=list("",cex =2, font=15),sub=list("Kernel CCA",cex =2, font=15),ylab="",xlab=list(" ",cex =2, font=15),
     type="o",lwd=3,cex.axis=2,cex.lab=2,axes=FALSE)
     axis(2, las=0,cex.axis=2,cex.lab=2)
      plot(KCCAIFCD, pch=21:22, main=list("",cex =2, font=15),sub=list(" Kernel CCA",cex =2, font=15),ylab="",xlab=list(" ",cex =2, font=15),
     type="o",lwd=3,cex.axis=2,cex.lab=2,axes=FALSE)
     axis(2, las=0,cex.axis=2,cex.lab=2)
      plot(HAKCCAIFID, pch=21:22, main=list("",cex =2, font=15), sub=list(" Hampel's  robust kernel CCA",cex =2, font=15),ylab="",xlab=list(" ",cex =2, font=15),
     type="o",lwd=3,cex.axis=2,cex.lab=2,axes=FALSE)
     axis(2, las=0,cex.axis=2,cex.lab=2)
      plot(HAKCCAIFCD, pch=21:22, main=list("",cex =2, font=15),sub=list(" Hampel's robust kernel CCA",cex =2, font=15), ylab="",xlab=list(" ",cex =2, font=15),
     type="o",lwd=3,cex.axis=2,cex.lab=2,axes=FALSE)
     axis(2, las=0,cex.axis=2,cex.lab=2)
      plot(HUKCCAIFID, pch=21:22, main=list("",cex =2, font=15),sub=list(" Huber's robust kernel CCA",cex =2, font=15),ylab="",xlab=list(" ",cex =2, font=15),
     type="o",lwd=3,cex.axis=2,cex.lab=2,axes=FALSE)
     axis(2, las=0,cex.axis=2,cex.lab=2)
      plot(HUKCCAIFCD, pch=21:22, main=list("",cex =2, font=15),sub=list(" Huber's robust kernel CCA",cex =2, font=15),ylab="",xlab=list(" ",cex =2, font=15),
     type="o",lwd=3,cex.axis=2,cex.lab=2,axes=FALSE)
     axis(2, las=0,cex.axis=2,cex.lab=2)
     }
 return(list(IFCCAID = CCAIFID,  IFCCACD = CCAIFCD,  IFKCCAID = KCCAIFID,  IFKCCACD = KCCAIFCD, IFHAKCCAID = HAKCCAIFID,  IFHAKCCACD = HAKCCAIFCD, IFHUKCCAID = HUKCCAIFID,  IFHUKCCACD = HUKCCAIFCD))
 }
