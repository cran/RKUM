snpfmrimth3D <-
function(n=500, gamma=0.00001, ncomps=1, jth=1)
     {
      sig1 <- 0.5 #signal level
      sig2 <- 1 #noise level
      sig3 <- 0.1 #noise level
    p1 <- 100
    p2 <- 100
      p3 <- 100
      pt <- 100 # dimension for correlated voxels,correlated SNPs and data.
     a <- rep(0, pt)
      a[1:p1] <-  ranuf(p1)
     b <- rep(0,pt)
    b[1:p2] <-  ranuf(p2) #using rndu.r file
      maf  <-  runif(pt,0.2, 0.4) # (0,1)/5+0.2  # form uiform distribution
     mu <- rnorm(n, 0, 1)
      mu <- sign(mu)*(abs(mu)+0.1);# normal distribution add some samll term 0.1
    mu <-mu*sig1 # Common feature
      ind <- rep(F,pt)        #logical(pt) # logical zeors
      ind <- as.integer(ind)
    ind[1:p2] <- 1 # take 1 1: p2  matrix( rnorm(N*M,mean=0,sd=1), N, M)
    m1 <- mu%*%t(a)+matrix(rnorm(n*pt, 0, 1), n,pt)*sig1 #fmri
    m2 <- mu%*%t(b)+matrix(rnorm(n*pt, 0, 1), n,pt)*sig2 # SNP
         # Methlation
      prop <- c(0.20,0.30,0.27,0.23)
      effect <- 2.5
      n.sample=n
      cluster.sample.prop = c(0.30,0.30,0.40)
      n.cluster <- length(cluster.sample.prop)
      delta.methyl=effect
      p.DMP=0.2
      sigma.methyl=NULL
      DMP <- sapply(1:n.cluster, function(x) rbinom(p3, 1, prob = p.DMP))
      d <- lapply(1:n.cluster,function(i) {
  effect <- DMP[,i]*delta.methyl
   mvnod(n=cluster.sample.prop[i]*n.sample, mu=effect, Sigma=diag(1, p3,p3))})
      sim.methyl <- do.call(rbind,d)
      m3 <- rlogit(sim.methyl) + matrix(rnorm(n*p3, 0, 1), n, p3)*sig3
     #tO CONSTRACT Normal distribution o binomial distribution
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
     temp <- runif(n,0,1)
  m2[,i] <- gt > temp
   temp <- runif(1, 0,1)
    m2[,i] <- m2[,i]+(gt > temp)
  }
  EX <- m1
    EY <- m2
      EZ <- m3
n1 <- (0.05*n)
      ss <- sample(1:n, n1)
    m1 <-mu[ss]%*%t(a)+matrix(rnorm(n1*pt, 0, 1), n1,pt)*1.5 #fmri
     m2 <- mu[ss]%*%t(b)+matrix(rnorm(n1*pt, 0, 1), n1,pt)*2
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
 EXout <- rbind(EX[-ss,],m1)
EYout <- rbind(EY[-ss,],m2)
      EZout <- rbind(EZ[-ss,], rlogit(sim.methyl[ss,])+ matrix(rnorm(n1*p3, 0, 1), n1, p3)*3)
      ReIF1 <- ifmkcca(EX, EY,  EZ, "linear", gamma, ncomps, jth)
      ReIFout1 <- ifmkcca(EXout, EYout, EZout, "linear", gamma, ncomps,jth)
     ReIF <-  ReIF1$ifmkcor
      ReIFout <-  ReIFout1$ifmkcor
       {
   par(mfrow=c(1,2))
   par(mar=c(5.1,5.1,4.1,2.1))
      plot(ReIF, pch=21:22, main=list("Ideal data",cex =2, font=15),sub=list("Multiple kernel CCA",cex =2, font=15), ylab="", xlab=list(" ",cex =2, font=15),
     type="o",lwd=3,cex.axis=2,cex.lab=2,axes=FALSE)
     axis(2, las=0,cex.axis=2,cex.lab=2)
      plot( ReIFout, pch=21:22, main=list("Contaminated data", cex =2, font=15), sub=list(" Multiple kernel CCA",cex =2, font=15), ylab="",xlab=list(" ",cex =2, font=15),
     type="o",lwd=3,cex.axis=2,cex.lab=2,axes=FALSE)
     axis(2, las=0,cex.axis=2,cex.lab=2)
      }
 return(list(IFim=  ReIF, IFcm =  ReIFout))
 }
