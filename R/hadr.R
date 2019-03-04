hadr <-
function(u)
  {
      qn<-  quantile(u, c(.5, .75, .85))
      c1<-qn[1]
      c2 <- qn[2]
      c3 <- qn[3]
U <- abs(u)
B <- (U > c1) & (U <= c2)#flat
C <- (U > c2) & (U <= c3)#descending
D <- (U > c3)# zero
psi <- u
psi[B] <- sign(u[B]) * c1
psi[C] <- sign(u[C]) * c1 * (c3 - U[C])/(c3 - c2)
psi[D] <- 0
return(psi/u)
 }
