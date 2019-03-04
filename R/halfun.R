halfun <-
function(u)
  {
      qn<-  quantile(u, c(.5, .75, .85))
      c1<-qn[1]
      c2 <- qn[2]
      c3 <- qn[3]
U <- abs(u)
A <- (U <= c1)
B <- (U > c1) & (U <= c2)
C <- (U > c2) & (U <= c3)
D <- (U > c3)
rho <- U
rho[A] <- (U[A] * U[A])/2
rho[B] <- c1 * (U[B] - c1/2)
rho[C] <- c1 * (c2 - c1/2) + c1 * (U[C] - c2) * (1 - (U[C] - c2)/(c3 - c2)/2)
rho[D] <- (c1 * (c2 - c1 + c3))/2
rho
  }
