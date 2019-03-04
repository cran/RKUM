hulfun <-
function(x)
  {
    n<-length(x)
    a <-median(x)
    y <- rep(0, n)
    for (i in 1: n)
    {
      if( x[i] <= (a))
      {
      y[i] <-x[i]^2/2
      }
      else
       y[i] <- a*x[i]-a^2/2
      }
    return(y)
   }
