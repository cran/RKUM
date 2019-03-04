hudr <-
function(x)
   {
    n<-length(x)
    a <- median(x)
    y <- rep(0, n)
    for (i in 1: n)
    {
       if( x[i] <= (a))
       {
        y[i] <-1
        }
        else
        y[i] <- a/x[i]
    }
  return(y)
 }
