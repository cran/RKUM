ranuf <-
function(p)
    {
   y <- runif(p,0,1)*3 - 1.5
  y <- sign(y)*(abs(y) + 0.5)
    return(y)
   }
