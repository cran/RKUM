gmedc <-
function(A,B=diag(nrow(A)))
  {
   ss<-medc(B,function(x) udtd(sqrt(x)))
   ev <- eigen(ss%*%A%*%ss)
   return(list(hfinv=ss,gvalues=ev$values,gvectors=ss%*%ev$vectors))
 }
