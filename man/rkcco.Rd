\name{rkcco}
\alias{rkcco}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Robust kernel cross-covariance opetator
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
# A function
}
\usage{
rkcco(X, Y, lossfu = "Huber", kernel = "rbfdot", gamma = 1e-05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
a  data matrix  index by row
}
  \item{Y}{
%%     ~~Describe \code{Y} here~~
a  data matrix  index by row
}
  \item{lossfu}{
%%     ~~Describe \code{lossfu} here~~
a loss function: square,  Hampel's or Huber's loss
}
  \item{kernel}{
%%     ~~Describe \code{kernel} here~~
a positive definite kernel
}
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
the  hyper-parameters
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{rkcmx }{ Robust kernel center matrix of  X dataset}
 \item{rkcmy }{ Robust kernel center matrix of  Y dataset}
 \item{rkcmx }{ Robust kernel covariacne operator of  X dataset}
 \item{rkcmy }{ Robust kernel covariacne operator of  Y dataset}
 \item{rkcmx }{ Robust kernel cross-covariacne operator of  X  and Y datasets}
}
\author{
%%  ~~who you are~~
Md Ashad Alam
<malam@tulane.edu>
}

\references{
%% ~put references to the literature/web site here ~
{Md Ashad Alam, Kenji Fukumizu and  Yu-Ping Wang (2018),}
{Influence Function and Robust Variant of Kernel Canonical Correlation Analysis,}
{Neurocomputing, Vol. 304 (2018) 12-29.}

 {M. Romanazzi (1992),}
{ Influence in canonical correlation analysis,}
{ Psychometrika
 vol 57(2) (1992) 237-259.}
 }


\seealso{
See also as  \code{\link{rkcca}} \code{\link{snpfmridata}}, \code{\link{ifrkcca}}
}
\examples{

##Dummy data:

X <- matrix(rnorm(2000),200); Y <- matrix(rnorm(2000),200)

rkcco(X,Y, "Huber","rbfdot", 1e-05)
}
