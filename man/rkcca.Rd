\name{rkcca}
\alias{rkcca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
 Robust kernel canonical correlation analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~

#A robust correlation
}
\usage{
rkcca(X, Y, lossfu = "Huber", kernel = "rbfdot", gamma = 1e-05, ncomps = 10)
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
  \item{ncomps}{
%%     ~~Describe \code{ncomps} here~~
the number of canonical vectors
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 An S3 object containing the following slots:
  \item{rkcor}{Robsut kernel canonical correlation}
  \item{ rxcoef}{Robsut  kernel canonical  coficient of X dataset}
  \item{ rycoef}{Robsut   kernel canonical  coficient of Y dataset}
  \item{ rxcv}{Robsut   kernel canonical  vector of X dataset}
  \item{ rycv}{Robsut  kernel  canonical  vector of Y dataset}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%% ...
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
See also as  \code{\link{ifcca}}, \code{\link{rkcca}}, \code{\link{ifrkcca}}
}


\examples{

##Dummy data:

X <- matrix(rnorm(1000),100); Y <- matrix(rnorm(1000),100)

rkcca(X,Y, "Huber",  "rbfdot", 1e-05,  10)
}

