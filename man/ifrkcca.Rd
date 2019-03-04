\name{ifrkcca}
\alias{ifrkcca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Influence Function of Robust Kernel Canonical Analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
##To define the robustness in statistics, different approaches have been pro-
 posed, for example, the minimax approach, the sensitivity curve, the  influence
 function (IF) and the finite sample breakdown point.  Due to its simplic-
ity, the IF is the most useful approach in statistical machine  learning.
}
\usage{
ifrkcca(X, Y, lossfu = "Huber", kernel = "rbfdot", gamma = 0.00001, ncomps = 10, jth = 1)
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
  \item{jth}{
%%     ~~Describe \code{jth} here~~
the influence function of  the jth canonical vector
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{ifrkcor}{Influence value of the data by robust kernel canonical correalation}
\item{ifrkxcv }{Influence value of  cnonical vector of X dataset}
\item{ifrkycv }{Influence value of  cnonical vector of Y dataset}
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
See also as \code{\link{rkcca}}, \code{\link{ifrkcca}}
}


\examples{

##Dummy data:

X <- matrix(rnorm(500),100); Y <- matrix(rnorm(500),100)

ifrkcca(X,Y, lossfu = "Huber", kernel = "rbfdot", gamma = 0.00001, ncomps = 10, jth = 2)
}
