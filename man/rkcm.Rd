\name{rkcm}
\alias{rkcm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Robsut Kernel Center Matrix
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
# A functioin
}
\usage{
rkcm(X, lossfu = "Huber", kernel = "rbfdot")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
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
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{rkcm}{ a square robust kernel center matrix}
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


 {Md Ashad Alam, Vince D. Calhoun and Yu-Ping Wang (2018),}
{ Identifying outliers using multiple kernel canonical
correlation analysis with application to imaging genetics,}
{Computational Statistics and Data Analysis, Vol. 125, 70- 85}
 }




 \seealso{
See also as  \code{\link{ifcca}}, \code{\link{rkcca}}, \code{\link{ifrkcca}}
}


\examples{

##Dummy data:

X <- matrix(rnorm(2000),200); Y <- matrix(rnorm(2000),200)

rkcm(X, "Huber","rbfdot")
}
