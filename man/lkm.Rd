\name{lkm}
\alias{lkm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Kernel Matrix Using Linear  Kernel
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The linear kernel is  used  by  the underlying Euclidean space to define the similarity measure.
Whenever the dimensionality is high, it may allow for more complexity in the function class than what we could measure and assess otherwise
}
\usage{
lkm(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
a data matrix
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The linear kernel is used by the underlying Euclidean space to define
the similarity measure.  Whenever the dimensionality of the data is high, it
may allow for more complexity in the function class than what we could
measure and assess otherwise.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{K}{a kernel matrix.}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
 }
\author{
%%  ~~who you are~~
Md Ashad Alam
<malam@tulane.edu>
}

\references{
%% ~put references to the literature/web site here ~

{Md. Ashad Alam, Hui-Yi Lin, HOng-Wen Deng, Vince Calhour Yu-Ping Wang (2018),}
{A kernel machine method for detecting higher order interactions in
multimodal datasets: Application to schizophrenia,}
{Journal of Neuroscience Methods, Vol. 309, 161-174.}


{Md Ashad Alam, Kenji Fukumizu and  Yu-Ping Wang (2018),}
{Influence Function and Robust Variant of Kernel Canonical Correlation Analysis,}
{Neurocomputing, Vol. 304 (2018) 12-29.}

 {Md Ashad Alam, Vince D. Calhoun and Yu-Ping Wang (2018),}
{ Identifying outliers using multiple kernel canonical
correlation analysis with application to imaging genetics,}
{Computational Statistics and Data Analysis, Vol. 125, 70- 85}
 }




%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also as \code{\link{gkm}}, \code{\link{ibskm}}

}


\examples{
##Dummy data:

X <- matrix(rnorm(500),100)
lkm(X)
}
