\name{gkm}
\alias{gkm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Kernel Matrix Using Guasian Kernel
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Many radial basis function kernels, such as the Gaussian kernel, map X into a infinte dimensional space.  While the Gaussian kernel has a free parameter (bandwidth), it still follows a number of theoretical properties such as boundedness, consistence, universality, robustness etc. It is the most applicable kernel of the  positive definite kernel based methods.
}
\usage{
gkm(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
 a data matrix.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Many radial basis function kernels, such as the Gaussian kernel, map
input sapce into a infinite dimensional space. The Gaussian kernel has a  a number of theoretical properties such as boundedness, consistence,
universality and robustness, etc.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{K}{a  Gram/ kernel  matrix}
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

 {M. Romanazzi (1992),}
{ Influence in canonical correlation analysis,}
{ Psychometrika
 vol 57(2) (1992) 237-259.}
 }



\examples{
##Dummy data:
X<-matrix(rnorm(1000),100)
gkm(X)
}
