\name{mdbw}
\alias{mdbw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bandwidth  of the Gaussian kernel
}
\description{
 A  median of the pairwise distance of the  data
}
\usage{
mdbw(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  a data matrix
}
}
\details{
While the Gaussian kernel has a free parameter (bandwidth), it still
follows a number of theoretical properties such as boundedness, consistenc, universality, robustness, etc. It is the most applicable one. In a Gaussian RBF kernel, we need to select an appropriate a bandwidth. It is well known that the parameter has a strong influence on the result of kernel methods. For the Gaussian kernel, we can use the median of the pairwise distance as a bandwidth.
}
\value{
 \item{s}{ a  median of the pairwise distance of the   X dataset}
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


{Md. Ashad Alam, Kenji Fukumizu and  Yu-Ping Wang (2018),}
{Influence Function and Robust Variant of Kernel Canonical Correlation Analysis,}
{Neurocomputing, Vol. 304 (2018) 12-29.}


{Md. Ashad Alam and Kenji Fukumizu (2015),}
{Higher-order regularized kernel canonical correlation analysis,}
{International Journal of Pattern Recognition and Artificial Intelligence, Vol. 29(4) 1551005.}


 {Arthu Gretton, Kenji. Fukumizu, C. H. Teo, L. Song, B. Scholkopf and A. Smola (2008),}
{ A Kernel
statistical test of independence,}
{ in Advances in Neural Information Processing Systems,
 Vol. 20  585–592.}
 }


\seealso{
See also as  \code{\link{lkm}}, \code{\link{gkm}}

}

\examples{

##Dummy data:

X <- matrix(rnorm(1000),100)

mdbw(X)
}

