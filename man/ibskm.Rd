\name{ibskm}
\alias{ibskm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Kernel Matrix Using Identity-by-state Kernel
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
For GWASs, a kernel captures the pairwise similarity across a number of SNPs in each gene.
Kernel projects the genotype data from original high dimensional space to a feature space. One of
the more popular kernels used for genomics similarity is the identity-by-state (IBS) kernel (non-
parametric function of the genotypes)
}
\usage{
ibskm(Z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
%%     ~~Describe \code{Z} here~~
 a data matrix
}
}
\details{
For genome-wide association study, a kernel captures the
pairwise similarity across a number of SNPs in each gene. Kernel projects
the genotype data from original high dimensional space to a feature
space. One popular kernel used for genomics similarity is the
identity-by-state (IBS) kernel, The IBS kernel does not need any assumption on the type of genetic interactions.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{ K}{a  Gram/ kernel  matrix}
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also as \code{\link{gkm}}, \code{\link{lkm}}
}


\examples{
##Dummy data:
X <- matrix(rnorm(200),50)
ibskm(X)
}
