\name{snpfmrimth3D}
\alias{snpfmrimth3D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
An example of imaging genetics and epi-genetics data  to  calcualte influential  observations from three view data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
 #A function
}
\usage{
snpfmrimth3D(n = 500, gamma = 1e-05, ncomps = 1, jth=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
the sample size
}
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
the  hyper-parameters
}
  \item{ncomps}{
%%     ~~Describe \code{ncoms} here~~
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
 \item{IFim}{ Influence value of multiple kernel canonical correlation analysis for the  ideal data }
\item{IFcm }{ Influence value of  multiple kernel canonical correlation analysis for the  contaminated data }
%% ...
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


 {Md Ashad Alam, Vince D. Calhoun and Yu-Ping Wang (2018),}
{ Identifying outliers using multiple kernel canonical
correlation analysis with application to imaging genetics,}
{Computational Statistics and Data Analysis, Vol. 125, 70- 85}
 }


\seealso{
See also as \code{\link{rkcca}}, \code{\link{snpfmridata}}, \code{\link{ifrkcca}}
}

\examples{

##Dummy data:

n<-100

snpfmrimth3D(n, 0.00001,  10, 1)
}
