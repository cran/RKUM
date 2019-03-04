\name{snpfmridata}
\alias{snpfmridata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
An example of imaging genetics data to calcualte influential  observations from two view data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
#A function
}
\usage{
snpfmridata(n = 300, gamma=0.00001, ncomps = 2, jth = 1)
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
\item{IFCCAID }{ Influence value of canonical correlation analysis for the  ideal data }
\item{IFCCACD }{ Influence value of canonical correlation analysis for the  contaminated data }
 \item{IFKCCAID }{ Influence value of kernel canonical correlation analysis for the  ideal data }
\item{IFKCCACD }{ Influence value of  kernel canonical correlation analysis for the  contaminated data }
\item{IFHACCAID }{ Influence value of  robsut (Hampel's loss) canonical correlation analysis for the  ideal data }
\item{IFHACCACD }{ Influence value of robsut (Hampel's loss) canonical correlation analysis for the  contaminated data }
\item{IFHUCCAID }{ Influence value of  robsut (Huber's loss) canonical correlation analysis for the  ideal data }
\item{IFHUCCACD }{ Influence value of robsut (Huber's loss) canonical correlation analysis for the  contaminated data }
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



%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also as  \code{\link{rkcca}}, \code{\link{ifrkcca}}, \code{\link{snpfmrimth3D}}
}

\examples{

##Dummy data:

n<-100

snpfmridata(n, 0.00001,  10, jth = 1)
}
