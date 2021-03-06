\name{ifcca}
\alias{ifcca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Influence Funciton of Canonical Correlation Analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
##To define the robustness in statistics, different approaches have been pro-
 posed, for example, the minimax approach, the sensitivity curve, the  influence
 function (IF) and the finite sample breakdown point.  Due to its simplic-
ity, the IF is the most useful approach in statistical machine  learning
}
\usage{
ifcca(X, Y, gamma = 1e-05, ncomps = 2, jth = 1)
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
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
the  hyper-parameters
}
  \item{ncomps}{
%%     ~~Describe \code{ncoms} here~~
   the number of canonical vectors
}
  \item{jth}{ the influence function of  the jth canonical vector
%%     ~~Describe \code{jth} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{iflccor}{Influence value of the data by linear canonical correalation}
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

ifcca(X,Y,  1e-05,  2, 2)
}
