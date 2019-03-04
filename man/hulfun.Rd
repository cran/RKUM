\name{hulfun}
\alias{hulfun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
 A Huber loss function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Tuning constants  of the  Huber loss fuction are fixed in different quintiles of the arguments.
}
\usage{
hulfun(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
a vector values
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Tuning constants  of the  Huber  fuction is fixed as a   median.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{}{a real number}
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
See also as \code{\link{halfun}}

}

