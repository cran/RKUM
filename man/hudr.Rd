\name{hudr}
\alias{hudr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Huber's psi function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The ratio of  the first derivative of the Huber loss fuction to  the argument. Tuning constants is fixed as a meadian vlue.
}
\usage{
hudr(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  vector  values
%%     ~~Describe \code{x} here~~
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~

%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  \item{y}{a real value}
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 See also as \code{\link{hadr}}
}
