\name{ECD}
\alias{ECD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ECD algorithm
%%  ~~function to do ... ~~
}
\description{
Function for solving envelope based on ECD algorithm
}
\usage{
ECD(S, Sx, d, maxiter = 2000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Sx}{Sample covariance matrix}
  \item{S}{Weighted covariance matrix across all clusters}
  \item{d}{Selected envelope dimension}
  \item{maxiter}{Maximum iteration}
}
\details{
Used in clemm_em function}

\value{
\item{Ghat}{Eestimation of envelope basis}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
