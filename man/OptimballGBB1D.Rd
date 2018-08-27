\name{OptimballGBB1D}
\alias{OptimballGBB1D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{1D algorithm for solving envelope
%%  ~~function to do ... ~~
}
\description{
Function for solving envelope based on 1D algorithm
}
\usage{
OptimballGBB1D(p1, Sx, Sk, d, opts=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p1}{The element p1_{ik} is the probability that i-th obersevation belongs to k-th cluster}
  \item{Sx}{Sample covariance matrix}
  \item{Sk}{Intra-cluster covariance}
  \item{d}{Selected envelope dimension}
  \item{opts}{Parameters for manifold optimization}
}
\details{
Used in clemm_em function
}
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
%\examples{

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
