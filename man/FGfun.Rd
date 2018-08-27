\name{FGfun}
\alias{FGfun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Objective function and its gradient
%%  ~~function to do ... ~~
}
\description{
Objective function and its gradient for mannifold optimization
}
\usage{
FGfun(W, p1, Sx, Sk)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{W}{Envelope subspace basis}
  \item{p1}{The element p1_{ik} is the probability that i-th obersevation belongs to k-th cluster}
  \item{Sx}{Sample covariance matrix}
  \item{Sk}{Intra-cluster covariance}
}
\details{
Used for manifold optimization algorithm to solve for envelope basis
}
\value{
\item{F}{The value of objective function}
\item{G}{The value of gradient of objective function}
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
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.

%## The function is currently defined as
%
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
