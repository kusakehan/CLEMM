\name{mixll_gamma}
\alias{mixll_gamma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{"Complete" loglikelihood
%%  ~~function to do ... ~~
}
\description{
Calculate "complete" loglikelihood approximated by function of Gamma
}
\usage{
mixll_gamma(p1, Gamma, M, U)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p1}{Estimation of mixing weights}
  \item{Gamma}{Envelope subspace basis}
  \item{M}{S_k in CLEMM, S in CLEMM-Shared}
  \item{U}{Sample covariance matrix S_x}
}
\details{
Being used in EM algorithm, see clemm-em and gmm-em
}
\value{
The value of "complete" loglikelihood
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
%##---- Should be DIRECTLY executable %!! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for %the standard data sets.

%## The function is currently defined %as
%function (x)
%{
%  }
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
