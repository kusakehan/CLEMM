\name{Mixturell}
\alias{Mixturell}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Observed data loglikelihood
}
\description{
Calculate the Observed data loglikelihood
}
\usage{
Mixturell(K, da, pi0, mu, Sigma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{da}{Data for clustering}
  \item{K}{Selected number of clusters}
  \item{pi0}{Mixing weights}
  \item{mu}{Clusters centers}
  \item{Sigma}{Clusters covariance}
%%     ~~Describe \code{x} here~~
}
\details{
Used in EM algorithm for checking covergence
}
\value{
The observed loglikelihood
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
